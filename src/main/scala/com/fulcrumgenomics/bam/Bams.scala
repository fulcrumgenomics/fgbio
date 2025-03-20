/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.SamOrder.Queryname
import com.fulcrumgenomics.bam.api._
import com.fulcrumgenomics.commons.collection.{BetterBufferedIterator, SelfClosingIterator}
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sorter}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequence
import htsjdk.samtools.util.{CloserUtil, CoordMath, Murmur3, SequenceUtil}

import java.io.Closeable
import scala.math.{max, min}


/** Class to store information about an alignment, as described in the SAM SA tag. */
private[bam] case class AlignmentInfo(refIndex: Int, start: Int, positiveStrand: Boolean, cigar: Option[Cigar], mapq: Int, nm: Int) {
  def negativeStrand: Boolean = !positiveStrand
  def refName(header: SAMFileHeader): String = header.getSequence(refIndex).getSequenceName
  def mapped: Boolean = cigar.isDefined
  def unmapped: Boolean = !mapped
  private def _cigar: Cigar = cigar.getOrElse {
    throw new IllegalStateException(s"Cannot get cigar for AlignmentInfo: ${this}")
  }
  def end: Int = start + _cigar.lengthOnTarget - 1
  def unclippedStart: Int = _cigar.unclippedStart(start)
  def unclippedEnd: Int = _cigar.unclippedEnd(end)
  /** Returns a formatted alignment as per the SA tag: `(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+` */
  def toSA(header: SAMFileHeader): String = {
    val strand  = if (positiveStrand) '+' else '-'
    val cigar   = this.cigar.getOrElse("*")
    f"${refName(header)},${start},${strand},${cigar},${mapq},${nm}"
  }
}

private[bam] object AlignmentInfo {
  def apply(rec: SamRecord, mate: Boolean = false): AlignmentInfo = {
    if (mate) {
      val mateRefIndex = if (rec.unpaired || rec.mateUnmapped) Int.MaxValue else rec.mateRefIndex
      val mateCigar    = if (rec.unpaired || rec.mateUnmapped) None else Some(rec.mateCigar.getOrElse {
        throw new IllegalStateException(s"Mate CIGAR (Tag 'MC') not found for $rec, consider using SetMateInformation.")
      })
      // NB: mateCigar has already checked for the existence of the MC tag, so using .get here is fine
      val mateStart    = if (rec.unpaired || rec.mateUnmapped) Int.MaxValue else rec.mateStart
      val mateStrand   = if (rec.unpaired || rec.mateUnmapped) true else rec.matePositiveStrand
      AlignmentInfo(mateRefIndex, mateStart, mateStrand, mateCigar, rec.mapq, 0)
    } else {
      val refIndex       = if (rec.unmapped)     Int.MaxValue else rec.refIndex
      val positiveStrand = rec.positiveStrand
      val start          = if (rec.unmapped)     Int.MaxValue else rec.start
      AlignmentInfo(refIndex, start, positiveStrand, Some(rec.cigar), rec.mapq, rec.getOrElse(SAMTag.NM.name(), 0))
    }
  }

  def apply(sa: String, header: SAMFileHeader): AlignmentInfo = {
    val parts    = sa.split(",")
    require(parts.length == 6, f"Could not parse SA tag: ${sa}")
    val refIndex = header.getSequenceIndex(parts(0))
    AlignmentInfo(refIndex, parts(1).toInt, parts(2) == "+", Some(Cigar(parts(3))), parts(4).toInt, parts(5).toInt)
  }
  /** Returns a formatted alignment as per the SA tag: `(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+` */
  def toSA(rec: SamRecord): String = AlignmentInfo(rec).toSA(rec.header)
}

/**
  * Class that represents all reads from a template within a BAM file.
  */
case class Template(r1: Option[SamRecord],
                    r2: Option[SamRecord],
                    r1Supplementals: Seq[SamRecord] = Nil,
                    r2Supplementals: Seq[SamRecord] = Nil,
                    r1Secondaries:   Seq[SamRecord] = Nil,
                    r2Secondaries:   Seq[SamRecord] = Nil) {

  /** Gets the query name of the template (assumes all records have the same read name). */
  val name: String = {
    if (r1.isDefined) r1.get.name
    else if (r2.isDefined) r2.get.name
    else Seq(r1Supplementals, r2Supplementals, r1Secondaries, r2Secondaries).find(_.nonEmpty).map(_.head.name)
      .getOrElse(throw new IllegalStateException("Template created with no reads!"))
  }

  /** Returns an iterator over all non-secondary non-supplementary reads. */
  def primaryReads: Iterator[SamRecord] = r1.iterator ++ r2.iterator

  /** If both the R1 and R2 are present and mapped to the same chromosome returns the pair orientation, else None. */
  def pairOrientation: Option[PairOrientation] = (r1, r2) match {
    case (Some(a), Some(b)) if a.mapped && b.mapped && a.refIndex == b.refIndex => Some(SamPairUtil.getPairOrientation(a.asSam))
    case _ => None
  }

  /** Returns an iterator over all records that are not the primary r1 and r2. */
  def allSupplementaryAndSecondary: Iterator[SamRecord] =
    r1Supplementals.iterator ++ r2Supplementals.iterator ++ r1Secondaries.iterator ++ r2Secondaries.iterator

  /** Returns an iterator of all reads for the template. */
  def allReads: Iterator[SamRecord] = r1.iterator ++ r2.iterator ++ allSupplementaryAndSecondary

  /** Returns an iterator of all read 1s for the template. */
  def allR1s: Iterator[SamRecord] = r1.iterator ++ r1Secondaries.iterator ++ r1Supplementals.iterator

  /** Returns an iterator of all read 2s for the template. */
  def allR2s: Iterator[SamRecord] = r2.iterator ++ r2Secondaries.iterator ++ r2Supplementals.iterator
  /**
    * Produces a copy of the template that has had mapping information removed.  Discards secondary and supplementary
    * records, and retains the primary records after un-mapping them.  Auxillary tags listed in [[Bams.AlignmentTags]]
    * are removed from all records.  Lastly the `OA` tag is added to reads to record the Original Alignment.
    */
  def unmapped: Template = {
    val x1 = r1.map(_.clone())
    val x2 = r2.map(_.clone())
    Seq(x1, x2).flatten.foreach { r =>
      if (r.mapped) {
        val nm = r.get[Int]("NM").getOrElse("")
        r("OA") = s"${r.refName},${r.start},${if (r.positiveStrand) "+" else "-"},${r.cigar},${r.mapq},$nm"
      }
      SAMUtils.makeReadUnmapped(r.asSam)
      SamRecordClipper.TagsToInvalidate.foreach(t => r(t) = null)
    }

    (x1, x2) match {
      case (Some(a), Some(b)) => SamPairUtil.setMateInfo(a.asSam, b.asSam)
      case _ => ()
    }

    Template(x1, x2)
  }


  /** Fixes mate information and sets mate cigar and mate score on all primary, secondary, and supplementary records. */
  def fixMateInfo(): Unit = {
    // Developer note: the mate score ("ms") tag is used by samtools markdup
    // Set all mate info on BOTH secondary and supplementary records, not just supplementary records.  We also need to
    // add the "pa" and "pm" tags with information about the primary alignments.  Finally, we need the MQ tag!
    val r1NonPrimary = r1Supplementals ++ r1Secondaries
    val r2NonPrimary = r2Supplementals ++ r2Secondaries
    for (primary <- r1; nonPrimary <- r2NonPrimary) {
      SamPairUtil.setMateInformationOnSupplementalAlignment(nonPrimary.asSam, primary.asSam, true)
      primary.get[Int]("AS").foreach(nonPrimary("ms") = _)
      nonPrimary(SAMTag.MQ.name()) = primary.mapq
      nonPrimary(Template.MatePrimarySamTag) = AlignmentInfo.toSA(primary)
      r2.foreach(r => nonPrimary(Template.ReadPrimarySamTag) = AlignmentInfo.toSA(r))
    }
    for (primary <- r2; nonPrimary <- r1NonPrimary) {
      SamPairUtil.setMateInformationOnSupplementalAlignment(nonPrimary.asSam, primary.asSam, true)
      primary.get[Int]("AS").foreach(nonPrimary("ms") = _)
      nonPrimary(SAMTag.MQ.name()) = primary.mapq
      nonPrimary(Template.MatePrimarySamTag) = AlignmentInfo.toSA(primary)
      r1.foreach(r => nonPrimary(Template.ReadPrimarySamTag) = AlignmentInfo.toSA(r))
    }
    for (first <- r1; second <- r2) {
      SamPairUtil.setMateInfo(first.asSam, second.asSam, true)
      first.get[Int]("AS").foreach(second("ms") = _)
      second.get[Int]("AS").foreach(first("ms") = _)
    }
  }

  /** The total count of records for the template. */
  lazy val readCount: Int = r1.size + r2.size + r1Supplementals.size + r2Supplementals.size + r1Secondaries.size + r2Secondaries.size
}

object Template {
  /** The local SAM tag to store the alignment information of the primary alignment (in the same format as the SA tag) */
  val ReadPrimarySamTag: String = "rp"
  /** The local SAM tag to store the alignment information of the mate's primary alignment (in the same format as the SA tag) */
  val MatePrimarySamTag: String = "mp"
  /**
    * Generates a Template for the next template in the buffered iterator. Assumes that the
    * iterator is queryname sorted or grouped.
    */
  def apply(recs: BetterBufferedIterator[SamRecord]): Template = {
    require(recs.hasNext)
    val name = recs.head.name
    apply(recs.takeWhile(_.name == name))
  }

  /** Generates a Template object from an iterator of SamRecords from the same template. */
  def apply(recs: Iterator[SamRecord]): Template = {
    var r1: Option[SamRecord] = None
    var r2: Option[SamRecord] = None
    var r1Supp: List[SamRecord] = Nil
    var r2Supp: List[SamRecord] = Nil
    var r1Sec : List[SamRecord] = Nil
    var r2Sec : List[SamRecord] = Nil

    var name: Option[String] = None

    recs foreach { r =>
      if (name.isEmpty) name = Some(r.name)
      else require(name.get == r.name, "Reads for same template have different query names.")

      val isR1 = !r.paired || r.firstOfPair
      if (isR1) {
        if      (r.secondary)     r1Sec  = r :: r1Sec
        else if (r.supplementary) r1Supp = r :: r1Supp
        else if (r1.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R1s for ${r.name}")
        else r1 = Some(r)
      }
      else {
        if      (r.secondary)     r2Sec  = r :: r2Sec
        else if (r.supplementary) r2Supp = r :: r2Supp
        else if (r2.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R2s for ${r.name}")
        else r2 = Some(r)
      }
    }

    Template(r1=r1, r2=r2, r1Supplementals=r1Supp, r2Supplementals=r2Supp, r1Secondaries=r1Sec, r2Secondaries=r2Sec)
  }
}

/**
  * Utility methods for working with BAMs.
  */
object Bams extends LazyLogging {
  /** The default maximum # of records to keep and sort in memory. */
  val MaxInMemory: Int = 1e6.toInt

  /** Auxillary tags that should be cleared or re-calculated when unmapping or changing the alignment of a SamRecord. */
  val AlignmentTags: Seq[String] = Seq("MD", "NM", "UQ")

  /** Generates a [[com.fulcrumgenomics.util.Sorter]] for doing disk-backed sorting of objects. */
  def sorter(order: SamOrder,
             header: SAMFileHeader,
             maxRecordsInRam: Int = MaxInMemory,
             tmpDir: DirPath = Io.tmpDir): Sorter[SamRecord,order.A] = {
    new Sorter(maxRecordsInRam, new SamRecordCodec(header), order.sortkey, tmpDir=tmpDir)
  }

  /** A wrapper to order objects of type [[TagType]] using the ordering given.  Used when sorting by tag where we wish
    * to sort on the transformed tag value. */
  private case class SortByTagKey[TagType](value: TagType)(implicit ordering: Ordering[TagType]) extends Ordered[SortByTagKey[TagType]] {
    override def compare(that: SortByTagKey[TagType]): Int = ordering.compare(this.value, that.value)
  }

  /**
    * Returns an iterator over the records in the given reader in such a way that all
    * reads with the same query name are adjacent in the iterator. Does NOT guarantee
    * a queryname sort, merely a query grouping.
    *
    * @param in a SamReader from which to consume records
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir an optional temp directory to use for temporary sorting files if needed
    * @return an Iterator with reads from the same query grouped together
    */
  def queryGroupedIterator(in: SamSource,
                           maxInMemory: Int = MaxInMemory,
                           tmpDir: DirPath = Io.tmpDir): BetterBufferedIterator[SamRecord] = {
    queryGroupedIterator(in.iterator, in.header, maxInMemory, tmpDir)
  }

  /**
    * Returns an iterator over the records in the given reader in such a way that all
    * reads with the same query name are adjacent in the iterator. Does NOT guarantee
    * a queryname sort, merely a query grouping.
    *
    * @param iterator an iterator from which to consume records
    * @param header the header associated with the records.
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir an optional temp directory to use for temporary sorting files if needed
    * @return an Iterator with reads from the same query grouped together
    */
  def queryGroupedIterator(iterator: Iterator[SamRecord],
                           header: SAMFileHeader,
                           maxInMemory: Int,
                           tmpDir: DirPath): SelfClosingIterator[SamRecord] = {
    if (header.getSortOrder == SortOrder.queryname || header.getGroupOrder == GroupOrder.query) {
      iterator match {
        case iter: SelfClosingIterator[SamRecord] => iter
        case _ => new SelfClosingIterator(iterator.bufferBetter, () => CloserUtil.close(iterator))
      }
    } else {
      querySortedIterator(iterator, header, maxInMemory, tmpDir)
    }
  }

  /** Returns an iterator over records in such a way that all reads with the same query name are adjacent in the
    * iterator. Although a queryname sort is guaranteed, the sort order may not be consistent with other queryname
    * sorting implementations, especially in other tool kits.
    *
    * @param in a SamReader from which to consume records
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    * @return an Iterator with reads from the same query grouped together
    */
  def querySortedIterator(in: SamSource,
                          maxInMemory: Int = MaxInMemory,
                          tmpDir: DirPath = Io.tmpDir): BetterBufferedIterator[SamRecord] = {
    querySortedIterator(in.iterator, in.header, maxInMemory, tmpDir)
  }

  /** Returns an iterator over records in such a way that all reads with the same query name are adjacent in the
    * iterator. Although a queryname sort is guaranteed, the sort order may not be consistent with other queryname
    * sorting implementations, especially in other tool kits.
    *
    * @param iterator an iterator from which to consume records
    * @param header the header associated with the records
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    * @return an Iterator with reads from the same query grouped together
    */
  def querySortedIterator(iterator: Iterator[SamRecord],
                          header: SAMFileHeader,
                          maxInMemory: Int,
                          tmpDir: DirPath): SelfClosingIterator[SamRecord] = {
    (SamOrder(header), iterator) match {
      case (Some(Queryname), _iterator: SelfClosingIterator[SamRecord]) => _iterator
      case (Some(Queryname), _) => new SelfClosingIterator(iterator.bufferBetter, () => CloserUtil.close(iterator))
      case (_, _) =>
        logger.info(parts = "Sorting into queryname order.")
        val progress = ProgressLogger(this.logger, "records", "sorted")
        val sort     = sorter(Queryname, header, maxInMemory, tmpDir)
        iterator.foreach { rec =>
          progress.record(rec)
          sort.write(rec)
        }
        new SelfClosingIterator(sort.iterator, () => sort.close())
    }
  }

  /**
    * Returns an iterator of Template objects generated from all the reads in the SamReader.
    * If the SamReader is not query sorted or grouped, a sort to queryname will be performed.
    *
    * @param in a SamReader from which to consume records
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir an optional temp directory to use for temporary sorting files if needed
    * @return an Iterator of Template objects
    */
  def templateIterator(in: SamSource,
                       maxInMemory: Int = MaxInMemory,
                       tmpDir: DirPath = Io.tmpDir): SelfClosingIterator[Template] = {
    templateIterator(in.iterator, in.header, maxInMemory, tmpDir)
  }

  /**
    * Returns an iterator of Template objects generated from all the reads in the SamReader.
    * If the SamReader is not query sorted or grouped, a sort to queryname will be performed.
    *
    * @param iterator an iterator from which to consume records
    * @param header the header associated with the records.
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir an optional temp directory to use for temporary sorting files if needed
    * @return an Iterator of Template objects
    */
  def templateIterator(iterator: Iterator[SamRecord],
                       header: SAMFileHeader,
                       maxInMemory: Int,
                       tmpDir: DirPath): SelfClosingIterator[Template] = {
    val queryIterator = queryGroupedIterator(iterator, header, maxInMemory, tmpDir)

    val iter = new Iterator[Template] {
      override def hasNext: Boolean = queryIterator.hasNext
      override def next(): Template = {
        require(hasNext, "next() called on empty iterator")
        Template(queryIterator)
      }
    }

    new SelfClosingIterator(iter, () => queryIterator.close())
  }

  /** Return an iterator over records sorted and grouped into [[Template]] objects. Although a queryname sort is
    * guaranteed, the sort order may not be consistent with other queryname sorting implementations, especially in other
    * tool kits. See [[templateIterator]] for a [[Template]] iterator which emits templates in a non-guaranteed sort
    * order.
    *
    * @see [[templateIterator]]
    *
    * @param in a SamReader from which to consume records
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    * @return an Iterator of queryname sorted Template objects
    */
  def templateSortedIterator(in: SamSource,
                             maxInMemory: Int = MaxInMemory,
                             tmpDir: DirPath = Io.tmpDir): SelfClosingIterator[Template] = {
    templateSortedIterator(in.iterator, in.header, maxInMemory, tmpDir)
  }

  /** Return an iterator over records sorted and grouped into [[Template]] objects. Although a queryname sort is
    * guaranteed, the sort order may not be consistent with other queryname sorting implementations, especially in other
    * tool kits. See [[templateIterator]] for a [[Template]] iterator which emits templates in a non-guaranteed sort
    * order.
    *
    * @see [[templateIterator]]
    *
    * @param iterator an iterator from which to consume records
    * @param header the header associated with the records
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    * @return an Iterator of queryname sorted Template objects
    */
  def templateSortedIterator(iterator: Iterator[SamRecord],
                             header: SAMFileHeader,
                             maxInMemory: Int,
                             tmpDir: DirPath): SelfClosingIterator[Template] = {
    val queryIterator = querySortedIterator(iterator, header, maxInMemory, tmpDir)

    val _iterator = new Iterator[Template] {
      override def hasNext: Boolean = queryIterator.hasNext
      override def next(): Template   = {
        require(hasNext, "next() called on empty iterator")
        Template(queryIterator)
      }
    }

    new SelfClosingIterator(_iterator, () => queryIterator.close())
  }

  /** Returns an iterator over records, grouped into [[Template]] objects, in a random order. Reads are first
    * sorted based on a hash of their read name in order to randomize the order while collating by read name.
    * The randomized reads are then iterated and grouped into [[Template]] objects.
    *
    * @param in the [[SamSource]] from which to pull reads and the SAM header
    * @param randomSeed a seed used to generate the hashes that drive the random sorting.  Different seeds will produce
    *                   different randomizations of the reads.  Using the same seed will result in the same ordering
    *                   of the data across invocations.
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    */
  def templateRandomIterator(in: SamSource,
                             randomSeed: Int = 42,
                             maxInMemory: Int = Bams.MaxInMemory,
                             tmpDir: DirPath = Io.tmpDir): Iterator[Template] = {
    val hasher = new Murmur3(randomSeed)
    val sorter = {
      val f = (r: SamRecord) => SamOrder.RandomQueryKey(hasher.hashUnencodedChars(r.name), r.name, r.asSam.getFlags)
      new Sorter(maxInMemory, new SamRecordCodec(in.header), f)
    }

    val sortProgress = ProgressLogger(logger, verb = "sorted", unit = 5e6.toInt)
    in.foreach { rec =>
      sorter += rec
      sortProgress.record()
    }

    val header = in.header.clone()
    SamOrder.RandomQuery.applyTo(header)
    Bams.templateIterator(sorter.iterator, header, maxInMemory, tmpDir)
  }


  /** Returns an iterator over the records in the given iterator such that the order of the records returned is
    * determined by the value of the given SAM tag, which can optionally be transformed.
    *
    * @param iterator an iterator from which to consume records
    * @param header the header to use for the sorted records
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir the temporary directory to use when spilling to disk
    * @param tag the SAM tag (two-letter key) to sort by
    * @param defaultValue the default value, if any, otherwise require all records to have the SAM tag.
    * @param transform the transform to apply to the value of the SAM tags (default is the identity)
    * @param ordering the ordering of [[B]].
    * @tparam A the type of the SAM tag
    * @tparam B the type of the SAM tag after any transformation, or just [[A]] if no transform is given.
    * @return an Iterator over records sorted by the given SAM tag, optionally transformed.
    */
  def sortByTransformedTag[A, B](iterator: Iterator[SamRecord],
                                 header: SAMFileHeader,
                                 maxInMemory: Int = MaxInMemory,
                                 tmpDir: DirPath = Io.tmpDir,
                                 tag: String,
                                 defaultValue: Option[A] = None,
                                 transform: A => B)
                                (implicit ordering: Ordering[B]): Iterator[SamRecord] = {
    val f: SamRecord => SortByTagKey[B] = r => SortByTagKey(
      transform(
        r.get[A](tag).orElse(defaultValue).getOrElse {
          throw new IllegalStateException(s"Missing value for tag '$tag' in record: $r")
        }
      )
    )
    val sort = new Sorter(maxInMemory, new SamRecordCodec(header), f, tmpDir=tmpDir)
    sort ++= iterator
    new SelfClosingIterator(sort.iterator, () => sort.close())
  }

  /** Returns an iterator over the records in the given iterator such that the order of the records returned is
    * determined by the value of the given SAM tag.
    *
    * @param iterator an iterator from which to consume records
    * @param header the header to use for the sorted records
    * @param maxInMemory the maximum number of records to keep and sort in memory if sorting is needed
    * @param tmpDir the temporary directory to use when spilling to disk
    * @param tag the SAM tag (two-letter key) to sort by
    * @param defaultValue the default value, if any, otherwise require all records to have the SAM tag.
    * @param ordering the ordering of [[A]].
    * @tparam A the type of the SAM tag
    * @return an Iterator over records sorted by the given SAM tag.
    */
  def sortByTag[A](iterator: Iterator[SamRecord],
                   header: SAMFileHeader,
                   maxInMemory: Int = MaxInMemory,
                   tmpDir: DirPath = Io.tmpDir,
                   tag: String,
                   defaultValue: Option[A] = None)
                  (implicit ordering: Ordering[A]): Iterator[SamRecord] = {
    this.sortByTransformedTag[A, A](iterator=iterator, header=header, maxInMemory=maxInMemory, tmpDir=tmpDir,
      tag=tag, defaultValue=defaultValue, transform = a => a)
  }

  /**
    * Ensures that any NM/UQ/MD tags on the read are accurate.  If the read is unmapped, any existing
    * values are removed.  If the read is mapped all three tags will have values regenerated.
    *
    * @param rec the SamRecord to update
    * @param ref a reference sequence if the record is mapped
    */
  def regenerateNmUqMdTags(rec: SamRecord, ref: => ReferenceSequence): Unit = {
    import SAMTag._
    if (rec.unmapped) {
      rec(NM.name()) =  null
      rec(UQ.name()) =  null
      rec(MD.name()) =  null
    }
    else {
      val refBases = ref.getBases
      SequenceUtil.calculateMdAndNmTags(rec.asSam, refBases, true, true)
      if (rec.quals != null && rec.quals.length != 0) {
        rec(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(rec.asSam, refBases, 0)
      }
    }
  }

  /**
    * Calculates the coordinates of the insert represented by this record and returns them
    * as a pair of 1-based closed ended coordinates.
    *
    * Invalid to call on a read that is not mapped in a pair to the same chromosome as it's mate.
    *
    * @return the start and end position of the insert as a tuple, always with start <= end
    */
  def insertCoordinates(rec: SamRecord): (Int, Int) = {
    require(rec.paired, s"Read ${rec.name} is not paired, cannot calculate insert coordinates.")
    require(rec.mapped, s"Read ${rec.name} must be mapped to calculate insert coordinates.")
    require(rec.mateMapped, s"Read ${rec.name} must have mapped mate to calculate insert coordinates.")
    require(rec.refIndex == rec.mateRefIndex, s"Read ${rec.name} must be mapped to same chrom as it's mate.")

    val isize      = rec.insertSize
    val firstEnd   = if (rec.negativeStrand) rec.end else rec.start
    val adjustment = if (isize < 0) 1 else -1
    val secondEnd  = firstEnd + isize + adjustment

    (min(firstEnd,secondEnd), max(firstEnd,secondEnd))
  }

  /** If the read is mapped in an FR pair, returns the distance of the position from the other end
    * of the template, other wise returns None.
    *
    * @param rec the SamRecord whose insert to calculate the position within
    * @param genomicPosition the genomic position of interest (NOT the position within the read)
    */
  def positionFromOtherEndOfTemplate(rec: SamRecord, genomicPosition: Int): Option[Int] = {
    if (rec.isFrPair && rec.insertSize != 0) {
      val isize        = rec.insertSize
      val thisEnd      = if (rec.negativeStrand) rec.end else rec.start
      val adjustment   = if (isize < 0) 1 else -1
      val otherEnd     = thisEnd + isize + adjustment
      val (start, end) = (min(thisEnd, otherEnd), max(thisEnd,otherEnd))
      require(genomicPosition >= start && genomicPosition <= end, s"genomicPosition is outside of template for $rec")

      if (isize < 0) Some(CoordMath.getLength(otherEnd, genomicPosition))
      else           Some(CoordMath.getLength(genomicPosition, otherEnd))
    }
    else {
      None
    }
  }

  /** Returns true if the header is queryname sorted or query grouped */
  def isQueryGrouped(header: SAMFileHeader): Boolean = {
    header.getSortOrder == SortOrder.queryname || header.getGroupOrder == GroupOrder.query
  }

  /** Requires that the header is queryname sorted or query grouped. */
  def requireQueryGrouped(header: SAMFileHeader, toolName: String = "<tool>", path: Option[PathToBam] = None): Unit = {
    require(isQueryGrouped(header),
      "Input was not queryname sorted or query grouped, found: " +
        s"SO:${header.getSortOrder} GO:${header.getGroupOrder}" +
        Option(header.getAttribute("SS")).map(ss => f" SS:$ss").getOrElse("") +
        f". Use `samtools sort -n -u in.bam | fgbio $toolName -i /dev/stdin`" +
        path.map(p => f"\nPath: $p").getOrElse("")
    )
  }

  /** Builds a [[Writer]] of [[SamRecord]]s that regenerates the NM, UQ, and MD tags using the given map of reference
    * sequences.
    *
    * @param writer the writer to write to
    * @param ref the path to the reference FASTA
    */
  def nmUqMdTagRegeneratingWriter(writer: SamWriter, ref: PathToFasta): Writer[SamRecord] with Closeable = {
    logger.debug("Reading the reference fasta into memory")
    val refMap = ReferenceSequenceIterator(ref, stripComments=true).map { ref => ref.getContigIndex -> ref}.toMap
    logger.debug(f"Read ${refMap.size}%,d contigs.")

    // Create the final writer based on if the full reference has been loaded, or not
    new Writer[SamRecord] with Closeable {
      override def write(rec: SamRecord): Unit = {
        Bams.regenerateNmUqMdTags(rec, refMap(rec.refIndex))
        writer += rec
      }
      def close(): Unit = writer.close()
    }
  }
}
