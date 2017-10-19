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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamRecordCodec, SamSource}
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sorter}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker
import htsjdk.samtools.util.{CloserUtil, CoordMath, SequenceUtil}

import scala.collection.mutable.ListBuffer
import scala.math.{max, min}

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

  /** Returns an iterator over all records that are not the primary r1 and r2. */
  def allSupplementaryAndSecondary: Iterator[SamRecord] =
    r1Supplementals.iterator ++ r2Supplementals.iterator ++ r1Secondaries.iterator ++ r2Secondaries.iterator

  /** Returns an iterator of all reads for the template. */
  def allReads: Iterator[SamRecord] = r1.iterator ++ r2.iterator ++ allSupplementaryAndSecondary

  /** The total count of records for the template. */
  lazy val readCount: Int = r1.size + r2.size + r1Supplementals.size + r2Supplementals.size + r1Secondaries.size + r2Secondaries.size
}

object Template {
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
    val r1Supp = ListBuffer[SamRecord]()
    val r2Supp = ListBuffer[SamRecord]()
    val r1Sec  = ListBuffer[SamRecord]()
    val r2Sec  = ListBuffer[SamRecord]()

    var name: Option[String] = None

    recs foreach { r =>
      if (name.isEmpty) name = Some(r.name)
      else require(name.get == r.name, "Reads for same template have different query names.")

      val isR1 = !r.paired || r.firstOfPair
      if (isR1) {
        if      (r.secondary)     r1Sec += r
        else if (r.supplementary) r1Supp += r
        else if (r1.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R1s for ${r.name}")
        else r1 = Some(r)
      }
      else {
        if      (r.secondary)     r2Sec += r
        else if (r.supplementary) r2Supp += r
        else if (r2.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R2s for ${r.name}")
        else r2 = Some(r)
      }
    }

    Template(r1=r1, r2=r2, r1Supplementals=r1Supp.toList, r2Supplementals=r2Supp.toList, r1Secondaries=r1Sec.toList, r2Secondaries=r2Sec.toList)
  }
}

/**
  * Utility methods for working with BAMs.
  */
object Bams extends LazyLogging {
  /** The default maximum # of records to keep and sort in memory. */
  val MaxInMemory: Int = 1e6.toInt

  /** Generates a [[Sorter]] for doing disk-backed sorting of objects. */
  def sorter(order: SamOrder,
             header: SAMFileHeader,
             maxRecordsInRam: Int = MaxInMemory,
             tmpDir: DirPath = Io.tmpDir): Sorter[SamRecord,order.A] = {
    // FIXME: tmpDir not used
    new Sorter(maxRecordsInRam, new SamRecordCodec(header), order.sortkey)
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
    }
    else {
      logger.info("Sorting into queryname order.")
      val progress = ProgressLogger(this.logger, "Queryname sorted")
      val sort     = sorter(SamOrder.Queryname, header, maxInMemory, tmpDir)
      iterator.foreach { rec =>
        sort += rec
        progress.record(rec)
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
                       tmpDir: DirPath = Io.tmpDir): Iterator[Template] = {
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
                       tmpDir: DirPath): Iterator[Template] = {
    val queryIterator = queryGroupedIterator(iterator, header, maxInMemory, tmpDir)

    new Iterator[Template] {
      override def hasNext: Boolean = queryIterator.hasNext
      override def next(): Template = {
        require(hasNext, "next() called on empty iterator")
        Template(queryIterator)
      }
    }
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
    // FIXME: tmpDir not used
    val sort = new Sorter(maxInMemory, new SamRecordCodec(header), f)
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
    * @param ref a reference sequence file walker to pull the reference information from
    */
  def regenerateNmUqMdTags(rec: SamRecord, ref: ReferenceSequenceFileWalker): Unit = {
    import SAMTag._
    if (rec.unmapped) {
      rec(NM.name()) =  null
      rec(UQ.name()) =  null
      rec(MD.name()) =  null
    }
    else {
      val refBases = ref.get(rec.refIndex).getBases
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
}
