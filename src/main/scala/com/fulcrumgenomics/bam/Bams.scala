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
import com.fulcrumgenomics.util.{BetterBufferedIterator, ProgressLogger}
import dagr.commons.util.LazyLogging
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.{SequenceUtil, SortingCollection}
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker

import scala.collection.JavaConversions.iterableAsScalaIterable
import scala.collection.mutable.ListBuffer

/**
  * Class that represents all reads from a template within a BAM file.
  */
case class Template(r1: Option[SAMRecord],
                    r2: Option[SAMRecord],
                    r1Supplementals: Seq[SAMRecord] = Nil,
                    r2Supplementals: Seq[SAMRecord] = Nil,
                    r1Secondaries:   Seq[SAMRecord] = Nil,
                    r2Secondaries:   Seq[SAMRecord] = Nil) {

  /** Gets the query name of the template (assumes all records have the same read name). */
  val name: String = {
    if (r1.isDefined) r1.get.getReadName
    else if (r2.isDefined) r2.get.getReadName
    else Seq(r1Supplementals, r2Supplementals, r1Secondaries, r2Secondaries).find(_.nonEmpty).map(_.head.getReadName)
      .getOrElse(throw new IllegalStateException("Template created with no reads!"))
  }

  /** Returns an iterator over all records that are not the primary r1 and r2. */
  def allSupplementaryAndSecondary: Iterator[SAMRecord] =
    r1Supplementals.iterator ++ r2Supplementals.iterator ++ r1Secondaries.iterator ++ r2Secondaries.iterator

  /** Returns an iterator of all reads for the template. */
  def allReads: Iterator[SAMRecord] = r1.iterator ++ r2.iterator ++ allSupplementaryAndSecondary

  /** The total count of records for the template. */
  lazy val readCount: Int = r1.size + r2.size + r1Supplementals.size + r2Supplementals.size + r1Secondaries.size + r2Secondaries.size
}

object Template {
  /**
    * Generates a Template for the next template in the buffered iterator. Assumes that the
    * iterator is queryname sorted or grouped.
    */
  def apply(recs: BetterBufferedIterator[SAMRecord]): Template = {
    require(recs.hasNext)
    val name = recs.head.getReadName
    apply(recs.takeWhile(_.getReadName == name))
  }

  /** Generates a Template object from an iterator of SAMRecords from the same template. */
  def apply(recs: Iterator[SAMRecord]): Template = {
    var r1: Option[SAMRecord] = None
    var r2: Option[SAMRecord] = None
    val r1Supp = ListBuffer[SAMRecord]()
    val r2Supp = ListBuffer[SAMRecord]()
    val r1Sec  = ListBuffer[SAMRecord]()
    val r2Sec  = ListBuffer[SAMRecord]()

    var name: Option[String] = None

    recs foreach { r =>
      if (name.isEmpty) name = Some(r.getReadName)
      else require(name.get == r.getReadName, "Reads for same template have different query names.")

      val isR1 = !r.getReadPairedFlag || r.getFirstOfPairFlag
      if (isR1) {
        if      (r.getNotPrimaryAlignmentFlag)    r1Sec += r
        else if (r.getSupplementaryAlignmentFlag) r1Supp += r
        else if (r1.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R1s for ${r.getReadName}")
        else r1 = Some(r)
      }
      else {
        if      (r.getNotPrimaryAlignmentFlag)    r2Sec += r
        else if (r.getSupplementaryAlignmentFlag) r2Supp += r
        else if (r2.isDefined) throw new IllegalArgumentException(s"Multiple non-secondary, non-supplemental R2s for ${r.getReadName}")
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
  /**
    * Builds a sorting collection for sorting SAMRecords.
    *
    * @param order the desired output sort order
    * @param header the header to use to encode records during sorting
    * @param maxInMemory the maximum number of records to keep and sort in memory
    * @param tmpDir an optional temp directory to use for temporary sorting files
    * @return a SortingCollection
    */
  def sortingCollection(order: SortOrder,
                        header: SAMFileHeader,
                        maxInMemory: Int = 1e6.toInt,
                        tmpDir: Option[DirPath] = None): SortingCollection[SAMRecord] = {
    SortingCollection.newInstance(
      classOf[SAMRecord],
      new BAMRecordCodec(header),
      order.getComparatorInstance,
      maxInMemory,
      tmpDir.map(_.toFile).orNull
    )
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
  def queryGroupedIterator(in: SamReader,
                           maxInMemory: Int = 1e6.toInt,
                           tmpDir: Option[DirPath] = None): BetterBufferedIterator[SAMRecord] = {
    if (in.getFileHeader.getSortOrder == SortOrder.queryname || in.getFileHeader.getGroupOrder == GroupOrder.query) {
      in.iterator().bufferBetter
    }
    else {
      logger.info("Sorting into queryname order.")
      val progress = new ProgressLogger(this.logger, "Queryname sorted")
      val sorter = sortingCollection(SortOrder.queryname, in.getFileHeader, maxInMemory, tmpDir)
      in.foreach { rec =>
        sorter.add(rec)
        progress.record(rec)
      }
      sorter.iterator().bufferBetter
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
  def templateIterator(in: SamReader,
                       maxInMemory: Int = 1e6.toInt,
                       tmpDir: Option[DirPath] = None): Iterator[Template] = {
    val queryIterator: BetterBufferedIterator[SAMRecord] = queryGroupedIterator(in, maxInMemory, tmpDir)

    new Iterator[Template] {
      override def hasNext: Boolean = queryIterator.hasNext
      override def next(): Template = {
        require(hasNext, "next() called on empty iterator")
        Template(queryIterator)
      }
    }
  }

  /**
    * Ensures that any NM/UQ/MD tags on the read are accurate.  If the read is unmapped, any existing
    * values are removed.  If the read is mapped all three tags will have values regenerated.
    *
    * @param rec the SAMRecord to update
    * @param ref a reference sequence file walker to pull the reference information from
    */
  def regenerateNmUqMdTags(rec: SAMRecord, ref: ReferenceSequenceFileWalker): Unit = {
    import SAMTag._
    if (rec.getReadUnmappedFlag) {
      rec.setAttribute(NM.name(), null)
      rec.setAttribute(UQ.name(), null)
      rec.setAttribute(MD.name(), null)
    }
    else {
      val refBases = ref.get(rec.getReferenceIndex).getBases
      SequenceUtil.calculateMdAndNmTags(rec, refBases, true, true)
      if (rec.getBaseQualities != null && rec.getBaseQualities.length != 0) {
        rec.setAttribute(SAMTag.UQ.name, SequenceUtil.sumQualitiesOfMismatches(rec, refBases, 0))
      }
    }
  }
}
