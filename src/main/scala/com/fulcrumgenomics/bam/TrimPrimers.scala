/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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
import com.fulcrumgenomics.bam.SamRecordClipper.ClippingMode
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{DelimitedDataParser, Io}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.util._
import Math.abs

import scala.collection.JavaConversions.iterableAsScalaIterable

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Trims primers from reads post-alignment.  Takes in a BAM file of aligned reads
    |and a tab-delimited file with five columns (chrom, left_start, left_stop, right_start,
    |and right_stop) which provide the 1-based inclusive start and end positions of the
    |primers for each amplicon.
    |
    |Paired end reads that map to a given amplicon position are trimmed so that the
    |alignment no-longer includes the primer sequences.
    |
    |All other aligned reads have the maximum primer length trimmed!
  """)
class TrimPrimers
( @arg(flag="i", doc="Input BAM file.")  val input: PathToBam,
  @arg(flag="o", doc="Output BAM file.") val output: PathToBam,
  @arg(flag="p", doc="File with primer locations.") val primers: FilePath,
  @arg(flag="H", doc="If true, hard clip reads, else soft clip.") val hardClip: Boolean = false,
  @arg(flag="S", doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
  @arg(flag="s", doc="Sort order of output BAM file (defaults to input sort order).") val sortOrder: Option[SortOrder] = None
)extends FgBioTool with LazyLogging {
  private val mode = if (hardClip) ClippingMode.Hard else ClippingMode.Soft

  Io.assertReadable(input)
  Io.assertReadable(primers)
  Io.assertCanWriteFile(output)

  /** A Locatable Amplicon class. */
  private case class Amplicon(chrom: String, leftStart: Int, leftEnd: Int, rightStart: Int, rightEnd: Int) extends Locatable {
    def leftPrimerLength: Int    = CoordMath.getLength(leftStart, leftEnd)
    def rightPrimerLength: Int   = CoordMath.getLength(rightStart, rightEnd)
    def longestPrimerLength: Int = Math.max(leftPrimerLength, rightPrimerLength)

    override def getContig: String = chrom
    override def getStart: Int = leftStart
    override def getEnd: Int = rightEnd
  }

  override def execute(): Unit = {
    val in = SamReaderFactory.make().open(input)
    val outSortOrder = sortOrder.getOrElse(in.getFileHeader.getSortOrder)
    val outHeader = in.getFileHeader.clone()
    outHeader.setSortOrder(outSortOrder)
    val out = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(outHeader, outSortOrder == SortOrder.queryname, output.toFile, null)

    val detector = loadPrimerFile(primers)
    val maxPrimerLength = detector.getAll.map(_.longestPrimerLength).max

    val iterator = queryNameOrderIterator(in)
    while (iterator.hasNext) {
      val reads = nextTemplate(iterator)
      val rec1 = reads.find(r => r.getReadPairedFlag && r.getFirstOfPairFlag  && !r.isSecondaryOrSupplementary)
      val rec2 = reads.find(r => r.getReadPairedFlag && r.getSecondOfPairFlag && !r.isSecondaryOrSupplementary)

      (rec1, rec2) match {
        case (Some(r1), Some(r2)) =>
          // FR mapped pairs get the full treatment
          if (!r1.getReadUnmappedFlag && !r2.getReadUnmappedFlag &&
            r1.getReferenceName == r2.getReferenceName &&
            SamPairUtil.getPairOrientation(r1) == PairOrientation.FR) {

            val (left, right) = if (r1.getReadNegativeStrandFlag) (r2, r1) else (r1, r2)
            val (start, end ) = (left.getUnclippedStart, right.getUnclippedEnd)
            val insert = new Interval(left.getContig, start, end)
            val amp = detector.getOverlaps(insert).find(amp => abs(amp.leftStart - start) <= slop && abs(amp.rightEnd - end) <= slop)
            amp match {
              case Some(amplicon) =>
                val leftClip  = amplicon.leftPrimerLength
                val rightClip = amplicon.rightPrimerLength
                reads.foreach { rec =>
                  val toClip = if (rec.getFirstOfPairFlag == left.getFirstOfPairFlag) leftClip else rightClip
                  SamRecordClipper.clip5PrimeEndOfRead(rec, toClip, mode)
                }
              case None =>
                reads.foreach(r => SamRecordClipper.clip5PrimeEndOfRead(r, maxPrimerLength, mode))
            }

            clipFullyOverlappedFrReads(r1, r2)
          }
          // Pairs without both reads mapped in FR orientation are just maximally clipped
          else {
            reads.foreach(r => SamRecordClipper.clip5PrimeEndOfRead(r, maxPrimerLength, mode))
          }

          SamPairUtil.setMateInfo(r1, r2, true)
          reads.filter(_.getSupplementaryAlignmentFlag).foreach { rec =>
            val mate = if (rec.getFirstOfPairFlag) r2 else r1
            SamPairUtil.setMateInformationOnSupplementalAlignment(rec, mate, true);
          }
        case _ =>
          // Just trim each read independently
          reads.foreach(r => SamRecordClipper.clip5PrimeEndOfRead(r, maxPrimerLength, mode))
      }

      reads.foreach(out.addAlignment)
    }

    out.close()
  }

  /** Gets an iterator in query name order over the records. */
  private def queryNameOrderIterator(in: SamReader): BufferedIterator[SAMRecord] = {
    if (in.getFileHeader.getSortOrder == SortOrder.queryname) {
      in.iterator().bufferBetter
    }
    else {
      logger.info("Sorting into queryname order.")
      val sorter = SortingCollection.newInstance[SAMRecord](classOf[SAMRecord], new BAMRecordCodec(in.getFileHeader), new SAMRecordQueryNameComparator(), 1e6.toInt);
      in.foreach(sorter.add)
      sorter.iterator().bufferBetter
    }
  }

  /** Fetches the next group of records that all share the same readname/template from the iterator. */
  private def nextTemplate(iterator: BufferedIterator[SAMRecord]): Seq[SAMRecord] = {
    val first    = iterator.next()
    val template = first.getReadName
    first :: iterator.takeWhile(_.getReadName == template).toList
  }

  /** Creates an overlap detector for all the amplicons from the input file. */
  private def loadPrimerFile(path: FilePath): OverlapDetector[Amplicon] = {
    val parser = DelimitedDataParser(path, '\t')
    val detector = new OverlapDetector[Amplicon](0,0)
    parser.foreach { row =>
      val amp = Amplicon(
        chrom = row[String]("chrom"),
        leftStart  = row[Int]("left_start"),
        leftEnd    = row[Int]("left_end"),
        rightStart = row[Int]("right_start"),
        rightEnd   = row[Int]("right_end")
      )

      detector.addLhs(amp, amp)
    }

    detector
  }

  /**
    * Adapted from Picard's AbstractAlignmentMerger - takes in mapped FR pairs only
    * and clips from the 3' ends of the reads if the reads are fully overlapped
    * and extend past each other's starts.
    */
  private def clipFullyOverlappedFrReads(r1: SAMRecord, r2: SAMRecord): Unit = {
    val (plus, minus) = if (r1.getReadNegativeStrandFlag) (r2,r1) else (r1, r2)

    if (plus.getAlignmentStart < minus.getAlignmentEnd) {
      val plusTrim  = plus.getAlignmentEnd   - minus.getAlignmentEnd
      val minusTrim = plus.getAlignmentStart - minus.getAlignmentStart

      if (plusTrim  > 0) SamRecordClipper.clip3PrimeEndOfAlignment(plus, plusTrim, mode)
      if (minusTrim > 0) SamRecordClipper.clip3PrimeEndOfAlignment(minus, minusTrim, mode)
    }
  }
}
