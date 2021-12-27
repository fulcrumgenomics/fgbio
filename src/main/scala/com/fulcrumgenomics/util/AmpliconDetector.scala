/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.util

import java.lang.Math.abs

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamRecord
import htsjdk.samtools.util.{Interval, OverlapDetector}
import scala.collection.compat._


object AmpliconDetector {
  def apply(amplicons: IterableOnce[Amplicon], slop: Int, unclippedCoordinates: Boolean): AmpliconDetector = {
    new AmpliconDetector(
      detector             = Amplicon.overlapDetector(amplicons=amplicons.iterator),
      slop                 = slop,
      unclippedCoordinates = unclippedCoordinates
    )
  }

  def apply(path: FilePath, slop: Int, unclippedCoordinates: Boolean): AmpliconDetector = {
    this.apply(amplicons=Metric.iterator[Amplicon](path), slop=slop, unclippedCoordinates=unclippedCoordinates)
  }
}

/** Detector to find the amplicon sequenced by a given read for a multiplex PCR (or similar) experiment.  The primer for
  * a read is assumed to be at the 5' end (sequencing order) of the read.  The various methods use the genomic
  * coordinates of the amplicons compared to the 5' end (sequencing order) of the read.  The genomic coordinates are
  * adjusted based on clipped (hard or soft) bases if `unclippedCoordinates` is true.
  *
  * @param detector the overlap detector across all amplicons
  * @param slop match to primer locations +/- this many bases
  * @param unclippedCoordinates true to adjust the genomic coordinates of the 5' end based on clipping (hard or soft)
  */
class AmpliconDetector(val detector: OverlapDetector[Amplicon],
                       val slop: Int,
                       unclippedCoordinates: Boolean) {

  /** The maximum primer length across all amplicons. */
  val maxPrimerLength: Int = detector.getAll.map(_.longestPrimerLength).max

  /** The class used to calculate coordinates, with or with clipping adjustment. */
  private val coordinates: SamRecordCoordinates = {
    if (unclippedCoordinates) new UnclippedCoordinates else new ClippedCoordinates
  }

  /** Finds the amplicon for a single read (not read pair).
    *
    * The primer is assumed to be located at the start of the read in 5' sequencing order.  Therefore, a positive strand
    * read will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
    * strand read will use its aligned end position to match against the amplicon's right-most coordinate.  The mate of
    * the read will not be used.
    *
    * @param refName the reference name
    * @param start the start position
    * @param end the end position
    * @param positiveStrand true for the positive strand, false for the negative strand
    * @return the amplicon if found
    */
  private[util] def findPrimer(refName: String, start: Int, end: Int, positiveStrand: Boolean): Option[Amplicon] = {
    val interval = new Interval(refName, start, end)
    val hits: Iterator[(Amplicon, Int)] = {
      if (positiveStrand) detector.getOverlaps(interval).map(amp => (amp, abs(amp.leftStart - start)))
      else detector.getOverlaps(interval).map(amp => (amp, abs(amp.rightEnd - end)))
    }
    hits.minByOption(_._2).find(_._2 <= slop).map(_._1)
  }

  /** Finds the amplicon that contains the primer for the given read.
    *
    * The primer is assumed to be located at the start of the read in 5' sequencing order.  Therefore, a positive strand
    * read will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
    * strand read will use its aligned end position to match against the amplicon's right-most coordinate.  The mate of
    * the read will not be used.
    *
    * If multiple amplicons are found, the amplicon with the smallest distance will be returned.
    *
    * @param rec the record to match.
    *
    * @return Some amplicon if a match is found, None otherwise
    */
  def findPrimer(rec: SamRecord): Option[Amplicon] = if (rec.unmapped) None else {
    val (start, end) = (coordinates.start(rec), coordinates.end(rec))
    findPrimer(refName=rec.refName, start=start, end=end, positiveStrand=rec.positiveStrand)
  }

  /** Finds the amplicon that contains the primer for the given read's mate.
    *
    * The primer is assumed to be located at the start of the mate in 5' sequencing order.  Therefore, a positive strand
    * mate will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
    * strand mate will use its aligned end position to match against the amplicon's right-most coordinate.  The current
    * read will not be used.
    *
    * If multiple amplicons are found, the amplicon with the smallest distance will be returned.
    *
    * @param rec the record for which to match the mate
    * @return Some amplicon if a match is found for the mate, None otherwise
    */
  def findMatePrimer(rec: SamRecord): Option[Amplicon] = if (rec.mateUnmapped) None else {
    val (start, end) = (coordinates.mateStart(rec), coordinates.mateEnd(rec))
    findPrimer(refName=rec.mateRefName, start=start, end=end, positiveStrand=rec.matePositiveStrand)
  }

  /** Finds the amplicon for a read pair.
    *
    * Both reads must be mapped to the same contig and in FR orientation.
    *
    * The primer is assumed to be located at the start of the read in 5' sequencing order.  Therefore, a positive strand
    * read will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
    * strand read will use its aligned end position to match against the amplicon's right-most coordinate.  Similarly
    * for its mate.  Both the read and its mate must match the same amplicon.
    *
    * If multiple amplicons are found, the amplicon with the smallest combined distance between coordinates will be
    * returned.
    *
    * @param r1 the first end of a pair
    * @param r2 the second end of a pair
    * @return Some amplicon if a match is found for the read pair, None otherwise
    */
  def find(r1: SamRecord, r2: SamRecord): Option[Amplicon] = if (!r1.isFrPair) None else {
    val (left, right) = if (r1.positiveStrand) (r1, r2) else (r2, r1)
    val (start, end)  = (coordinates.start(left), coordinates.end(right))
    val insert        = new Interval(left.refName, start, end)
    detector.getOverlaps(insert)
      .map(amp => (amp, abs(amp.leftStart - start), abs(amp.rightEnd - end)))
      .filter(hit => hit._2 <= slop && hit._3 <= slop)
      .minByOption(hit => (hit._2 + hit._3))
      .map(_._1)
  }
}

/** Base trait for classes that compute the genomic coordinates for a record. */
private sealed trait SamRecordCoordinates {
  def start(rec: SamRecord): Int
  def end(rec: SamRecord): Int
  def mateStart(rec: SamRecord): Int
  @inline final def mateEnd(rec: SamRecord): Int = this.mateEndOption(rec=rec).getOrElse(invalidNoMateCigar(rec))
  protected def mateEndOption(rec: SamRecord): Option[Int]
  @inline final protected def invalidNoMateCigar(rec: SamRecord): Nothing = {
    throw new IllegalArgumentException(s"Record missing mate cigar (MC): ${rec.name}")
  }
}

/** Computes the genomic coordinates based on the aligned bases. */
private class ClippedCoordinates extends SamRecordCoordinates {
  @inline def start(rec: SamRecord): Int     = rec.start
  @inline def end(rec: SamRecord): Int       = rec.end
  @inline def mateStart(rec: SamRecord): Int = rec.mateStart
  @inline protected def mateEndOption(rec: SamRecord): Option[Int] = rec.mateEnd
}

/** Computes the genomic coordinates adjusting for clipped bases */
private class UnclippedCoordinates extends SamRecordCoordinates {
  @inline def start(rec: SamRecord): Int     = rec.unclippedStart
  @inline def end(rec: SamRecord): Int       = rec.unclippedEnd
  @inline def mateStart(rec: SamRecord): Int = rec.mateUnclippedStart.getOrElse(invalidNoMateCigar(rec))
  @inline protected def mateEndOption(rec: SamRecord): Option[Int] = rec.mateUnclippedEnd
}
