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

import com.fulcrumgenomics.FgBioDef.FilePath
import htsjdk.samtools.util.{CoordMath, OverlapDetector}

/** Convenience methods for [[Amplicon]] */
object Amplicon {
  /** Builds a [[OverlapDetector]] for the given amplicons. */
  def overlapDetector(amplicons: Iterator[Amplicon]): OverlapDetector[Amplicon] = {
    val detector = new OverlapDetector[Amplicon](0,0)
    amplicons.foreach(amp => detector.addLhs(amp, amp))
    detector
  }

  /** Builds a [[OverlapDetector]] for the given file of amplicons. */
  def overlapDetector(path: FilePath): OverlapDetector[Amplicon] = {
    Amplicon.overlapDetector(amplicons=Metric.iterator[Amplicon](path))
  }
}


// Developer note: the left/right start/end members are snakecase, which is the convention for  "Metric" files
// (originally required by [[TrimPrimers]]).  The are made private so that code that uses [[Amplicon]] uses the
// camel case instead (leftStart/leftEnd/rightStart/rightEnd).
/** A Locatable Amplicon class.
  * @param chrom the chromosome for the amplicon
  * @param left_start the 1-based start position of the left-most primer
  * @param left_end the 1-based end position inclusive of the left-most primer
  * @param right_start the 1-based start position of the right-most primer
  * @param right_end the 1-based end position inclusive of the right-most primer
  */
case class Amplicon
( chrom: String,
  private val left_start: Int,
  private val left_end: Int,
  private val right_start: Int,
  private val right_end: Int,
  private val id: Option[String] = None
) extends GenomicSpan with Metric {
  @inline def contig: String  = chrom
  @inline def leftStart: Option[Int]  = if (left_start > -1) Some(left_start) else None 
  @inline def leftEnd: Option[Int]    = if (left_end > -1) Some(left_end)  else None
  @inline def rightStart: Option[Int] = if (right_start > -1) Some(right_start) else None
  @inline def rightEnd: Option[Int]   = if (right_end > -1) Some(right_end) else None 
  val (s, e, longest_primer_length, left_primer_length, right_primer_length, left_primer_location, right_primer_location) = (leftStart, leftEnd, rightStart, rightEnd) match {
    case (Some(left_start), Some(left_end), Some(right_start), Some(right_end)) =>
      require(left_start <= left_end, f"leftStart is > leftEnd: $this")
      require(right_start <= right_end, f"rightStart is > rightEnd: $this")
      require(left_start <= right_start, f"leftStart is > rightStart: $this")

      def leftPrimerLength: Int       = CoordMath.getLength(left_start, left_end)
      def rightPrimerLength: Int      = CoordMath.getLength(right_start, right_end)
      def longestPrimerLength: Int    = Math.max(leftPrimerLength, rightPrimerLength)
      def leftPrimerLocation: Option[String]  = Some(f"$chrom:$left_start-$left_end")
      def rightPrimerLocation: Option[String] = Some(f"$chrom:$right_start-$right_end")
      def ampliconLocation: String    = f"$chrom:$left_start-$right_end"
      def identifier: String          = this.id.getOrElse(ampliconLocation)
      (left_start, right_end, Math.max(leftPrimerLength, rightPrimerLength), leftPrimerLength, rightPrimerLength, leftPrimerLocation, rightPrimerLocation)
    case (Some(left_start), Some(left_end), None, None) =>
      require(left_start <= left_end, f"leftStart is > leftEnd: $this")
      def leftPrimerLength: Int       = CoordMath.getLength(left_start, left_end)
      def leftPrimerLocation: Option[String]  = Some(f"$chrom:$left_start-$left_end")
      (left_start, left_end, leftPrimerLength, leftPrimerLength, 0, leftPrimerLocation, None)
    case (None, None, Some(right_start), Some(right_end)) => 
      require(right_start <= right_end, f"rightStart is > rightEnd: $this")

      def rightPrimerLength: Int      = CoordMath.getLength(right_start, right_end)
      def rightPrimerLocation: Option[String] = Some(f"$chrom:$right_start-$right_end")
      (right_start, right_end, rightPrimerLength, 0, rightPrimerLength, None, rightPrimerLocation)
    case _ =>
      throw new Exception(s"At least (left_start and left_end) or (right_start and right_end) need to be set in every row of the primer file.")
  }
  @inline def start: Int = s
  @inline def end: Int = e
  def leftPrimerLength: Int       = left_primer_length
  def rightPrimerLength: Int      = right_primer_length
  def longestPrimerLength: Int    = longest_primer_length
  def leftPrimerLocation: Option[String]  = left_primer_location
  def rightPrimerLocation: Option[String] = right_primer_location
  def ampliconLocation: String    = f"$chrom:$start-$end"
  def identifier: String          = this.id.getOrElse(ampliconLocation)
}
