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
  require((leftStart == -1) == (leftEnd == -1))
  require((rightStart == -1) == (rightEnd == -1))
  require(leftStart != -1 || rightStart != -1)
  require(leftStart == -1 || leftStart <= leftEnd, f"leftStart is > leftEnd: $this")
  require(rightStart == -1 || rightStart <= rightEnd, f"rightStart is > rightEnd: $this")
  require(leftStart == -1 || rightStart == -1 || leftStart <= rightStart, f"leftStart is > rightStart: $this")

  @inline def leftStart: Int  = left_start
  @inline def leftEnd: Int    = left_end
  @inline def rightStart: Int = right_start
  @inline def rightEnd: Int   = right_end
  @inline def contig: String  = chrom
  @inline def start: Int      = if (leftStart == -1) rightStart else leftStart
  @inline def end: Int        = if (rightEnd == -1) leftEnd else rightEnd

  def leftPrimerLength: Int       = CoordMath.getLength(leftStart, leftEnd)
  def rightPrimerLength: Int      = CoordMath.getLength(rightStart, rightEnd)
  def longestPrimerLength: Int    = Math.max(leftPrimerLength, rightPrimerLength)
  def leftPrimerLocation: String  = f"$chrom:$leftStart-$leftEnd"
  def rightPrimerLocation: String = f"$chrom:$rightStart-$rightEnd"
  def ampliconLocation: String    = f"$chrom:$leftStart-$rightEnd"
  def identifier: String          = this.id.getOrElse(ampliconLocation)
}
