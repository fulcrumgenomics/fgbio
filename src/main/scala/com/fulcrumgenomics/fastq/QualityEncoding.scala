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

package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.{NumericCounter, SimpleCounter}

/** Trait for describing quality score encoding systems. */
sealed trait QualityEncoding {
  /** The valid numeric range of the quality scores, prior to any ascii encoding. */
  def numericRange: Range

  /** The shift to go from a quality score to the ascii encoding of the quality score. */
  def asciiOffset: Int

  /** The valid range of ascii characters (as ints). */
  lazy val asciiRange: Range = Range.inclusive(numericRange.start + asciiOffset, numericRange.end + asciiOffset)

  /** Must be implemented to convert from a numeric quality score in the given format to a numeric quality score
    * in standard format (i.e. -10*log10(pError)).
    */
  def toStandardNumeric(q: Byte): Byte

  /** Translates a character/ascii-encoded quality to a numeric quality in standard phred scaling. */
  @inline final def toStandardNumeric(ch: Char): Byte = toStandardNumeric((ch.toByte - asciiOffset).toByte)

  /** Converts a String of quality scores in the given format into standard numeric quality scores. */
  def toStandardNumeric(qs: String): Array[Byte] = {
    val bs = new Array[Byte](qs.length)
    forloop (from=0, until=bs.length) { i => bs(i) = toStandardNumeric(qs(i)) }
    bs
  }

  /** Converts a character from one ascii encoding another */
  def toStandardAscii(ch: Char): Char = (toStandardNumeric(ch) + QualityEncoding.Standard.asciiOffset).toChar

  /** Converts a String in one ascii encoding to another. */
  def toStandardAscii(qs: String): String = qs.map(toStandardAscii)
}

object QualityEncoding {
  /** Quality encoding for the legacy Solexa scores used pre version 1.3.
    * Q[Solexa] = −10 * log10[(P/1-P)]
    */
  case object Solexa extends QualityEncoding {
    override val numericRange: Range = Range.inclusive(-5, 62)
    override val asciiOffset : Int   = 64
    private  val mappingToStandard   = numericRange.toArray.map(q => (10 * Math.log10(1 + Math.pow(10, q / 10.0))).round.toByte)
    override def toStandardNumeric(q: Byte): Byte = mappingToStandard(q+5)
  }

  /** Quality encoding for Illumina, 1.3 <= version < 1.8; Phred scaling + 64. */
  case object Illumina extends QualityEncoding {
    override val numericRange: Range = Range.inclusive(0, 62)
    override val asciiOffset : Int   = 64
    override def toStandardNumeric(q: Byte): Byte = q
  }

  /** Quality encoding for Sanger & Illumina 1.8+; Phred scaling + 33. */
  case object Standard extends QualityEncoding {
    override val numericRange: Range = Range.inclusive(0, 93)
    override val asciiOffset : Int   = 33
    override def toStandardNumeric(q: Byte):Byte = q
    override def toStandardAscii(ch: Char): Char = ch
    override def toStandardAscii(qs: String): String = qs
  }

  /** All possible values of quality encoding. */
  val all: List[QualityEncoding] = Standard :: Illumina :: Solexa :: Nil
}

/**
  * Class to detect the quality encoding in use in FASTQ files.
  */
class QualityEncodingDetector {
  private val counter = new SimpleCounter[Char]

  /** Adds an ascii quality character to the set of evidence. */
  def add(ch: Char): Unit = counter.count(ch)

  /** Adds a String of ascii quality characters to the evidence. */
  def add(chs: String): Unit = forloop(from=0, until=chs.length) { i => add(chs.charAt(i)) }

  /**
    * Samples strings from the iterator until at least `count` qualities have been sampled.
    * @param ss an iterator of Strings where each string is an ascii string of quality scores
    * @param count the number of quality scores (not strings!) to attempt to sample
    */
  def sample(ss: Iterator[String], count: Int = 100000): Unit = {
    val goal = counter.total + count
    while (counter.total < goal && ss.hasNext) add(ss.next())
  }

  /** Returns the set of quality formats that are compatible with the qualities seen thus far. */
  def compatibleEncodings: List[QualityEncoding] = {
    QualityEncoding.all.filter(enc => counter.forall { case (ch, count) => enc.asciiRange.contains(ch.toInt)} )
  }

  /** Returns true if the observed qualities are compatible with the encoding, false otherwise. */
  def isCompatible(encoding: QualityEncoding): Boolean = compatibleEncodings.contains(encoding)

  /** Returns true if more than one quality encoding is compatible given the evidence so far. */
  def isAmbiguous: Boolean = compatibleEncodings.size > 1

  /** Returns the set of compatible encodings ranked by how close the mean value is to Q under each encoding. */
  def rankedCompatibleEncodings(q: Int): List[QualityEncoding] = {
    val compatible = compatibleEncodings
    if (compatible.size < 2 || this.counter.isEmpty) {
      compatible
    }
    else {
      compatible.map { enc =>
        val x = new NumericCounter[Int]()
        this.counter.foreach { case (ch, count) => x.count(enc.toStandardNumeric(ch), count) }
        (enc, math.abs(x.mean() - q))
      }.sortBy(_._2).map(_._1)
    }
  }
}
