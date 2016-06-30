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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.util
import htsjdk.samtools.SAMUtils

import scala.math._

object PhredScore {
  /** The maximum Phred Value. */
  val MaxValue: Int = SAMUtils.MAX_PHRED_SCORE

  /** The ASCII representation of a zero phred score. */
  val PhredZeroChar: Char = phredScoreToChar(0)

  /** The precision to use when converting to an Integer. */
  val Precision = 0.001

  val ZeroProbability: Double = Double.PositiveInfinity

  val OneProbability: Double = 0.0

  /** The value of log(10.0). */
  private val LogOf10: Double= log(10.0)

  /** Converts a probability (log-space) to a phred score. */
  def toPhredScore(probability: LogDouble): Double = -10.0 * (probability.logValue / LogOf10)

  /** Converts a probability to a phred-score as an Int. If the score is at most `Precision` away
    * from the next higher score, it will return the higher score*/
  def toInt(probability: LogDouble): Int = {
    val upper = (toPhredScore(probability=probability)+Precision).toInt
    val score = toPhredScore(probability=probability).toInt
    if (upper > score) upper
    else score
  }

  /** Converts to the ASCII representation of this phred score as a Byte */
  def toByte(probability: LogDouble): Byte = SAMUtils.phredToFastq(toInt(probability = probability)).toByte

  /** Converts to the ASCII representation of this phred score as a Char */
  def toChar(probability: LogDouble): Char = SAMUtils.phredToFastq(Math.min(MaxValue, toInt(probability = probability)))

  /** Converts the phred score to the ASCII representation. */
  def phredScoreToChar(phredScore: Int): Char = SAMUtils.phredToFastq(phredScore)

  /** Converts a probability to a phred-score as an integer. */
  def toString(probability: LogDouble): String = toInt(probability=probability).toString

  /** Converts a phred score to a probability (log-space). */
  def fromPhredScore(phredScore: Double): LogDouble = new LogDouble(logValue = LogOf10 * phredScore / -10)

  implicit class DoubleWrapper(phredScore: Double) {
    def fromPhredScore: LogDouble = PhredScore.fromPhredScore(phredScore=phredScore)
  }
  implicit class IntWrapper(phredScore: Int) {
    def fromPhredScore: LogDouble = PhredScore.fromPhredScore(phredScore = phredScore)
  }
  implicit class LogDoubleWrapper(probability: LogDouble) {
    def toPhredScore: Double       = PhredScore.toPhredScore(probability=probability)
    def toPhredScoreInt: Int       = PhredScore.toInt(probability=probability)
    def toPhredScoreByte: Byte     = PhredScore.toByte(probability=probability)
    def toPhredScoreChar: Char     = PhredScore.toChar(probability=probability)
    def toPhredScoreString: String = PhredScore.toString(probability=probability)
  }
}

object LogDouble {

  /** log(X) == 0. */
  private[util] val Zero = new LogDouble(logValue=Double.NegativeInfinity)

  /** log(X) == 1. */
  private[util] val One = new LogDouble(logValue=0.0)

  /** Converts the double to log-space. */
  def toLogDouble(value: Double): LogDouble = {
    if (value < 0) throw new IllegalArgumentException("Cannot use LogDouble to store values less than zero.")
    new LogDouble(logValue=log(value))
  }

  /** Gets the mean of the given phred values. */
  def mean(logDoubles: LogDouble*): LogDouble = logDoubles.foldLeft(Zero)((sum, qual) => sum + qual) / toLogDouble(logDoubles.length)

  /** Converts a Double to log-space. */
  implicit class DoubleValue(v: Double) {
    def toLogDouble: LogDouble = LogDouble.toLogDouble(value=v)
  }

  /** Converts a string to log-space. */
  implicit class StringValue(doubleValue: String) {
    def toLogDouble: LogDouble = LogDouble(doubleValue=doubleValue)
  }

  /** String constructor when a double (not in log-space) is given. */
  def apply(doubleValue: String): LogDouble = {
    new util.LogDouble(logValue=log(if (doubleValue.contains(".")) doubleValue.toDouble else doubleValue.toInt))
  }

  /** Convenience constructor. */
  def apply(logValue: Double): LogDouble = new LogDouble(logValue=logValue)

  /** Implementation of addition. */
  def add(a: Double, b: Double): Double = {
    if (a.isNegInfinity) b
    else if (b.isNegInfinity) a
    else if (b < a) add(b, a)
    else a + log1p(exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 + exp(b - a))`
  }

  /** Implementation of subtraction. */
  def sub(a: Double, b: Double): Double = {
    if (b.isNegInfinity) a
    else if (a == b) Double.NegativeInfinity
    else if (a < b) throw new IllegalArgumentException("Subtraction will be less than zero.")
    else a + log1p(0.0 - exp(b - a)) // using log1p, otherwise it would be `a + log(1.0 - exp(b - a))`
  }

  def addAll(values: Double*): Double = {
    val minValue = values.min
    if (minValue.isNegInfinity) minValue
    else {
      val minValueIndex = values.indexOf(minValue)
      var sum: Double = 0.0
      var i = 0
      while (i < values.length) {
        if (minValueIndex != i) sum += exp(values(i) - minValue)
        i += 1
      }
      minValue + log1p(sum)
    }
  }
}

/**
  * Represents a Double in log-space, to prevent underflow and overflow.
  * @param logValue the Double in log-space
  */
class LogDouble(val logValue: Double) extends AnyVal with Ordered[LogDouble] {


  /** Returns the value as a Double. */
  def linearValue: Double = exp(logValue)

  /** Division. */
  def /(that: LogDouble): LogDouble = {
    if (that.logValue == LogDouble.Zero.logValue) throw new IllegalStateException("Trying to divide by 0")
    new LogDouble(logValue=this.logValue - that.logValue)
  }

  /** Multiplication. */
  def * (that: LogDouble): LogDouble = new LogDouble(logValue=this.logValue + that.logValue)

  /** Addition. */
  def + (that: LogDouble): LogDouble = new LogDouble(logValue=LogDouble.add(this.logValue, that.logValue))

  /** Subtraction. */
  def - (that: LogDouble): LogDouble = new LogDouble(logValue=LogDouble.sub(this.logValue, that.logValue))

  /** 1.0 - this; useful for probabilities. */
  def oneMinus(): LogDouble           = {
    if (LogDouble.One < this) LogDouble.Zero
    else if (this < LogDouble.Zero) LogDouble.One
    else LogDouble.One - this
  }

  /** Compares two doubles in log space. */
  def compare(that: LogDouble): Int = this.logValue.compare(that.logValue)

  override def toString: String = s"$linearValue:$logValue"
}