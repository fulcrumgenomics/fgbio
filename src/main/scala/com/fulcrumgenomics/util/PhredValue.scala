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

import dagr.commons.CommonsDef._

/** Methods and values to help convert to and from Phred-scaled probabilities.
  *
  * A number of implicit conversions are defined in this object.  Integers are implicitly are interpreted as
  * literal phred values, while Floats and Doubles converted to Phred values using `toPhred`.
  */
object PhredValue {
  /** Phred-scaled zero probability. */
  val ZeroProbability = PhredValue(Double.PositiveInfinity)
  /** Phred-scaled probability of one. */
  val OneProbability  = PhredValue(0.0)
  /** Phred value of 0. */
  val ZeroPhred       = PhredValue(0.0)
  /** The precision to use when converting to an Integer. */
  val Precision       = 0.0000000001 // Q100
  /** converts a probability to a phred score, with no validation. */
  def toPhred(pr: Double) = -10.0*Math.log10(pr)
  /** converts to a probability from a phred score, with no validation. */
  def toPr(phred: Double) = Math.pow(10.0, phred / -10.0)
  /** Gets the mean of the given phred values. */
  def mean(phred: PhredValue*): PhredValue = phred.foldLeft(ZeroProbability)((sum, qual) => sum + qual) / phred.length.toDouble
  /*
  implicit def phredValueWrapper[T](v: T)(implicit numeric: Numeric[T]): PhredValue = new PhredValue(value=numeric.toDouble(v))
  implicit class NumericPhredValue[T](v: T)(implicit numeric: Numeric[T]) {
    def toPhredValue: PhredValue = new PhredValue(value=numeric.toDouble(v))
  }
  */
  implicit def intToPhredWrapper(v: Int): PhredValue = new PhredValue(value=v)
  implicit def floatToPhredWrapper(v: Float): PhredValue = {
    //if (v < 0 || v > 1) throw new IllegalArgumentException(s"Float value is out of range to convert to phred: $v")
    PhredValue(toPhred(v))
  }
  implicit def doubleToPhredWrapper(v: Double): PhredValue = {
    //if (v < 0 || v > 1) throw new IllegalArgumentException(s"Double value is out of range to convert to phred: $v")
    PhredValue(toPhred(v))
  }
  implicit class IntPhredValue(v: Int) {
    def toPhredValue: PhredValue = new PhredValue(value=v.toDouble)
  }
  implicit class FloatPhredValue(v: Float) {
    def toPhredValue: PhredValue = {
      //if (v < 0 || v > 1) throw new IllegalArgumentException(s"Float value is out of range to convert to phred: $v")
      new PhredValue(value=toPhred(v.toDouble))
    }
  }
  implicit class DoublePhredValue(v: Double) {
    def toPhredValue: PhredValue = {
      //if (v < 0 || v > 1) throw new IllegalArgumentException(s"Double value is out of range to convert to phred: $v")
      new PhredValue(value=toPhred(v.toDouble))
    }
  }
}

/** Represents Phred-scaled probabilities, which are `-10 * log10(Pr)`. */
case class PhredValue(var value: Double) extends Ordered[PhredValue] {
  import PhredValue._
  def this(v: String) = this(if (v.contains(".")) v.toDouble else v.toInt)
  //if (value < 0) throw new IllegalArgumentException("PhredValue was less than zero")
  def /(that: PhredValue): PhredValue = {
    if (that.value == Double.PositiveInfinity) throw new IllegalStateException("Trying to divide by 0")
    PhredValue(this.value - that.value)
  }
  def /=(that: PhredValue): PhredValue = {
    if (that.value == Double.PositiveInfinity) throw new IllegalStateException("Trying to divide by 0")
    yieldAndThen(this)(this.value -= that.value)
  }
  def * (that: PhredValue): PhredValue = PhredValue(this.value + that.value)
  def *=(that: PhredValue): PhredValue = yieldAndThen(this)(this.value += that.value)
  def + (that: PhredValue): PhredValue = PhredValue(add(this.value, that.value))
  def +=(that: PhredValue): PhredValue = yieldAndThen(this)(add(this.value, that.value))
  def - (that: PhredValue): PhredValue = PhredValue(sub(this.value, that.value))
  def -=(that: PhredValue): PhredValue = yieldAndThen(this)(sub(this.value, that.value))
  def inv(): PhredValue                = PhredValue(toPhred(1.0 - toPr(this.value)))
  def invEq: PhredValue                = yieldAndThen(this)(this.value = toPhred(1.0 - toPr(this.value)))
  def compare(that: PhredValue): Int   = this.value.compare(that.value)
  private def add(a: Double, b: Double): Double = {
    // from phred to log10(pr)
    (a / -10.0, b / -10.0) match {
      case (aLog10Pr, bLog10Pr) if aLog10Pr == Double.NegativeInfinity => b
      case (aLog10Pr, bLog10Pr) if bLog10Pr == Double.NegativeInfinity => a
      case (aLog10Pr, bLog10Pr) if bLog10Pr < aLog10Pr => add(b, a)
      case (aLog10Pr, bLog10Pr) =>
        // sum in log 10 space and convert back to phred.  NB: could try to use log1p for precision.
        -10.0 * (bLog10Pr + Math.log10(1.0 + Math.pow(10.0, aLog10Pr - bLog10Pr)))
    }
  }
  private def sub(a: Double, b: Double): Double = {
    // from phred to log10(pr)
    (a / -10.0, b / -10.0) match {
      case (aLog10Pr, bLog10Pr) if bLog10Pr == Double.NegativeInfinity => a
      case (aLog10Pr, bLog10Pr) if bLog10Pr > aLog10Pr => throw new IllegalArgumentException("Subtraction will be less than zero.")
      case (aLog10Pr, bLog10Pr) => // bLog10Pr <= aLog10Pr
        // subtract in log 10 space and convert back to phred.  NB: could try to use log1p for precision.
      -10.0 * (aLog10Pr + Math.log10(1.0 - Math.pow(10.0, bLog10Pr - aLog10Pr)))
    }
  }
  override def toString: String = this.toInt.toString
  def toInt: Int = {
    // check that we have not lost some precision
    val upper = (this.value + Precision).toInt
    val current = this.value.toInt
    if (upper != current) upper else current
  }
  def toByte: Byte = this.toInt.toByte
}