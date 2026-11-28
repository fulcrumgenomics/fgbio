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

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.MathUtil._

/** Tests for the MathUtil class. */
class MathUtilTest extends UnitSpec {
  "MathUtil.mean(Array[Byte])" should "compute averages correctly" in {
    mean(Array[Byte](5)) shouldBe 5
    mean(Array[Byte](5, 5, 5)) shouldBe 5
    mean(Array[Byte](10, -10)) shouldBe 0
    mean(Array[Byte](0)) shouldBe 0
    mean(Array[Byte](1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) shouldBe 1
    mean((0 to 100).map(_.toByte).toArray) shouldBe 50
  }

  it should "throw an exception if given no input" in {
    an[NoSuchElementException] should be thrownBy mean(Array[Byte]())
  }

  "MathUtil.minWithIndex" should "find the minimum value" in {
    minWithIndex(Array(5.0, 10.0, -5.0, 100.0)) shouldBe (-5.0, 2)
    minWithIndex(Array(5.0, 10.0, -5.0, 100.0, Double.NegativeInfinity)) shouldBe (-5.0, 2)
    minWithIndex(Array(5.0, 10.0, -5.0, 10.00, Double.PositiveInfinity, Double.NaN, Double.NegativeInfinity)) shouldBe (-5.0, 2)
    minWithIndex(Array(5.0, Double.NegativeInfinity, -5.0, Double.NaN), allowNegativeInfinity=true) shouldBe (Double.NegativeInfinity, 1)
    minWithIndex(Array(Double.NegativeInfinity, Double.NegativeInfinity, Double.NaN), allowNegativeInfinity=true) shouldBe (Double.NegativeInfinity, 0)
    minWithIndex(Array(Double.NaN, Double.MaxValue), allowNegativeInfinity=true) shouldBe (Double.MaxValue, 1)
  }

  it should "throw exceptions on invalid inputs" in {
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double]())
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double](Double.NegativeInfinity))
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double](Double.NegativeInfinity, Double.NaN))
  }

  it should "report -1 as the index when requireUniqueMinimum is true" in {
    minWithIndex(Array(1.0, 2.0, 3.0, 4.0, 1.0)) shouldBe (1.0, 0)
    minWithIndex(Array(1.0, 2.0, 3.0, 4.0, 1.0), requireUniqueMinimum=true) shouldBe (1.0, -1)
  }

  "MathUtil.maxWithIndex" should "find the maximum value" in {
    maxWithIndex(Array(5.0, 10.0, -5.0, 100.0)) shouldBe (100, 3)
    maxWithIndex(Array(5.0, 100, -5.0, 100.0)) shouldBe (100, 1)
    maxWithIndex(Array(Double.MinValue)) shouldBe (Double.MinValue, 0)
    maxWithIndex(Array(Double.NaN, Double.MinValue)) shouldBe (Double.MinValue, 1)
  }

  it should "throw exceptions on invalid inputs" in {
    an[NoSuchElementException] should be thrownBy maxWithIndex(Array[Double]())
    an[NoSuchElementException] should be thrownBy maxWithIndex(Array[Double](Double.NaN))
  }

  it should "report -1 as the index when requireUniqueMaximum is true" in {
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0)) shouldBe (20.0, 1)
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0), requireUniqueMaximum=true) shouldBe (20.0, -1)
  }

  it should "detect a near-tie regardless of the order the values appear in" in {
    // A near-tie is a property of the values, not of the order they happen to be stored in. Comparing each value
    // against a *running* maximum only detects the tie when the tied value appears after the maximum.
    val value  = 20.0
    val oneUlp = Math.ulp(value)

    maxWithIndex(Array(1.0, value, value - oneUlp), requireUniqueMaximum=true) shouldBe (value, -1)
    maxWithIndex(Array(1.0, value - oneUlp, value), requireUniqueMaximum=true) shouldBe (value, -1)
  }

  it should "treat values a single ulp apart as tied at the magnitudes log-likelihoods occupy" in {
    // The former absolute tolerance was the ulp of 1.0. Near -500 a single ulp is ~5.7e-14, roughly 250x larger,
    // so that tolerance could not match distinct likelihoods. ULP comparison accommodates summation roundoff.
    Seq(-1.0, -10.0, -500.0, -5000.0).foreach { value =>
      val oneUlp = Math.ulp(value)
      withClue(s"value=$value ulp=$oneUlp: ") {
        maxWithIndex(Array(value, value - oneUlp, -1e9, -1e9), requireUniqueMaximum=true) shouldBe (value, -1)
      }
    }
  }

  it should "still report a unique maximum when values are clearly separated at those magnitudes" in {
    Seq(-1.0, -10.0, -500.0, -5000.0).foreach { value =>
      withClue(s"value=$value: ") {
        maxWithIndex(Array(value, value - 1.0, -1e9, -1e9), requireUniqueMaximum=true) shouldBe (value, 0)
      }
    }
  }

  it should "honour an explicit maxUlps" in {
    val value  = -500.0
    val oneUlp = Math.ulp(value)

    // No ulps of tolerance: only exactly equal values tie.
    maxWithIndex(Array(value, value - oneUlp, -1e9, -1e9), requireUniqueMaximum=true, maxUlps=0) shouldBe (value, 0)

    // A wider tolerance ties values several ulps apart.
    maxWithIndex(Array(value, value - (oneUlp * 3), -1e9, -1e9), requireUniqueMaximum=true, maxUlps=8) shouldBe (value, -1)
  }

  it should "report a tie when every eligible value is negative infinity" in {
    // Math.abs(-Inf - -Inf) is NaN, so an absolute comparison alone never matched these; a position with no evidence
    // for any value must not be resolved to the first index.
    maxWithIndex(Array(Double.NegativeInfinity, Double.NegativeInfinity), requireUniqueMaximum=true) shouldBe (Double.NegativeInfinity, -1)
    minWithIndex(Array(Double.NegativeInfinity, Double.NegativeInfinity), allowNegativeInfinity=true, requireUniqueMinimum=true) shouldBe (Double.NegativeInfinity, -1)
  }

  it should "ignore NaNs, and negative infinities when they are not allowed, when detecting ties" in {
    maxWithIndex(Array(Double.NaN, 1.0, Double.NaN), requireUniqueMaximum=true) shouldBe (1.0, 1)
    minWithIndex(Array(Double.NaN, 1.0, Double.NaN), requireUniqueMinimum=true) shouldBe (1.0, 1)
    minWithIndex(Array(Double.NegativeInfinity, 1.0, Double.NegativeInfinity), requireUniqueMinimum=true) shouldBe (1.0, 1)
  }

  it should "reject a negative maxUlps" in {
    an[IllegalArgumentException] should be thrownBy maxWithIndex(Array(1.0, 2.0), requireUniqueMaximum=true, maxUlps = -1)
    an[IllegalArgumentException] should be thrownBy minWithIndex(Array(1.0, 2.0), requireUniqueMinimum=true, maxUlps = -1)
    an[IllegalArgumentException] should be thrownBy maxWithIndex(Array(1.0, 2.0), maxUlps = -1)
    an[IllegalArgumentException] should be thrownBy minWithIndex(Array(1.0, 2.0), maxUlps = -1)
  }

  it should "not treat negative zero as tied with a distant value of the opposite sign" in {
    // -0.0 is not `< 0`, but its bit representation is Long.MinValue. Comparing the sign of the value rather than
    // of the representation would put these in the same branch and overflow the ulp subtraction.
    maxWithIndex(Array(1.0, -0.0), requireUniqueMaximum=true) shouldBe (1.0, 0)
    minWithIndex(Array(-0.0, 1.0), requireUniqueMinimum=true) shouldBe (-0.0, 0)
  }

  it should "detect a near-tie regardless of order when finding the minimum" in {
    val value  = 20.0
    val oneUlp = Math.ulp(value)

    minWithIndex(Array(100.0, value, value + oneUlp), requireUniqueMinimum=true) shouldBe (value, -1)
    minWithIndex(Array(100.0, value + oneUlp, value), requireUniqueMinimum=true) shouldBe (value, -1)
  }

  it should "treat values a single ulp apart as tied when finding the minimum at large magnitudes" in {
    Seq(1.0, 10.0, 500.0, 5000.0).foreach { value =>
      val oneUlp = Math.ulp(value)
      withClue(s"value=$value ulp=$oneUlp: ") {
        minWithIndex(Array(value, value + oneUlp, 1e9, 1e9), requireUniqueMinimum=true) shouldBe (value, -1)
      }
    }
  }

  it should "apply inclusive ULP boundaries to positive and negative extrema in either order" in {
    val anchors = Seq(-5000.0, -500.0, -20.0, -1.0, 1.0, 20.0, 500.0, 5000.0,
      -8.0 * java.lang.Double.MIN_VALUE, 8.0 * java.lang.Double.MIN_VALUE,
      -java.lang.Double.MIN_NORMAL, java.lang.Double.MIN_NORMAL)
    anchors.foreach { value =>
      Seq(0, 1, 4, 5).foreach { distance =>
        val next = Iterator.iterate(value)(v => Math.nextUp(v)).drop(distance).next()
        Seq(Array(value, next), Array(next, value)).foreach { values =>
          withClue(s"values=${values.mkString(",")} distance=$distance: ") {
            val min = values.min
            val max = values.max
            minWithIndex(values) shouldBe (min, values.indexOf(min))
            maxWithIndex(values) shouldBe (max, values.indexOf(max))
            minWithIndex(values, requireUniqueMinimum=true) shouldBe (min, if (distance <= 4) -1 else values.indexOf(min))
            maxWithIndex(values, requireUniqueMaximum=true) shouldBe (max, if (distance <= 4) -1 else values.indexOf(max))
            Seq(0, 1, 4, 5).foreach { tolerance =>
              val minIndex = if (distance <= tolerance) -1 else values.indexOf(min)
              val maxIndex = if (distance <= tolerance) -1 else values.indexOf(max)
              minWithIndex(values, requireUniqueMinimum=true, maxUlps=tolerance) shouldBe (min, minIndex)
              maxWithIndex(values, requireUniqueMaximum=true, maxUlps=tolerance) shouldBe (max, maxIndex)
            }
          }
        }
      }
    }
  }

  it should "compare only to the final extremum rather than chaining approximate ties" in {
    val low    = 20.0
    val middle = Iterator.iterate(low)(v => Math.nextUp(v)).drop(4).next()
    val high   = Iterator.iterate(low)(v => Math.nextUp(v)).drop(8).next()
    Seq(low, middle, high).permutations.foreach { values =>
      minWithIndex(values.toArray, requireUniqueMinimum=true) shouldBe (low, -1)
      maxWithIndex(values.toArray, requireUniqueMaximum=true) shouldBe (high, -1)
    }
    Seq(low, low, high).permutations.foreach { values =>
      maxWithIndex(values.toArray, requireUniqueMaximum=true) shouldBe (high, values.indexOf(high))
    }
    Seq(low, high, high).permutations.foreach { values =>
      minWithIndex(values.toArray, requireUniqueMinimum=true) shouldBe (low, values.indexOf(low))
    }
  }

  it should "treat signed zeros as exactly equal and preserve the first zero without uniqueness" in {
    Seq(Array(-0.0, 0.0), Array(0.0, -0.0)).foreach { values =>
      minWithIndex(values, requireUniqueMinimum=true, maxUlps=0)._2 shouldBe -1
      maxWithIndex(values, requireUniqueMaximum=true, maxUlps=0)._2 shouldBe -1
      val (min, minIndex) = minWithIndex(values)
      val (max, maxIndex) = maxWithIndex(values)
      minIndex shouldBe 0
      maxIndex shouldBe 0
      java.lang.Double.doubleToRawLongBits(min) shouldBe java.lang.Double.doubleToRawLongBits(values(0))
      java.lang.Double.doubleToRawLongBits(max) shouldBe java.lang.Double.doubleToRawLongBits(values(0))
    }
  }

  it should "count ULPs next to and across zero without an absolute tolerance" in {
    val smallest = java.lang.Double.MIN_VALUE
    Seq(-0.0, 0.0).foreach { zero =>
      Seq(-1.0, 1.0).foreach { sign =>
        Seq(1, 4, 5).foreach { distance =>
          val subnormal = sign * distance * smallest
          Seq(Array(zero, subnormal), Array(subnormal, zero)).foreach { values =>
            val min = values.min
            val max = values.max
            Seq(0, 1, 4, 5).foreach { tolerance =>
              minWithIndex(values, requireUniqueMinimum=true, maxUlps=tolerance) shouldBe
                (min, if (distance <= tolerance) -1 else values.indexOf(min))
              maxWithIndex(values, requireUniqueMaximum=true, maxUlps=tolerance) shouldBe
                (max, if (distance <= tolerance) -1 else values.indexOf(max))
            }
          }
        }
      }
    }
    Seq((1, 1), (2, 2), (2, 3)).foreach { case (negativeSteps, positiveSteps) =>
      val low  = -negativeSteps * smallest
      val high = positiveSteps * smallest
      val tied = negativeSteps + positiveSteps <= MathUtil.maxUlps
      Seq(Array(low, high), Array(high, low)).foreach { values =>
        minWithIndex(values, requireUniqueMinimum=true) shouldBe (low, if (tied) -1 else values.indexOf(low))
        maxWithIndex(values, requireUniqueMaximum=true) shouldBe (high, if (tied) -1 else values.indexOf(high))
      }
    }
  }

  it should "keep distant values distinct without overflowing the ULP distance" in {
    Seq(Array(-Double.MaxValue, Double.MaxValue), Array(-1.0, 1.0), Array(0.0, 1e-20)).foreach { values =>
      minWithIndex(values, requireUniqueMinimum=true, maxUlps=Int.MaxValue) shouldBe (values(0), 0)
      maxWithIndex(values, requireUniqueMaximum=true, maxUlps=Int.MaxValue) shouldBe (values(1), 1)
    }
  }

  it should "distinguish finite values from infinities and tie equal positive infinities" in {
    maxWithIndex(Array(Double.MaxValue, Double.PositiveInfinity), requireUniqueMaximum=true) shouldBe (Double.PositiveInfinity, 1)
    minWithIndex(Array(Double.MaxValue, Double.PositiveInfinity), requireUniqueMinimum=true) shouldBe (Double.MaxValue, 0)
    maxWithIndex(Array(Double.NegativeInfinity, -Double.MaxValue), requireUniqueMaximum=true) shouldBe (-Double.MaxValue, 1)
    minWithIndex(Array(Double.NegativeInfinity, -Double.MaxValue), allowNegativeInfinity=true, requireUniqueMinimum=true) shouldBe (Double.NegativeInfinity, 0)
    maxWithIndex(Array(Double.PositiveInfinity, Double.NaN, Double.PositiveInfinity), requireUniqueMaximum=true, maxUlps=0) shouldBe (Double.PositiveInfinity, -1)
    minWithIndex(Array(Double.PositiveInfinity, Double.NaN, Double.PositiveInfinity), requireUniqueMinimum=true, maxUlps=0) shouldBe (Double.PositiveInfinity, -1)
    Seq(Double.NegativeInfinity, -Double.MaxValue, -0.0, 0.0, Double.MaxValue, Double.PositiveInfinity).foreach { value =>
      maxWithIndex(Array(value), requireUniqueMaximum=true) shouldBe (value, 0)
      minWithIndex(Array(value), allowNegativeInfinity=true, requireUniqueMinimum=true) shouldBe (value, 0)
    }
  }
}
