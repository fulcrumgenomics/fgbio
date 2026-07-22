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

  it should "use epsilon for equality comparison when requireUniqueMinimum is true" in {
    // Use a custom epsilon that's large enough to avoid floating point representation issues
    val testEpsilon = 0.01
    val smallDiff = testEpsilon / 2  // 0.005 - within epsilon
    val largeDiff = testEpsilon * 2  // 0.02 - outside epsilon

    // Values within epsilon should be considered equal (returns -1)
    minWithIndex(Array(1.0, 2.0, 1.0 + smallDiff), requireUniqueMinimum=true, epsilon=testEpsilon) shouldBe (1.0, -1)

    // Values outside epsilon should be considered distinct (returns first index)
    minWithIndex(Array(1.0, 2.0, 1.0 + largeDiff), requireUniqueMinimum=true, epsilon=testEpsilon) shouldBe (1.0, 0)

    // Custom epsilon: large epsilon treats different values as equal
    minWithIndex(Array(1.0, 2.0, 1.5), requireUniqueMinimum=true, epsilon=1.0) shouldBe (1.0, -1)

    // Custom epsilon: zero epsilon treats any different values as distinct
    minWithIndex(Array(1.0, 2.0, 1.0 + smallDiff), requireUniqueMinimum=true, epsilon=0.0) shouldBe (1.0, 0)
  }

  "MathUtil.maxWithIndex" should "find the maximum value" in {
    maxWithIndex(Array(5.0, 10.0, -5.0, 100.0)) shouldBe (100, 3)
    maxWithIndex(Array(5.0, 100, -5.0, 100.0)) shouldBe (100, 1)
    maxWithIndex(Array(Double.MinValue)) shouldBe (Double.MinValue, 0)
    maxWithIndex(Array(Double.NaN, Double.MinValue)) shouldBe (Double.MinValue, 1)
  }

  it should "throw exceptions on invalid inputs" in {
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double]())
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double](Double.NaN))
  }

  it should "report -1 as the index when requireUniqueMaximum is true" in {
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0)) shouldBe (20.0, 1)
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0), requireUniqueMaximum=true) shouldBe (20.0, -1)
  }

  it should "use epsilon for equality comparison when requireUniqueMaximum is true" in {
    // Use a custom epsilon that's large enough to avoid floating point representation issues
    val testEpsilon = 0.01
    val smallDiff = testEpsilon / 2  // 0.005 - within epsilon
    val largeDiff = testEpsilon * 2  // 0.02 - outside epsilon

    // Values within epsilon should be considered equal (returns -1)
    maxWithIndex(Array(1.0, 20.0, 20.0 - smallDiff), requireUniqueMaximum=true, epsilon=testEpsilon) shouldBe (20.0, -1)

    // Values outside epsilon should be considered distinct (returns first index)
    maxWithIndex(Array(1.0, 20.0, 20.0 - largeDiff), requireUniqueMaximum=true, epsilon=testEpsilon) shouldBe (20.0, 1)

    // Custom epsilon: large epsilon treats different values as equal
    maxWithIndex(Array(1.0, 20.0, 19.5), requireUniqueMaximum=true, epsilon=1.0) shouldBe (20.0, -1)

    // Custom epsilon: zero epsilon treats any different values as distinct
    maxWithIndex(Array(1.0, 20.0, 20.0 - smallDiff), requireUniqueMaximum=true, epsilon=0.0) shouldBe (20.0, 1)
  }

  it should "detect a near-tie regardless of the order the values appear in" in {
    // A near-tie is a property of the values, not of the order they happen to be stored in. Comparing each value
    // against a *running* maximum only detects the tie when the tied value appears after the maximum.
    val value  = 20.0
    val oneUlp = Math.ulp(value)

    maxWithIndex(Array(1.0, value, value - oneUlp), requireUniqueMaximum=true, epsilon=oneUlp) shouldBe (value, -1)
    maxWithIndex(Array(1.0, value - oneUlp, value), requireUniqueMaximum=true, epsilon=oneUlp) shouldBe (value, -1)
  }

  it should "treat values a single ulp apart as tied at the magnitudes log-likelihoods occupy" in {
    // The default epsilon is the ulp of 1.0. Consensus log-likelihoods are sums of log-probabilities and routinely
    // sit near -500, where a single ulp is ~5.7e-14 -- roughly 250x larger than that epsilon -- so no two distinct
    // values there could ever compare equal. A one-ulp separation is summation noise, not evidence.
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

  it should "detect a near-tie regardless of order when finding the minimum" in {
    val value  = 20.0
    val oneUlp = Math.ulp(value)

    minWithIndex(Array(100.0, value, value + oneUlp), requireUniqueMinimum=true, epsilon=oneUlp) shouldBe (value, -1)
    minWithIndex(Array(100.0, value + oneUlp, value), requireUniqueMinimum=true, epsilon=oneUlp) shouldBe (value, -1)
  }

  it should "treat values a single ulp apart as tied when finding the minimum at large magnitudes" in {
    Seq(1.0, 10.0, 500.0, 5000.0).foreach { value =>
      val oneUlp = Math.ulp(value)
      withClue(s"value=$value ulp=$oneUlp: ") {
        minWithIndex(Array(value, value + oneUlp, 1e9, 1e9), requireUniqueMinimum=true) shouldBe (value, -1)
      }
    }
  }
}
