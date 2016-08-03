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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import MathUtil._

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
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double]())
    an[NoSuchElementException] should be thrownBy minWithIndex(Array[Double](Double.NaN))
  }

  it should "report -1 as the index when requireUniqueMaximum is true" in {
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0)) shouldBe (20.0, 1)
    maxWithIndex(Array(1.0, 20.0, 20.0, 20.0, 5.0), requireUniqueMaximum=true) shouldBe (20.0, -1)
  }
}
