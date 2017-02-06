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

package com.fulcrumgenomics.math

import com.fulcrumgenomics.math.FishersExactTest.Alternative._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.NumericTypes.LogProbability

import scala.util.Random

class FishersExactTestTest extends UnitSpec {
  "FishersExactTest" should "calculate values the same as R within specified precision" in {
    Seq(
      (0,       0,       0,       0,       1.0),
      (100000,  100000,  100000,  100000,  1.0),
      (1,       2,       3,       4,       1.0),
      (0,       0,       100000,  100000,  1.0),
      (100000,  100000,  100000,  0,       0.0),          //below  R's  or  Java's  precision
      (200000,  100000,  1,       2,       0.259263),
      (100,     100,     100,     0,       3.730187e-23),
      (13736,   9047,    41,      1433,    0.0),          //below  R's  or  Java's  precision
      (66,      14,      64,      4,       0.0424333),
      (351169,  306836,  153739,  2379,    0.0),          //below  R's  or  Java's  precision
      (116449,  131216,  289,     16957,   0.0),          //below  R's  or  Java's  precision
      (137,     159,     9,       23,      0.06088506),
      (129,     90,      21,      20,      0.3919603),
      (14054,   9160,    16,      7827,    0.0),          //below  R's  or  Java's  precision
      (32803,   9184,    32117,   3283,    0.0),          //below  R's  or  Java's  precision
      (2068,    6796,    1133,    0,       0.0),          //below  R's  or  Java's  precision
      (100,     400,     30,      0,       2.989967e-20)

    ).foreach { case (s1, f1, s2, f2, pvalue) =>
      val precision = if (pvalue == 0 || pvalue == 1) 1e-10 else pvalue * 0.01
      FishersExactTest.pValue(s1, f1, s2, f2, TwoSided) shouldBe pvalue +- precision
    }
  }

  it should "compute the same two-sided pvalue with the data in all four reasonable orientations" in {
    val (a1, a2, b1, b2) = (100, 10, 0, 10)
    val expected  = 1.591789e-09
    val precision = expected / 10000
    FishersExactTest.pValue(a1, a2, b1, b2, TwoSided) shouldBe expected +- precision
    FishersExactTest.pValue(a1, a2, b1, b2, TwoSided) shouldBe FishersExactTest.pValue(b1, b2, a1, a2, TwoSided) +- precision
    FishersExactTest.pValue(a1, a2, b1, b2, TwoSided) shouldBe FishersExactTest.pValue(a1, b1, a2, b2, TwoSided) +- precision
    FishersExactTest.pValue(a1, a2, b1, b2, TwoSided) shouldBe FishersExactTest.pValue(b1, a1, b2, a2, TwoSided) +- precision
  }

  it should "compute the same values as R for alternative=greater" in {
    FishersExactTest.pValue(20, 20, 20, 20, Greater) shouldBe 0.5884  +- 0.0001
    FishersExactTest.pValue( 0, 20, 20, 20, Greater) shouldBe 1.0     +- 0.0001
    FishersExactTest.pValue(20,  0, 20, 20, Greater) shouldBe 3.288e-05 +- 0.01e-05
    FishersExactTest.pValue(20, 20,  0, 20, Greater) shouldBe 3.288e-05 +- 0.01e-05
    FishersExactTest.pValue(20, 20, 20,  0, Greater) shouldBe 1.0     +- 0.0001

  }

  it should "compute the same values as R for alternative=less" in {
    FishersExactTest.pValue(20, 20, 20, 20, Less) shouldBe 0.5884    +- 0.0001
    FishersExactTest.pValue( 0, 20, 20, 20, Less) shouldBe 3.288e-05 +- 0.01e-05
    FishersExactTest.pValue(20,  0, 20, 20, Less) shouldBe 1.0       +- 0.0001
    FishersExactTest.pValue(20, 20,  0, 20, Less) shouldBe 1.0       +- 0.0001
    FishersExactTest.pValue(20, 20, 20,  0, Less) shouldBe 3.288e-05 +- 0.01e-05
  }
}
