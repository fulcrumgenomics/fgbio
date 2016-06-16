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
import com.fulcrumgenomics.util.PhredValue._

/**
  * Tests for PhredValue.
  */
class PhredValueTest extends UnitSpec {
  "PhredValue" should "convert probabilities to Phreds" in {
    toPhred(0) shouldBe Double.PositiveInfinity
    toPhred(0.1) shouldBe 10.0
    toPhred(0.5) shouldBe 3.0103 +- 0.001
    toPhred(1.0) shouldBe 0.0
  }

  it should "convert Phreds to probabilities" in {
    toPr(Double.PositiveInfinity) shouldBe 0
    toPr(10.0) shouldBe 0.1
    toPr(3.0103) shouldBe 0.5 +- 0.001
    toPr(0.0) shouldBe 1.0
  }

  it should "divide Phreds" in {
    PhredValue(10) / PhredValue(10) shouldBe PhredValue(0.0)
    PhredValue(Double.PositiveInfinity) / PhredValue(10) shouldBe ZeroProbability
    (PhredValue(3.0103) / PhredValue(10.0)).value shouldBe toPhred(5.0) +- 0.001
    (PhredValue(3.0103) / 5.0).value shouldBe 10.0 +- 0.001
    an[IllegalStateException] should be thrownBy PhredValue(10) / ZeroProbability
    (PhredValue(10) / 0.1) shouldBe OneProbability
  }

  it should "multiply Phreds" in {
    PhredValue(10) *  PhredValue(10) shouldBe  PhredValue(20.0)
    PhredValue(10) * ZeroProbability shouldBe ZeroProbability
    ZeroProbability * PhredValue(10) shouldBe ZeroProbability
    PhredValue(0.0) * PhredValue(10.0) shouldBe  PhredValue(10.0)
    (PhredValue(10) * 5.0).value shouldBe 3.0103 +- 0.001
    PhredValue(0.0) * 0.1 shouldBe PhredValue(10.0)
    PhredValue(0.0) * 1.0 shouldBe PhredValue(0.0)
  }

  it should "add probabilities" in {
    (PhredValue(10) + PhredValue(10)).value shouldBe 6.9897 +- 0.001
    (PhredValue(10) + PhredValue(20)).value shouldBe 9.5860 +- 0.001
    (PhredValue(20) + PhredValue(10)).value shouldBe 9.5860 +- 0.001
    (PhredValue(10) + ZeroProbability).value shouldBe 10
    (ZeroProbability + PhredValue(10)).value shouldBe 10
    (PhredValue(10) + 0.0).value shouldBe 10
  }

  it should "subtract Phreds" in {
    (PhredValue(10) - PhredValue(10)).value shouldBe ZeroProbability.value
    (PhredValue(10) - PhredValue(20)).value shouldBe toPhred(0.1-0.01)
    an[IllegalArgumentException] should be thrownBy (PhredValue(20) - PhredValue(10)).value
    (PhredValue(10) - ZeroProbability).value shouldBe 10
    (PhredValue(10) - 0.0).value shouldBe 10
  }

  it should "invert Phreds (1-pr)" in {
    PhredValue(10).inv().value shouldBe 0.4575749 +- 0.001
    PhredValue(20).inv().value shouldBe 0.0436 +- 0.001
    PhredValue( 0.4575749).inv().value shouldBe 10.0 +- 0.001
    PhredValue(0.04364805).inv().value shouldBe 20.0 +- 0.001
    OneProbability.value shouldBe 0.0
    ZeroProbability.value shouldBe Double.PositiveInfinity
    ZeroProbability.inv().value shouldBe 0.0 +- 0.001
    OneProbability.inv().value shouldBe Double.PositiveInfinity
  }

  it should "display Phreds as integers" in {
    PhredValue(10).toInt shouldBe 10
    PhredValue(100).toInt shouldBe 100
    PhredValue(10-Precision).toInt shouldBe 10
    PhredValue(10-(Precision*10)).toInt shouldBe 9
    PhredValue(9.0001).toInt shouldBe 9
    
    PhredValue(10).toString shouldBe "10"
    PhredValue(100).toString shouldBe "100"
    PhredValue(10-Precision).toString shouldBe "10"
    PhredValue(10-(Precision*10)).toString shouldBe "9"
    PhredValue(9.0001).toString shouldBe "9"
  }
}
