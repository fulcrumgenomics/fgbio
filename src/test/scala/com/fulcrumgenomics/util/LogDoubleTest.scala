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
import com.fulcrumgenomics.util.PhredScore._
import com.fulcrumgenomics.util.LogDouble._
import scala.math._

/**
  * Tests for PhredValue.
  */
class LogDoubleTest extends UnitSpec {
  "PhredScore" should "convert probabilities to phred scores" in {
    toPhredScore(0.0.toLogDouble) shouldBe Double.PositiveInfinity
    toPhredScore(0.1.toLogDouble) shouldBe 10.0 +- Precision
    toPhredScore(0.5.toLogDouble) shouldBe 3.0103 +- Precision
    toPhredScore(1.0.toLogDouble) shouldBe 0.0
  }

  it should "convert phred scores to probabilities" in {
    Double.PositiveInfinity.fromPhredScore.linearValue shouldBe 0.0 +- Precision
    10.0.fromPhredScore.linearValue shouldBe 0.1 +- Precision
    3.0103.fromPhredScore.linearValue shouldBe 0.5 +- Precision
    0.0.fromPhredScore.linearValue shouldBe 1.0
  }

  it should "display phred scores as integers" in {
    10.0.fromPhredScore.toPhredScoreInt shouldBe 10
    100.fromPhredScore.toPhredScoreInt shouldBe 100
    (10-Precision).fromPhredScore.toPhredScoreInt shouldBe 10
    (10-(Precision*10)).fromPhredScore.toPhredScoreInt shouldBe 9
    9.0001.fromPhredScore.toPhredScoreInt shouldBe 9

    10.0.fromPhredScore.toPhredScoreString shouldBe "10"
    100.0.fromPhredScore.toPhredScoreString shouldBe "100"
    (10-Precision).fromPhredScore.toPhredScoreString shouldBe "10"
    (10-(Precision*10)).fromPhredScore.toPhredScoreString shouldBe "9"
    (9.0001).fromPhredScore.toPhredScoreString shouldBe "9"
  }

  "LogDouble" should "convert doubles to log-space doubles" in {
    toLogDouble(0).logValue shouldBe Double.NegativeInfinity
    toLogDouble(0.1).logValue shouldBe -2.302585 +- 0.00001
    toLogDouble(0.5).logValue shouldBe -0.6931472 +- 0.00001
    toLogDouble(1.0).logValue shouldBe 0.0
  }

  it should "throw an exception if a negative value is to be converted to a log-space doubles" in {
    an[Exception] should be thrownBy toLogDouble(-1)
  }

  it should "convert log-space doubles to doubles" in {
    toLogDouble(0.0).linearValue shouldBe 0.0 +- 0.00001
    toLogDouble(0.1).linearValue shouldBe 0.1 +- 0.00001
    toLogDouble(0.5).linearValue shouldBe 0.5 +- 0.00001
    toLogDouble(1.0).linearValue shouldBe 1.0 +- 0.00001
    toLogDouble(100.0).linearValue shouldBe 100.0 +- 0.00001
  }

  it should "divide doubles in log space" in {
    (LogDouble(10.0) / LogDouble(10)).logValue shouldBe LogDouble(0.0).logValue +- 0.00001
    (LogDouble(Double.NegativeInfinity) / LogDouble(10)).logValue shouldBe Zero.logValue +- 0.00001
    (toLogDouble(100.0) / toLogDouble(10.0)).linearValue shouldBe toLogDouble(10.0).linearValue +- 0.00001
    (toLogDouble(0.1) / toLogDouble(5.0)).linearValue shouldBe toLogDouble(0.02).linearValue +- 0.00001
    an[IllegalStateException] should be thrownBy toLogDouble(10) / Zero
    (toLogDouble(10) / 10.0.toLogDouble).logValue shouldBe 0.0
  }

  it should "multiply doubles in log space" in {
    (LogDouble(10) *  LogDouble(10)).logValue shouldBe  LogDouble(20.0).logValue
    (LogDouble(10) * Zero).logValue shouldBe Zero.logValue
    (Zero * LogDouble(10)).logValue shouldBe Zero.logValue
    (LogDouble(0.0) * LogDouble(10.0)).linearValue shouldBe LogDouble(10.0).linearValue +- 0.00001
    (toLogDouble(10) * 5.0.toLogDouble).linearValue shouldBe toLogDouble(50).linearValue +- 0.00001
    (LogDouble(0.0) * 0.1.toLogDouble).logValue shouldBe toLogDouble(0.1).logValue
    (toLogDouble(0.0) * 1.0.toLogDouble).logValue shouldBe toLogDouble(0.0).logValue
  }

  it should "add doubles in log space" in {
    (LogDouble(10) + LogDouble(10)).linearValue shouldBe exp(10)*2 +- 0.00001
    (LogDouble(10) + LogDouble(20)).linearValue shouldBe exp(10)+exp(20) +- 0.00001
    (LogDouble(20) + LogDouble(10)).linearValue shouldBe exp(20)+exp(10) +- 0.00001
    (LogDouble(10) + Zero).linearValue shouldBe exp(10)
    (Zero + LogDouble(10)).linearValue shouldBe exp(10)
    (LogDouble(10) + 0.0.toLogDouble).linearValue shouldBe exp(10)
  }

  it should "subtract doubles in log space" in {
    (LogDouble(10) - LogDouble(10)).linearValue shouldBe Zero.linearValue
    (10.fromPhredScore - 10.fromPhredScore).linearValue shouldBe Zero.linearValue
    (10.fromPhredScore - 20.fromPhredScore).linearValue shouldBe (0.1-0.01) +- 0.00001
    an[IllegalArgumentException] should be thrownBy (20.fromPhredScore - 10.fromPhredScore).linearValue
    (toLogDouble(10) - Zero).linearValue shouldBe 10.0 +- 0.00001
    (toLogDouble(10) - 0.0.toLogDouble).linearValue shouldBe 10.0 +- 0.00001
    an[IllegalArgumentException] should be thrownBy (toLogDouble(1.0) - toLogDouble(1.0000000000000004))
    an[IllegalArgumentException] should be thrownBy (toLogDouble(1.0) - toLogDouble(-0.0000000000000004))
  }

  it should "1 - doubles in log space" in {
    10.fromPhredScore.oneMinus().linearValue shouldBe 0.9 +- 0.00001
    20.fromPhredScore.oneMinus().linearValue shouldBe 0.99 +- 0.00001
    toLogDouble(0.9).oneMinus().linearValue shouldBe  0.1 +- 0.00001
    toLogDouble(0.99).oneMinus().linearValue shouldBe 0.01 +- 0.00001
    One.logValue shouldBe 0.0
    One.linearValue shouldBe 1.0
    Zero.linearValue shouldBe 0.0
    Zero.logValue shouldBe Double.NegativeInfinity
    Zero.oneMinus().linearValue shouldBe One.linearValue +- 0.00001
    One.oneMinus().linearValue shouldBe Zero.linearValue
    toLogDouble(1.0000000000000004).oneMinus().linearValue shouldBe 0.0 +- 0.00001
    LogDouble(log(0.0000000000000004)).oneMinus().linearValue shouldBe 1.0 +- 0.00001
  }

  it should "compute the mean of doubles in log space" in {
    mean(toLogDouble(0.1)).linearValue shouldBe 0.1  +- 0.00001
    mean(toLogDouble(0.1), toLogDouble(0.1)).linearValue shouldBe 0.1  +- 0.00001
    mean(toLogDouble(0.4), toLogDouble(0.2)).linearValue shouldBe 0.3  +- 0.00001
  }

  it should "compare doubles in log space" in {
    toLogDouble(0.01) should be < toLogDouble(0.1)
    toLogDouble(0.1) should be > toLogDouble(0.01)
    toLogDouble(0.1) should be >= toLogDouble(0.1)
  }

  it should "addAll" in {
    LogDouble.addAll(log(10), log(10), log(10)) shouldBe log(30)
  }
}
