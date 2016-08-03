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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.NumericTypes.PhredScore

class ConsensusCallerTest extends UnitSpec {
  "ConsensusCaller.adjustBaseQualities" should "cap base qualities" in {
    val caller   = new ConsensusCaller(maxRawBaseQuality=10.toByte, rawBaseQualityShift=0.toByte, errorRatePostLabeling=PhredScore.MaxValue, errorRatePreLabeling=PhredScore.MaxValue)
    val adjusted = Array[Byte](20, 15, 10, 5).map(caller.adjustedErrorProbability).map(PhredScore.fromLogProbability).toSeq
    adjusted shouldBe Seq(10, 10, 10, 5)
  }

  it should "shift base qualities" in {
    val caller   = new ConsensusCaller(maxRawBaseQuality=PhredScore.MaxValue, rawBaseQualityShift=10.toByte, errorRatePostLabeling=PhredScore.MaxValue, errorRatePreLabeling=PhredScore.MaxValue)
    val adjusted = Array[Byte](20, 15, 10, 5).map(caller.adjustedErrorProbability).map(PhredScore.fromLogProbability).toSeq
    adjusted shouldBe Seq(10, 5, 2, 2)
  }

  it should "scale base qualities using the post-umi error rate" in {
    val caller   = new ConsensusCaller(maxRawBaseQuality=PhredScore.MaxValue, rawBaseQualityShift=0.toByte, errorRatePostLabeling=10.toByte, errorRatePreLabeling=PhredScore.MaxValue)
    val adjusted = Array[Byte](20, 15, 10, 5).map(caller.adjustedErrorProbability).map(PhredScore.fromLogProbability).toSeq
    adjusted shouldBe Seq(9, 8, 7, 4)
  }

  it should "cap, shift, and scale base qualities" in {
    val caller   = new ConsensusCaller(maxRawBaseQuality=10.toByte, rawBaseQualityShift=5.toByte, errorRatePostLabeling=10.toByte, errorRatePreLabeling=PhredScore.MaxValue)
    val adjusted = Array[Byte](20, 15, 10, 5).map(caller.adjustedErrorProbability).map(PhredScore.fromLogProbability).toSeq
    adjusted shouldBe Seq(7, 7, 4, 1)
  }

  it should "calculate consensus base and quality given a single base pileup" in {
    val caller = new ConsensusCaller(errorRatePreLabeling=PhredScore.MaxValue, errorRatePostLabeling=PhredScore.MaxValue)
    val builder = caller.builder()
    builder.add('A'.toByte, 20.toByte)
    val (base, qual) = builder.call()
    base shouldBe 'A'
    qual shouldBe 20

    val lls = builder.logLikelihoods
    lls(0) should be > lls(1)
    lls(0) should be > lls(2)
    lls(0) should be > lls(3)
  }

  it should "produce a no-call if two bases are of equal likelihood" in {
    val caller = new ConsensusCaller(errorRatePreLabeling=PhredScore.MaxValue, errorRatePostLabeling=PhredScore.MaxValue)
    val builder = caller.builder()
    builder.call() shouldBe ('N'.toByte, 2.toByte)

    builder.add('A'.toByte, 20.toByte)
    builder.add('C'.toByte, 20.toByte)
    builder.call() shouldBe ('N'.toByte, 2.toByte)
  }

  it should "calculate consensus base and quality, and observation counts, given a massive pileup" in {
    val caller = new ConsensusCaller(errorRatePreLabeling=50.toByte, errorRatePostLabeling=50.toByte)
    val builder = caller.builder()
    (0 to 999).foreach(i => builder.add('C'.toByte, 20.toByte))
    builder.call() shouldBe ('C', 50)
    builder.contributions shouldBe 1000
    builder.observations('A'.toByte) shouldBe 0
    builder.observations('C'.toByte) shouldBe 1000
    builder.observations('G'.toByte) shouldBe 0
    builder.observations('T'.toByte) shouldBe 0

    (0 to 9).foreach(i => builder.add('T'.toByte, 20.toByte))
    builder.call() shouldBe ('C', 50)
    builder.contributions shouldBe 1010
    builder.observations('A'.toByte) shouldBe 0
    builder.observations('C'.toByte) shouldBe 1000
    builder.observations('G'.toByte) shouldBe 0
    builder.observations('T'.toByte) shouldBe 10
    an[Exception] should be thrownBy builder.observations('Z'.toByte)
  }

  it should "calculate consensus base and quality given conflicting evidence" in {
    val caller = new ConsensusCaller(errorRatePreLabeling=50.toByte, errorRatePostLabeling=50.toByte)
    val builder = caller.builder()
    builder.add('A'.toByte, 30.toByte)
    builder.add('C'.toByte, 28.toByte)
    val (base, qual) = builder.call()
    base shouldBe 'A'
    qual.toInt should be <= 5
  }

  it should "support calling multiple pileups from the same builder" in {
    val caller = new ConsensusCaller(errorRatePreLabeling=50.toByte, errorRatePostLabeling=50.toByte)
    val builder = caller.builder()
    builder.add('A'.toByte, 20.toByte)
    builder.call shouldBe ('A'.toByte, 20.toByte)
    builder.contributions shouldBe 1

    builder.reset()
    builder.add('C'.toByte, 20.toByte)
    builder.call shouldBe ('C'.toByte, 20.toByte)
    builder.contributions shouldBe 1
  }
}
