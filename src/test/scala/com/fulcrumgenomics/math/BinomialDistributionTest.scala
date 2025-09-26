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

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.NumericTypes.LogProbability
import org.apache.commons.math3.distribution.{BinomialDistribution => CommonsBinomial}

class BinomialDistributionTest extends UnitSpec {
  val binom = new BinomialDistribution

  "BinomialDistribution.coefficient" should "calculate the correct number of combinations given k=0 or k=n" in {
    binom.coefficient(100, 0)   shouldBe 1
    binom.coefficient(100, 100) shouldBe 1
  }

  it should "fail if n or k is less than 0" in {
    an[Exception] shouldBe thrownBy { binom.coefficient(100, -2) }
    an[Exception] shouldBe thrownBy { binom.coefficient(-2,  0) }
  }

  it should "fail if k > n" in {
    an[Exception] shouldBe thrownBy { binom.coefficient(n=100, k=101) }
  }

  Seq((20, 5, BigInt(15504)), (75, 25, BigInt("52588547141148893628"))).foreach { case (n, k, coeff) =>
    it should s"correctly calculate coefficient(n=$n, k=$k) = $coeff" in {
      binom.coefficient(n=n, k=k) shouldBe coeff
    }
  }

  it should "work for small values of n and k" in {
    binom.coefficient(20, 5) shouldBe BigInt("15504")
  }

  it should "work for moderate values of n and k" in {
    binom.coefficient(75, 25) shouldBe BigInt("52588547141148893628")
  }

  it should "work for very large values" in {
    binom.coefficient(200, 100) shouldBe BigInt("90548514656103281165404177077484163874504589675413336841320")
  }

  "BinomialCoefficient.cumulativeProbability" should "match commons math for small numbers" in {
    val commons = new CommonsBinomial(20, 0.05)
    Range.inclusive(0, 20).foreach { k =>
      binom.cumulativeProbability(k, n=20, p=0.05).toDouble shouldBe commons.cumulativeProbability(k) +- 1e-10
    }
  }

  it should "compute values where commons math runs out of precision" in {
    val (k, n, p) = (39, 55, 3/150.0)

    val commons = new CommonsBinomial(n, p)

    // Calculate the cumulative using log probabilities
    val range    = Range(0,k)
    val probs    = range.map(i => commons.logProbability(i)).toArray
    val sum      = LogProbability.or(probs)
    val oneMinus = LogProbability.not(sum)
    val linear   = LogProbability.expProb(oneMinus)

    // Results using three methods
    val commonsResult    = 1 - commons.cumulativeProbability(k-1)
    val commonsLogResult = linear
    val result           = binom.cumulativeProbability(k=k-1, n=n, p=p, lower=false).toDouble

    commonsResult    shouldBe 0.0 // The real answer is not 0.0, here to show this method underflows
    commonsLogResult shouldBe 0.0 // The real answer is not 0.0, here to show this method underflows
    result           shouldBe 1.193E-53 +- 0.001E-53 // Calculated @ https://www.wolframalpha.com/input/?i=cumulative+binomial+probability
  }

  it should "fail if given probabilities < 0 or > 1" in {
    an[Exception] shouldBe thrownBy { binom.cumulativeProbability(k=5, n=10, p= -0.01) }
    an[Exception] shouldBe thrownBy { binom.cumulativeProbability(k=5, n=10, p= 1.00000000000001) }
  }

  it should "yield answers from lower=true and lower=false than sum to 1" in {
    val (k, n, p) = (5, 10, 0.3)
    val lower = binom.cumulativeProbability(k, n, p, lower=true)
    val upper = binom.cumulativeProbability(k, n, p, lower=false)
    (lower + upper).toDouble shouldBe 1.0
  }
}
