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

import java.lang.Math.{log, exp, max, min}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.util.MathUtil
import org.apache.commons.math3.distribution.HypergeometricDistribution

/**
  * Implementation of Fisher's Exact Test for 2x2 contingency tables. Follow's the
  * implementation from R's "fisher.test".
  */
object FishersExactTest {
  /** The alternatives for p-value computation. */
  object Alternative extends Enumeration {
    val TwoSided, Greater, Less = Value
    type Alternative = Alternative.Value
  }

  import Alternative._

  /**
    * Computes the 2-sided p-value of the Fisher's Exact Test on 2x2 contingency table.
    * @param succ1 The number of successes in the test condition
    * @param fail1 The number of failures in the test condition
    * @param succ2 The number of successes under the null hypothesis
    * @param fail2 The number of failures under the null hypothesis
    * @param alternative either `TwoSided` to test for any difference in proportion,
    *                           `Less` to test if the test condition had a lower proportion of successes than the null, or
    *                           `Greater` to test if the test condition had a higher proportion of successes than the null
    */
  def pValue(succ1: Int, fail1: Int, succ2: Int, fail2: Int, alternative: Alternative): Double = {
    val m = succ1 + fail1
    val n = succ2 + fail2
    val k = succ1 + succ2
    val x = succ1
    val lo = max(0, k - n)
    val hi = min(k, m)
    val support = Range.inclusive(lo, hi)

    // special case, support has only one value
    if (lo == hi) {
      1.0
    }
    else {
      val dist = new HypergeometricDistribution(null, m + n, m, k)

      // Floating point imprecision may yield values slightly > 1, hence all the min(x,1)s below
      alternative match {
        case TwoSided =>
          val logds     = support.toArray.map(dist.logProbability)
          val threshold = logds(succ1 - lo)
          min(sumOfLogs(logds.filter(_ <= threshold)), 1.0)
        case Less =>
          val logps = Range.inclusive(0, x).toArray.map(dist.logProbability)
          min(sumOfLogs(logps), 1.0)
        case Greater =>
          val logps = Range.inclusive(0, x-1).toArray.map(dist.logProbability)
          1.0 - min(sumOfLogs(logps), 1.0)
      }
    }
  }

  /** Calculates the non-log sum of the log values. */
  private def sumOfLogs(ds: Array[Double]): Double = if (ds.length == 0) 0 else exp(logSumOfLogs(ds))

  /** Calculates the sum of log values. (i.e log(sum(exp(val))). */
  private[math] def logSumOfLogs(logVals: Array[Double]): Double = {
    val (maxValue, maxIndex) = MathUtil.maxWithIndex(logVals)

    if (maxValue.isNegInfinity) {
      maxValue
    }
    else {
      var sum = 1.0
      forloop (from=0, until=logVals.length) { i =>
        if (i != maxIndex && logVals(i) != Double.NegativeInfinity) sum += exp(logVals(i) - maxValue)
      }

      require(!sum.isNaN && !sum.isPosInfinity, "log values must be non-infinite and non-NAN")
      maxValue + (if(sum != 1.0) log(sum) else 0.0)
    }
  }
}
