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

/**
  * Some simple utility methods for various mathematical operations that are implemented
  * in efficient albeit non-idiomatic scala.
  */
object MathUtil {

  /** The default number of ulps within which two values are treated as equal.
    *
    * Comparing in ulps scales with the magnitude of the values, unlike an absolute tolerance. Four ulps are chosen
    * to accommodate roundoff between mathematically tied likelihoods in [[com.fulcrumgenomics.umi.ConsensusCaller]].
    * That caller uses Kahan-compensated summation of same-sign (non-positive) log-probability terms: tied bases
    * can accumulate the same terms in different orders. Allowing roughly two ulps of error per sum gives roughly
    * four ulps between them. This is a practical allowance, not a universal rigorous error bound; the rationale
    * depends on retaining Kahan summation and does not apply to arbitrary sums, especially with cancellation.
    *
    * ULP distance counts representable steps, treating signed zeros as equal. Near zero this includes subnormal
    * values on either side, but there is no additional absolute tolerance. Equal infinities match only each other;
    * NaNs never match.
    */
  val maxUlps: Int = 4

  /** The number of representable steps between `a` and `b`, capped at `Long.MaxValue`.
    *
    * Signed zeros have distance zero. NaNs and unequal non-finite values are reported as maximally distant.
    */
  private def ulpsBetween(a: Double, b: Double): Long = {
    if (a == b) 0L
    else if (java.lang.Double.isNaN(a) || java.lang.Double.isNaN(b)) Long.MaxValue
    else if (a.isInfinite || b.isInfinite) Long.MaxValue
    else {
      val aBits = java.lang.Double.doubleToLongBits(a)
      val bBits = java.lang.Double.doubleToLongBits(b)
      // Same-sign finite representations lie in one half of the long range, so their difference cannot overflow.
      if ((aBits < 0) == (bBits < 0)) Math.abs(aBits - bBits)
      else {
        // Across zero, add distances to zero rather than subtracting signed representations, which can overflow.
        val aMagnitude = aBits & Long.MaxValue
        val bMagnitude = bBits & Long.MaxValue
        if (aMagnitude > Long.MaxValue - bMagnitude) Long.MaxValue else aMagnitude + bMagnitude
      }
    }
  }

  /** True if more than one eligible entry of `ds` is approximately equal to `target`.
    *
    * `NaN`s are always ineligible, and `-Infinity` is ineligible unless `allowNegativeInfinity` is set, mirroring
    * the entries that were candidates for the extremum in the first place. The scan stops as soon as a second
    * match is seen; these functions sit on the consensus calling hot path, so neither the predicate nor the count
    * is worth materialising.
    */
  private def hasApproximateTie(ds: Array[Double],
                                target: Double,
                                maxUlps: Int,
                                allowNegativeInfinity: Boolean): Boolean = {
    var seen = false
    var idx  = 0
    while (idx < ds.length) {
      val v = ds(idx)
      if (!java.lang.Double.isNaN(v) && (Double.NegativeInfinity != v || allowNegativeInfinity) &&
          ulpsBetween(v, target) <= maxUlps) {
        if (seen) return true
        seen = true
      }
      idx += 1
    }
    false
  }

  /** Calculates the arithmetic mean of an array of bytes. The computation is performed
    * in integer space and the result will therefore always round down.
    *
    * @param bytes the array of bytes to average
    * @return the floor of the arithmetic mean
    */
  def mean(bytes: Array[Byte]): Byte = {
    if (bytes.length == 0) throw new NoSuchElementException("Cannot compute the mean of a zero length array.")
    var sum: Int = 0
    var i: Int = 0
    while (i < bytes.length) {
      sum += bytes(i)
      i += 1
    }

    (sum / bytes.length).toByte
  }

  /**
    * Finds the minimum element in the array and returns it with its index. By default
    * negative infinity values are excluded - meaning that minimum non-negative-infinity
    * value will be returned.
    *
    * If multiple slots in the array contain the same, minimum, value the index of the
    * first one will be returned.
    *
    * @param ds an array of doubles of length >= 1
    * @param allowNegativeInfinity if true allow the result to be Double.NegativeInfinity
    *                              otherwise return the lowest non-infinite value
    * @param requireUniqueMinimum When false the earliest index at which the minimum value occurs is
    *                             reported. When true, if there are multiple indices with the minimum
    *                             value within `maxUlps` then -1 will be returned for the index.
    * @param maxUlps The inclusive ULP distance from the final minimum within which values tie when uniqueness
    *                is required. Defaults to [[MathUtil.maxUlps]]; zero requires exact equality, including signed
    *                zeros and equal infinities. Must be non-negative, even when uniqueness is not required.
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN and possibly NegativeInfinity)
    * @throws java.lang.IllegalArgumentException if `maxUlps` is negative
    */
  def minWithIndex(ds: Array[Double], allowNegativeInfinity: Boolean=false, requireUniqueMinimum: Boolean=false, maxUlps: Int = MathUtil.maxUlps): (Double, Int) = {
    require(maxUlps >= 0, s"maxUlps must be non-negative, found $maxUlps")
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the min of a zero length array.")
    var min      = Double.MaxValue
    var minIndex = -1
    var assigned = false
    val len      = ds.length
    var idx      = 0
    while (idx < len) {
      val v = ds(idx)
      if (!java.lang.Double.isNaN(v) && (Double.NegativeInfinity != v || allowNegativeInfinity) &&
          (!assigned || v < min)) {
        min = v
        minIndex = idx
        assigned = true
      }

      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("All values are NaNs or negative infinity.")

    // Ties are detected against the final minimum, not a running one: a near-tie is a property of the values, and
    // comparing against a running minimum would only detect it when the tied value happens to appear later.
    if (requireUniqueMinimum && hasApproximateTie(ds, min, maxUlps, allowNegativeInfinity)) (min, -1)
    else (min, minIndex)
  }


  /**
    * Finds the maximum element in the array and returns it with its index.
    *
    * If multiple slots in the array contain the same, maximum, value the index of the
    * first one will be returned.
    *
    * @param ds an array of doubles of length >= 1
    * @param requireUniqueMaximum When false the earliest index at which the maximum value occurs is
    *                             reported. When true, if there are multiple indices with the maximum
    *                             value within `maxUlps` then -1 will be returned for the index.
    * @param maxUlps The inclusive ULP distance from the final maximum within which values tie when uniqueness
    *                is required. Defaults to [[MathUtil.maxUlps]]; zero requires exact equality, including signed
    *                zeros and equal infinities. Must be non-negative, even when uniqueness is not required.
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN)
    * @throws java.lang.IllegalArgumentException if `maxUlps` is negative
    */
  def maxWithIndex(ds: Array[Double], requireUniqueMaximum: Boolean=false, maxUlps: Int = MathUtil.maxUlps): (Double,Int) = {
    require(maxUlps >= 0, s"maxUlps must be non-negative, found $maxUlps")
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the max of a zero length array.")
    var max      =  Double.MinValue
    var maxIndex = -1
    var assigned = false
    val len = ds.length
    var idx = 0
    while (idx < len) {
      val v = ds(idx)
      if (!java.lang.Double.isNaN(v) && (!assigned || v > max)) {
        max = v
        maxIndex = idx
        assigned = true
      }
      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("Array contained only NaNs.")

    // Ties are detected against the final maximum, not a running one: a near-tie is a property of the values, and
    // comparing against a running maximum would only detect it when the tied value happens to appear later.
    // Negative infinities are eligible here: unlike the minimum, they are never excluded from the maximum.
    if (requireUniqueMaximum && hasApproximateTie(ds, max, maxUlps, allowNegativeInfinity=true)) (max, -1)
    else (max, maxIndex)
  }
}
