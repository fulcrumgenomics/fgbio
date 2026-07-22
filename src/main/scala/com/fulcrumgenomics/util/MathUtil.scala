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

  /** The difference between 1.0 and the next larger representable number.
    *
    * @see https://doc.rust-lang.org/std/primitive.f64.html#associatedconstant.EPSILON
    */
  val epsilon: Double = 1.0 / Math.pow(2, 52)

  /** The default number of ulps within which two values are treated as equal.
    *
    * [[epsilon]] is an *absolute* tolerance, and is only meaningful for values close to 1.0. Callers that compare
    * accumulated log-probabilities work at magnitudes that grow with the size of the input: a value near -500 is
    * routine, and there a single ulp is already ~5.7e-14, some 250x larger than [[epsilon]]. An absolute tolerance
    * therefore silently stops matching anything as the magnitude grows. Comparing in ulps scales with the values
    * being compared, so the same tolerance behaves consistently at any magnitude.
    */
  val maxUlps: Int = 4

  /** The number of representable doubles between `a` and `b`, or `Long.MaxValue` when that is not meaningful.
    *
    * Values of opposite sign, and non-finite values that are not equal, are reported as maximally distant; the
    * absolute [[epsilon]] comparison is what covers values straddling zero.
    */
  private def ulpsBetween(a: Double, b: Double): Long = {
    if (a == b) 0L
    else if (java.lang.Double.isNaN(a) || java.lang.Double.isNaN(b)) Long.MaxValue
    else if (a.isInfinite || b.isInfinite) Long.MaxValue
    else if ((a < 0) != (b < 0)) Long.MaxValue
    else Math.abs(java.lang.Double.doubleToLongBits(a) - java.lang.Double.doubleToLongBits(b))
  }

  /** True if `a` and `b` are within either an absolute `epsilon` or `maxUlps` ulps of one another. */
  private def approximatelyEqual(a: Double, b: Double, epsilon: Double, maxUlps: Int): Boolean = {
    Math.abs(a - b) <= epsilon || ulpsBetween(a, b) <= maxUlps
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
    *                             value then -1 will be returned for the index.
    * @param epsilon The absolute epsilon for comparison of equality; values within [[maxUlps]] ulps of one
    *                another are also treated as equal
    * @param maxUlps The number of ulps within which two values are treated as equal
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN and possibly NegativeInfinity)
    */
  def minWithIndex(ds: Array[Double], allowNegativeInfinity: Boolean=false, requireUniqueMinimum: Boolean=false, epsilon: Double = MathUtil.epsilon, maxUlps: Int = MathUtil.maxUlps): (Double, Int) = {
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the min of a zero length array.")
    def eligible(v: Double): Boolean = {
      !java.lang.Double.isNaN(v) && (Double.NegativeInfinity != v || allowNegativeInfinity)
    }

    var min      = Double.MaxValue
    var minIndex = -1
    var assigned = false
    val len      = ds.length
    var idx      = 0
    while (idx < len) {
      val v = ds(idx)
      if (eligible(v) && (!assigned || v < min)) {
        min = v
        minIndex = idx
        assigned = true
      }

      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("All values are NaNs or negative infinity.")

    // Ties are counted against the final minimum, not a running one: a near-tie is a property of the values, and
    // comparing against a running minimum would only detect it when the tied value happens to appear later.
    if (requireUniqueMinimum && countApproximatelyEqual(ds, min, epsilon, maxUlps, eligible) > 1) (min, -1)
    else (min, minIndex)
  }

  /** Counts the eligible entries of `ds` that are approximately equal to `target`. */
  private def countApproximatelyEqual(ds: Array[Double],
                                      target: Double,
                                      epsilon: Double,
                                      maxUlps: Int,
                                      eligible: Double => Boolean): Int = {
    var count = 0
    var idx   = 0
    while (idx < ds.length) {
      val v = ds(idx)
      if (eligible(v) && approximatelyEqual(v, target, epsilon, maxUlps)) count += 1
      idx += 1
    }
    count
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
    *                             value then -1 will be returned for the index.
    * @param epsilon The absolute epsilon for comparison of equality; values within [[maxUlps]] ulps of one
    *                another are also treated as equal
    * @param maxUlps The number of ulps within which two values are treated as equal
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN)
    */
  def maxWithIndex(ds: Array[Double], requireUniqueMaximum: Boolean=false, epsilon: Double = MathUtil.epsilon, maxUlps: Int = MathUtil.maxUlps): (Double,Int) = {
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the max of a zero length array.")
    def eligible(v: Double): Boolean = !java.lang.Double.isNaN(v)

    var max      =  Double.MinValue
    var maxIndex = -1
    var assigned = false
    val len = ds.length
    var idx = 0
    while (idx < len) {
      val v = ds(idx)
      if (eligible(v) && (!assigned || v > max)) {
        max = v
        maxIndex = idx
        assigned = true
      }
      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("Array contained only NaNs.")

    // Ties are counted against the final maximum, not a running one: a near-tie is a property of the values, and
    // comparing against a running maximum would only detect it when the tied value happens to appear later.
    if (requireUniqueMaximum && countApproximatelyEqual(ds, max, epsilon, maxUlps, eligible) > 1) (max, -1)
    else (max, maxIndex)
  }
}
