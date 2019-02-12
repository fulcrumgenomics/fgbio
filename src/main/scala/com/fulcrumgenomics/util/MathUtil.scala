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

import java.util.NoSuchElementException

/**
  * Some simple utility methods for various mathematical operations that are implemented
  * in efficient albeit non-idiomatic scala.
  */
object MathUtil {
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
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN and possibly NegativeInfinity)
    */
  def minWithIndex(ds: Array[Double], allowNegativeInfinity: Boolean=false, requireUniqueMinimum: Boolean=false): (Double,Int) = {
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the min of a zero length array.")
    var min      = Double.MaxValue
    var minIndex = -1
    var assigned = false
    val len      = ds.length
    var idx      = 0
    while (idx < len) {
      val v = ds(idx)
      if (!java.lang.Double.isNaN(v) && (Double.NegativeInfinity != v || allowNegativeInfinity)) {
        if (!assigned || v < min) {
          min = v
          minIndex = idx
          assigned = true
        }
        else if (v == min && requireUniqueMinimum) {
          minIndex = -1
        }
      }

      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("All values are NaNs or negative infinity.")
    (min, minIndex)
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
    * @throws java.util.NoSuchElementException if either the input array is zero length or the array
    *                                contains only invalid values (NaN)
    */
  def maxWithIndex(ds: Array[Double], requireUniqueMaximum: Boolean=false): (Double,Int) = {
    if (ds.length == 0) throw new NoSuchElementException("Cannot find the max of a zero length array.")
    var max      =  Double.MinValue
    var maxIndex = -1
    var assigned = false
    val len = ds.length
    var idx = 0
    while (idx < len) {
      val v = ds(idx)
      if (!java.lang.Double.isNaN(v)) {
        if (!assigned || v > max) {
          max = v
          maxIndex = idx
          assigned = true
        }
        else if (v == max && requireUniqueMaximum) {
          maxIndex = -1
        }
      }
      idx += 1
    }

    if (!assigned) throw new NoSuchElementException("Array contained only NaNs.")
    (max, maxIndex)
  }
}
