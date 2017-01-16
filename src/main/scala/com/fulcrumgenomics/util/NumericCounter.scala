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
 *
 */

package com.fulcrumgenomics.util

import dagr.commons.CommonsDef._
import scala.collection.mutable

object NumericCounter {
  /** Generates a counter that has counted all the numerics provided. */
  def apply[T](ts: TraversableOnce[T])(implicit numeric: Numeric[T]) : NumericCounter[T] = {
    val counter = new NumericCounter[T]
    ts.foreach(counter.count)
    counter
  }

  /** Generates a counter that has counted all the numerics provided. */
  def from[T](ts: TraversableOnce[(T, Long)])(implicit numeric: Numeric[T]) : NumericCounter[T] = {
    val counter = new NumericCounter[T]
    ts.foreach { case (value, count) => counter.count(value, count) }
    counter
  }
}

/**
  * Super-simple class for counting occurrences of any [Numeric]. Will return
  * zero for any item that has not been counted yet.  Implements some useful methods
  * to compute statistics when the objects being counted are numeric types.
  */
class NumericCounter[T](implicit numeric: Numeric[T]) extends SimpleCounter[T] {
  import scala.collection.JavaConversions.mapAsScalaMap
  import numeric._

  /** Create counts such that they are ordered. */
  override protected def makeMap(): mutable.Map[T, Long] = new java.util.TreeMap[T,Long]()

  /** The sum of the products of each numeric and count tuple. */
  private var _totalMass: T = fromInt(0)

  /** Increment the count of object T by the specified value. */
  override def count(t: T, n: Long): Unit = {
    _totalMass += (t * fromInt(n.toInt))
    super.count(t, n)
  }

  /** Returns the mean as a [Double], or zero if there are no counts. */
  def mean(): Double = {
    if (0 == size) return 0.0
    this._totalMass.toDouble() / total
  }

  /** Returns the standard deviation as a [Double], or zero if there are no counts.  */
  def stddev(m: Double = mean()): Double = {
    if (0 == size) return 0.0
    val sum = iterator.map { case (k, v) => v * Math.pow(k.toDouble() - m, 2) }.sum
    if (0 == total) 0
    else if (1== total) Math.sqrt(sum / this._totalMass.toDouble())
    else Math.sqrt(sum / (total - 1.0))
  }

  /** Returns the mean and standar deviation as a [(Double, Double)], or zero for each if there are no counts. */
  def meanAndStddev(): (Double, Double) = {
    val m = mean()
    (m, stddev(m=m))
  }

  /** Returns the key with the greatest count (mode of the distribution), or zero there are no counts.  */
  def mode(): T = {
    if (0 == size) return fromInt(0)
    iterator.maxBy(_._2)._1
  }

  /** Returns the median key as a [Double], or zero there are no counts.  */
  def median(): Double = {
    val count: Long = total

    if (count == 0) 0.0
    else if (count == 1) iterator.next()._1.toDouble()
    else {
      // Find the break point for finding the median value.  If we have an even # of items, we will need to average
      // two values.
      val (midLow: Double, midHigh: Double) = {
        if (count % 2 == 0) {
          val x = count / 2.0
          (x, x + 1)
        }
        else {
          val x = Math.ceil(count / 2.0)
          (x, x)
        }
      }

      // Get the key(s) with the median value.  This is complicated by the fact that we have (key -> counts).
      var midLowValue: Option[T] = None
      var midHighValue: Option[T] = None
      var total: Long = 0
      iterator.exists { case (k, v) =>
        total += v
        if (midLowValue.isEmpty && total >= midLow) midLowValue = Some(k)
        if (midHighValue.isEmpty && total >= midHigh) midHighValue = Some(k)
        midLowValue.nonEmpty && midHighValue.nonEmpty
      }

      (midLowValue, midHighValue) match {
        case (Some(low), Some(high)) => (low + high).toDouble() / 2.0
        case _ => unreachable()
      }
    }
  }

  /** Gets the median absolute deviation. */
  def mad(m: Double = median()): Double = {
    val deviations = new NumericCounter[Double]()
    this.iterator.foreach { case (value, count) =>
      val deviation = Math.abs(value.toDouble() - m)
      deviations.count(deviation, count)
    }
    deviations.median()
  }

  /** Gets sum of values stored in this counter (sum over all (count * value)). */
  def totalMass: T = _totalMass
}

