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

import scala.collection.mutable

object SimpleCounter {
  /** Generates a counter that has counted all the objects provided. */
  def apply[T](ts: TraversableOnce[T]): SimpleCounter[T] = {
    val counter = new SimpleCounter[T]
    ts.foreach(counter.count)
    counter
  }
}

/**
  * Super-simple class for counting occurrences of any kind of object. Will return
  * zero for any item that has not been counted yet.  
  */
class SimpleCounter[T] extends Iterable[(T, Long)] {
  private val counts = mutable.Map[T,Long]()

  /** Increment the count of object T by 1. */
  def count(t: T): Unit = count(t, 1)

  /** Increment the count of object T by the specified value. */
  def count(t: T, n: Long): Unit = counts(t) = counts.getOrElse(t, 0L) + n

  /** Gets the count of the provided object. */
  def countOf(t: T): Long = this.counts.getOrElse(t, 0)

  /** Iterates over all the objects and their counts. */
  override def iterator: Iterator[(T, Long)] = this.counts.iterator
}
