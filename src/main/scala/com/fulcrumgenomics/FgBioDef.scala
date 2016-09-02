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

package com.fulcrumgenomics

import com.fulcrumgenomics.util.BetterBufferedIterator
import dagr.commons.CommonsDef

/**
  * Place to put common function, type and implicit definitions that can be
  * imported into other classes easily.
  */
class FgBioDef extends CommonsDef {
  /** Implicit class that provides a method to wrap an iterator into a BetterBufferedIterator. */
  implicit class BetterBufferedIteratorScalaWrapper[A](val iterator: Iterator[A]) {
    def bufferBetter = new BetterBufferedIterator(iterator)
  }

  /** Implicit class that provides a method to wrap a Java iterator into a BetterBufferedIterator. */
  implicit class BetterBufferedIteratorJavaWrapper[A](val iterator: java.util.Iterator[A]) {
    def bufferBetter = new BetterBufferedIterator(scala.collection.JavaConversions.asScalaIterator(iterator))
  }

  /**
    * An implementation of a for loop that has performance similar to writing a custom
    * while loop for primitive types, which suffer great performance loss when iterating
    * via foreach(), and especially with zipWithIndex.foreach().
    *
    * Equivalent to: for(int i=from; i<until; i+=by) f(i)
    *
    * @param from an initial integer values
    * @param until a value one higher than the last value that should be accepted
    * @param by a step value to increment by each iteration
    * @param f a function that takes the index and is called each iteration
    */
  @inline def forloop(from: Int, until: Int, by: Int=1)(f: Int => Unit): Unit = {
    val comp: (Int,Int) => Boolean = if (by > 0) (_ < _) else (_ > _)
    var i = from
    while (comp(i, until)) {
      f(i)
      i += by
    }
  }

  /**
    * A more generic for loop, that gives performance similar, but slightly worse than,
    * FgBioDef#forloop(Int,Int,Int).  Equivalent to:
    *
    *  for (i=from; check(i); i=next(i)) f(i)
    *
    * @param from an initial value
    * @param check a function that checks to see if a value should be evaluated (not cause termination)
    * @param next a function that takes a value and produced the next value
    * @param f a function called on all values
    * @tparam T the type of the index
    */
  @inline def forloop[@specialized T](from: T)(check: T => Boolean)(next: T => T)(f: T => Unit): Unit = {
    var t = from
    while (check(t)) {
      f(t)
      t = next(t)
    }
  }
}

/** Singleton object that extends from the class. */
object FgBioDef extends FgBioDef
