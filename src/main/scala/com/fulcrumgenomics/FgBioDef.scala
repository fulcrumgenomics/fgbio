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

import scala.collection.Parallelizable
import scala.collection.parallel.{ForkJoinTaskSupport, ParIterableLike, TaskSupport}
import scala.concurrent.forkjoin.ForkJoinPool

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
    def bufferBetter = new BetterBufferedIterator(new JavaIteratorAdapter(iterator))
  }

  /** A trait that combine's scala's Iterator with Java's Iterator. */
  trait DualIterator[A] extends java.util.Iterator[A] with Iterator[A]

  /** Class that wraps a Java iterator into a Scala iterator. */
  private final class JavaIteratorAdapter[A](private[FgBioDef] val underlying: java.util.Iterator[A]) extends DualIterator[A] {
    override def hasNext: Boolean = underlying.hasNext
    override def next(): A = underlying.next()
  }

  /** Class that wraps a Scala iterator into a Java iterator. */
  private final class ScalaIteratorAdapter[A](private[FgBioDef] val underlying: Iterator[A]) extends DualIterator[A] {
    override def hasNext: Boolean = underlying.hasNext
    override def next(): A = underlying.next()
  }

  /** Implicit that wraps a Java iterator as a Scala iterator. */
  implicit def javaIteratorAsScalaIterator[A](iterator: java.util.Iterator[A]): DualIterator[A] = iterator match {
    case iter: DualIterator[A]       => iter
    case iter: java.util.Iterator[A] => new JavaIteratorAdapter[A](iter)
  }

  /** Implicit that wraps a Scala iterator as a Java iterator. */
  implicit def scalaIteratorAsJavaIterator[A](iterator: Iterator[A]): DualIterator[A] = iterator match {
    case iter: DualIterator[A] => iter
    case iter: Iterator[A]     => new ScalaIteratorAdapter[A](iter)
  }

  /** Implicit that converts a Java Iterable into a scala Iterator. */
  implicit def javaIterableToIterator[A](iterable: java.lang.Iterable[A]): Iterator[A] = {
    new JavaIteratorAdapter(iterable.iterator())
  }

  /** Implicit class that provides methods for creating Java collections from an Iterator. */
  implicit class IteratorToJavaCollectionsAdapter[A](private val iterator: Iterator[A]) {
    def toJavaList: java.util.List[A]           = fill(new java.util.ArrayList[A]())
    def toJavaSet : java.util.Set[A]            = fill(new java.util.HashSet[A]())
    def toJavaSortedSet: java.util.SortedSet[A] = fill(new java.util.TreeSet[A]())

    private def fill[C <: java.util.Collection[A]](coll: C): C = {
      iterator.foreach(coll.add)
      coll
    }
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

  /**
    * Implicit that provides additional methods to any collection that is Parallelizable.
    * Introduces [[parWith()]] methods that create parallel versions of the collection
    * with various configuration options.
    *
    * @param parallelizable any parallelizable collection
    * @tparam A the type of the elements in the collection
    * @tparam B the type of the parallel representation of the collection
    */
  implicit class ParSupport[A, B <: ParIterableLike[_, _, _]](private val parallelizable: Parallelizable[A,B]) {
    /** Creates a parallel collection with the provided TaskSupport. */
    def parWith(taskSupport: TaskSupport): B = {
      val par = parallelizable.par
      par.tasksupport = taskSupport
      par
    }

    /** Creates a parallel collection with the provided ForkJoinPool. */
    def parWith(pool: ForkJoinPool): B = parWith(taskSupport=new ForkJoinTaskSupport(pool))

    /** Creates a parallel collection with the desired level of parallelism and FIFO semantics. */
    def parWith(parallelism: Int, fifo: Boolean = true): B = {
      parWith(new ForkJoinPool(parallelism, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, fifo))
    }
  }
}

/** Singleton object that extends from the class. */
object FgBioDef extends FgBioDef
