/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.async.AsyncIterator

import java.util.concurrent.ForkJoinPool
import scala.collection.AbstractIterator
import scala.collection.parallel.{ForkJoinTaskSupport, TaskSupport}

/**
  * Companion object to help with manufacturing ParIterators.
  */
object ParIterator {
  /**
    * Constructs a ParIterator from the underlying iterator, that will perform operations in parallel on chunks
    * of `chunkSize` using `threads` threads.
    */
  def apply[A](iter: Iterator[A], chunkSize: Int = 1024, threads: Int): ParIterator[A] = {
    val pool = new ForkJoinPool(threads, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
    val support = new ForkJoinTaskSupport(pool)
    apply(iter, chunkSize, support)
  }

  /**
    * Constructs a ParIterator from the underlying iterator, that will perform operations in parallel on chunks
    * of `chunkSize` using the provided `TaskSupport`.
    */
  def apply[A](iter: Iterator[A], chunkSize: Int, support: TaskSupport): ParIterator[A] = {
    new ParIterator[A](iter.grouped(chunkSize), chunkSize, support)
  }
}


/**
  * A parallel iterator that operates by grouping or chunking the underlying iterator and then performing
  * operations on the chunks in parallel.  Optionally supports wrapping the iterator in an [[AsyncIterator]]
  * to draw items through the parallel processing and make them available for consumption while the next chunk(s)
  * are processed.
  *
  * For example:
  *
  * ```
  * val iter = Io.readLines(reallyBigFile)
  * ParIterator(iter, chunkSize=4096, threads=8)
  *   .filter(_.contains("foo"))
  *   .map(_.toUpperCase)
  *   .map(_.reverse)
  *   .toAsync(cache=4096 * 8)
  *   .foreach(println)
  *
  * ```
  */
class ParIterator[A] private (private val iter: Iterator[Seq[A]],
                              val chunkSize: Int = 1024,
                              private val taskSupport: TaskSupport) extends AbstractIterator[A] {
  private var currentChunk: Seq[A] = if (iter.hasNext) iter.next() else Seq.empty

  /** True if there are more items available, false otherwise. */
  override def hasNext: Boolean = currentChunk.nonEmpty || iter.hasNext

  /** Retrieve the next item.  Will throw an exception if there are no more items. */
  override def next(): A = {
    while (currentChunk.isEmpty && iter.hasNext) currentChunk = iter.next()
    val v = currentChunk.head
    currentChunk = currentChunk.tail
    v
  }

  /** Returns an iterator over all remaining chunks including the current chunk. */
  private def chunks: Iterator[Seq[A]] = Iterator(currentChunk) ++ iter

  /** Wraps an async iterator around this parallel iterator to draw items through the iterator. This should only
    * be invoked after and other transforming operations such as map/filter/etc.
    *
    * @param cache how many elements to accumulate in the async iterator's cache - for best memory usage this
    *              should be a multiple of chunkSize.
    */
  def toAsync(cache: Int = 4 * chunkSize): AsyncIterator[A] = AsyncIterator(this, Some(cache))

  //////////////////////////////////////////////////////////////////////////////
  // Everything below here is just overrides of transform methods on Iterator
  //////////////////////////////////////////////////////////////////////////////

  override def map[B](f: A => B): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).map(f).seq), chunkSize, taskSupport)

  override def flatMap[B](f: A => IterableOnce[B]): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).flatMap(f).seq), chunkSize, taskSupport)

  override def filter(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).filter(p).seq), chunkSize, taskSupport)

  override def filterNot(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).filterNot(p).seq), chunkSize, taskSupport)

  override def withFilter(p: A => Boolean): ParIterator[A] =
    new ParIterator[A](chunks.map(c => c.parWith(taskSupport).withFilter(p).seq), chunkSize, taskSupport)

  override def collect[B](pf: PartialFunction[A, B]): ParIterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).collect(pf).seq), chunkSize, taskSupport)

  override def flatten[B](implicit ev: A => IterableOnce[B]): Iterator[B] =
    new ParIterator[B](chunks.map(c => c.parWith(taskSupport).flatten.seq), chunkSize, taskSupport)
}
