/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

import java.util.concurrent.locks.ReentrantLock
import java.util.concurrent.{BlockingQueue, LinkedBlockingQueue}

import scala.collection.JavaConversions._
import scala.collection.mutable

/**
  * A blocking queue that has a close method.
  *
  * The following methods are blocking: [[put]], [[putAll]], [[putSome]], [[take]], [[takeToMax]], [[takeAtLeast]], and [[close]].
  *
  * Calling [[close]] method will cause all blocking methods to return after the next time the acquire the lock.
  *
  * @param capacity the capacity of the queue, otherwise no capacity will be enforced.
  */
class CloseableBlockingQueue[E](val capacity: Int = Int.MaxValue) {
  private[util] val lock: ReentrantLock = new ReentrantLock()
  private val queue: BlockingQueue[E] = new LinkedBlockingQueue[E](capacity)
  private var closed: Boolean = false

  /** True if [[close]] has been called, false otherwise. */
  def isClosed: Boolean = closed

  /** Puts an element into the queue, waiting if necessary for space to become available.
    *
    * @param e the element to enqueue.
    * @return true if the element was enqueued, false if the queue was closed while blocking.
    */
  def put(e: E): Boolean = {
    var success = false
    while (!success && !closed) {
      lock.lock()
      if (queue.remainingCapacity() > 0) {
        queue.put(e)
        success = true
      }
      lock.unlock()
    }
    success
  }

  /** Puts all the elements into the queue, waiting if necessary for space to become available.
    *
    * If the queue is closed such that this method returns false, some elements may have been enqueued, and as such,
    * this method is not guaranteed to be atomic.
    *
    * @param e the elements to enqueue.
    * @return true if the elements were enqueued, false if the queue was closed while blocking.
    */
  def putAll(e: Iterable[E]): Boolean = {
    val iter = e.iterator
    while (iter.hasNext && !closed) {
      lock.lock()
      while (iter.hasNext && queue.remainingCapacity() > 0) {
        queue.put(iter.next())
      }
      lock.unlock()
    }
    !iter.hasNext // return false if there were some remaining
  }

  /** Puts at least one element, or the minimum amount if given, into the queue, waiting if necessary for space to
    * become available.
    *
    * If the number of elements enqueued is less than `min`, the queue was closed.
    *
    * @param e the elements to enqueue.
    * @param min the minimum number to enqueue, or 1 if not specified.
    * @return the number of elements enqueued.
    */
  def putSome(e: Iterable[E], min: Int = 1): Int = {
    val iter = e.iterator
    var n = 0
    var m = n + 1
    // if min == 0 and n == 0, try the first time
    while (iter.hasNext && n < min && !closed) {
      lock.lock()
      // TODO: use addAll
      while (iter.hasNext && queue.remainingCapacity() > 0) {
        queue.put(iter.next())
        n += 1
      }
      m = min // use `min` for subsequent iterations of this loop
      lock.unlock()
    }
    n
  }

  /** Takes an element from the queue, waiting if necessary until an element becomes available.
    *
    * @return the element dequeued, or None if the queue was closed while waiting.
    */
  def take(): Option[E] = {
    var ret: Option[E] = None
    while (ret.isEmpty && !closed) {
      lock.lock()
      if (queue.nonEmpty) {
        ret = Some(queue.take())
      }
      lock.unlock()
    }
    ret
  }

  /** Takes one or more elements from the queue, up to the maximum number, waiting if necessary until an element becomes
    * available.
    *
    * The first time at least one element is available, all available elements, up to the maximum number, will be take
    * and stored in `buffer`.
    *
    * @param buffer the buffer to store the dequeued elements.
    * @param max the maximum number to take, or as many as possible if not specified.
    * @return the number of elements dequeued, with zero indicating the queue was closed.
    */
  def takeToMax(buffer: mutable.AbstractBuffer[E], max: Int = Int.MaxValue): Int = {
    var n = 0
    while (n == 0 && !closed) {
      lock.lock()
      if (queue.nonEmpty) {
        n = queue.drainTo(buffer, max)
      }
      lock.unlock()
    }
    n
  }

  /** Takes at least the specified number of elements from the queue, waiting if necessary until an element becomes
    * available.
    *
    * @param buffer the buffer to store the dequeued elements.
    * @param min the minimum number to take, or 1 if not specified.
    * @return the number of elements dequeued, with zero indicating the queue was closed.
    */
  def takeAtLeast(buffer: mutable.AbstractBuffer[E], min: Int = 1): Int = {
    var n = 0
    while (n < min && !closed) {
      lock.lock()
      if (queue.nonEmpty) {
        n += queue.drainTo(buffer)
      }
      lock.unlock()
    }
    n
  }

  /** Closes this queue, such that any subsequent call to enqueue or dequeue elements will not change this queue.
    * Furthermore, any blocked calls will return after obtaining the lock.
    */
  def close(): Unit = {
    lock.lock()
    closed = true
    lock.unlock()
  }

  def isEmpty: Boolean = queue.isEmpty

  def nonEmpty: Boolean = queue.nonEmpty

  /** The number elements currently in this queue. */
  def size(): Int = queue.size()

  /** The remaining capacity of this queue (for more elements). */
  def remainingCapacity: Int = queue.remainingCapacity()

  /** Useful for testing. */
  private[util] def toSeq: Seq[E] = {
    lock.lock()
    val s = queue.toSeq
    if (queue.isEmpty) queue.addAll(s)
    lock.unlock()
    s
  }
}
