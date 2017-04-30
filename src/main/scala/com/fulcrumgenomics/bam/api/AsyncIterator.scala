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

package com.fulcrumgenomics.bam.api

import java.util.concurrent.TimeUnit

import scala.collection.mutable
import scala.concurrent.{Await, Future}
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.Duration

class AsyncIterator[A](private val iter: Iterator[A], bufferSize: Int = 1000) extends BufferedIterator[A] {
  private var buffer: mutable.Buffer[A] = mutable.Buffer.empty
  private var bufferPosition: Int = 0
  private var possibleFuture = Option(toTheFuture)
  private val timeout = Duration.create(60, TimeUnit.SECONDS)
  
  revolve()
  
  /** Generates a new future to read the next batch from the underlying iterator. */
  def toTheFuture: Future[mutable.Buffer[A]] = Future { iter.take(bufferSize).toBuffer }
  
  /** If we've hit the end of the buffer, replace it with the buffer from the current future. */
  private def revolve(): Unit = {
    if (bufferPosition == buffer.length) {
      for (future <- possibleFuture) {
        val xs              = Await.result(future, timeout)
        this.buffer         = xs
        this.bufferPosition = 0
      }
      
      possibleFuture = if (iter.hasNext) Some(toTheFuture) else None 
    }
  }
  
  /** Retrieves the current element from the iterator. */
  override def head: A = {
    if (!hasNext) throw new NoSuchElementException("next() called on empty iterator")
    buffer(bufferPosition)
  }

  override def hasNext: Boolean = bufferPosition < buffer.length || possibleFuture.isDefined
  
  override def next(): A = {
    val retval = head
    this.buffer(bufferPosition) = null.asInstanceOf[A]
    this.bufferPosition += 1
    revolve()
    retval
  }
}
