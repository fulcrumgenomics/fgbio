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

import java.io.{Closeable, InputStream}
import java.util.concurrent.atomic.AtomicInteger
import dagr.commons.CommonsDef._

import scala.io.Source

object AsyncStreamSink {
  private val n: AtomicInteger = new AtomicInteger(1)
  private def nextName: String = "AsyncStreamSinkThread-" + n.getAndIncrement
}

/**
  * Class that is capable of reading output from a sub-process (or arbitrary other InputStream)
  * and processing it through any method that accepts a String, while managing the extra thread
  * necessary to ensure nothing blocks.
  */
class AsyncStreamSink(in: InputStream, private val sink: String => Unit) extends Closeable {
  private val source = Source.fromInputStream(in).withClose(() => in.close())
  private val thread = new Thread( new Runnable { def run() : Unit = { source.getLines().foreach(sink) } } )
  this.thread.setName(AsyncStreamSink.nextName)
  this.thread.setDaemon(true)
  this.thread.start()

  /** Sets the name on the underlying thread so that it's more easily identifiable. */
  def setName(name: String): AsyncStreamSink = yieldAndThen(this) { this.thread.setName(name) }

  /** Give the thread 500 seconds to wrap up what it's doing and then interrupt it. */
  def close() : Unit = close(500)

  /** Give the thread `millis` milliseconds to finish what it's doing, then interrupt it. */
  def close(millis: Long) : Unit = {
    this.thread.join(millis) // Give it some time to wrap up what it's doing
    this.thread.interrupt()  // then interrupt it
  }
}
