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
package com.fulcrumgenomics.fastq

import java.io.Closeable

import com.fulcrumgenomics.util.Io
import dagr.commons.CommonsDef._

object SequentialFastqSource {
  /** Creates a new fastq source from a Path. */
  def apply(paths: Seq[PathToFastq]): SequentialFastqSource = {
    new SequentialFastqSource(paths.map(path => FastqSource(Io.toInputStream(path))).iterator)
  }
}

/** Reads multiple FASTQ sources one after the other. This is useful for when we have more than one FASTQ for a
  * set of reads. */
class SequentialFastqSource(sources: Iterator[FastqSource]) extends Iterator[FastqRecord] with Closeable {
  val _sources = sources.toSeq
  val iters: Iterator[FastqRecord] = _sources.iterator.flatten
  def hasNext(): Boolean = iters.hasNext
  def next(): FastqRecord = iters.next()
  def close(): Unit = _sources.foreach(_.close())
}
