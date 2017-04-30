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
package com.fulcrumgenomics.fastq

import java.io.{BufferedWriter, Closeable}
import java.nio.file.Path

import com.fulcrumgenomics.Writer
import com.fulcrumgenomics.util.Io

/**
  * Companion object for manufacturing FastqWriter instances.
  */
object FastqWriter {
  /** Constructs a FastqWriter that will write to the provided path. */
  def apply(path: Path): FastqWriter = apply(Io.toWriter(path))

  /** Constructs a FastqWriter from a Writer. */
  def apply(writer: java.io.Writer): FastqWriter = {
    if (writer.isInstanceOf[BufferedWriter]) new FastqWriter(writer.asInstanceOf[BufferedWriter])
    else new FastqWriter(new BufferedWriter(writer))
  }
}

/**
  * Implements a writer for fastq records.
  */
class FastqWriter private(val out: BufferedWriter) extends Closeable with Writer[FastqRecord] {
  /** Writes a single record to the output. */
  override def write(rec: FastqRecord): Unit = {
    out.append('@').append(rec.header).append('\n')
    out.append(rec.bases).append('\n')
    out.append('+').append('\n')
    out.append(rec.quals).append('\n')
  }

  /** Closes the underlying writer. */
  override def close(): Unit = out.close()
}
