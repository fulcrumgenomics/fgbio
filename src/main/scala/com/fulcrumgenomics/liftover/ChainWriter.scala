/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

package com.fulcrumgenomics.liftover

import java.io.{BufferedWriter, Closeable}
import java.nio.file.Path

import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.util.Io

/**
  * Companion object for manufacturing [[ChainWriter]] instances.
  */
object ChainWriter {
  /** Constructs a [[ChainWriter]] that will write to the provided path. */
  def apply(path: Path): ChainWriter = apply(Io.toWriter(path))

  /** Constructs a [[ChainWriter]] from a Writer. */
  def apply(writer: java.io.Writer): ChainWriter = writer match {
    case bufferedWriter: BufferedWriter => new ChainWriter(bufferedWriter)
    case otherWriter                    => new ChainWriter(new BufferedWriter(otherWriter))
  }
}

/**
  * Implements a writer for [[Chain]] records.
  */
class ChainWriter private(val out: BufferedWriter)  extends Closeable with Writer[Chain] {

  /** Writes a single record to the output. */
  override def write(chain: Chain): Unit = {
    out.write(chain.header.toLine)
    out.write('\n')
    chain.blocks.foreach { block =>
      out.write(block.toLine)
      out.write('\n')
    }
    out.write('\n')
  }

  /** Closes the underlying writer. */
  override def close(): Unit = out.close()
}
