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

import java.io.{Closeable, _}

import com.fulcrumgenomics.commons.CommonsDef.{BetterBufferedIteratorScalaWrapper, FilePath, yieldAndThen}
import com.fulcrumgenomics.util.Io

import scala.io.Source

/**
  * Provides factory methods for creating a [[ChainSource]] from multiple types.
  */
object ChainSource {
  /** Creates a new chain source from a sequence of lines. */
  def apply(lines: Seq[String]): ChainSource = new ChainSource(lines.iterator)

  /** Creates a new chain source from an iterator of lines. */
  def apply(lines: Iterator[String]): ChainSource = new ChainSource(lines)

  /** Creates a new chain source from an input stream. */
  def apply(stream: InputStream): ChainSource = new ChainSource(Source.fromInputStream(stream).getLines(), Some(stream))

  /** Creates a new chain source from a source. */
  def apply(source: Source): ChainSource = new ChainSource(source.getLines(), Some(source))

  /** Creates a new chain source from a File. */
  def apply(file: File): ChainSource = apply(path=file.toPath)

  /** Creates a new chain source from a Path. */
  def apply(path: FilePath): ChainSource = apply(Io.toInputStream(path))
}


class ChainSource(val lines: Iterator[String],
                  private[this] val source: Option[{ def close(): Unit }] = None
                 ) extends Iterator[Chain] with Closeable {

  private val iter = lines.bufferBetter

  private var nextRecord: Option[Chain] = fetchNextRecord()

  /** True if calling `next()` will yield another record, false otherwise. */
  override def hasNext: Boolean = this.nextRecord.isDefined

  /** Returns the next record if available, or throws an exception if none is available. */
  override def next(): Chain =  yieldAndThen(nextRecord.get) { this.nextRecord = fetchNextRecord() }

  /** Sets the current record to None, and then attempts to read the next record from the input. */
  private def fetchNextRecord() : Option[Chain] = {
    val lines = iter.takeWhile(_.nonEmpty).toSeq // get all the lines!
    iter.next() // consume the empty line
    lines match {
      case Nil => None
      case headerLine :: blockLines =>
        val header      = ChainHeader(headerLine)
        var sourceStart = header.source.start
        var targetStart = header.target.start
        val blocks = blockLines.map { line =>
          val block: ChainAlignment = ChainAlignment(line, sourceStart, targetStart)
          sourceStart += block.size + block.sourceStart
          targetStart += block.size + block.targetStart
          block
        }
        Some(Chain(header = header, blocks = blocks.toIndexedSeq))
    }
  }

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  override def close(): Unit = this.source.foreach(_.close())
}
