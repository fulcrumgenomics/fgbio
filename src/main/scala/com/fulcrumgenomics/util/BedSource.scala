/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics
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

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.tribble.bed.BEDCodec.StartOffset
import htsjdk.tribble.bed.{BEDCodec, BEDFeature}

import java.io.{Closeable, File, InputStream}
import scala.io.Source
import scala.util.{Failure, Success, Try}

/** A class for sourcing BED features from a stream of ASCII string data. */
class BedSource private(
  private val lines: Iterator[String],
  private val sd: Option[SequenceDictionary] = None,
  private val source: Option[{ def close(): Unit}] = None
) extends Iterator[BEDFeature] with Closeable {

  /** The underlying codec used to parse the lines of BED data. */
  private val codec = new BEDCodec(StartOffset.ONE)

  /** The current line count. */
  private var lineNumber = 1L

  /** The underlying buffered iterator of BED data. */
  private val iter: BetterBufferedIterator[String] = lines match {
    case iter: BetterBufferedIterator[String] => iter
    case iter                                 => iter.bufferBetter
  }

  /** The header of this BED file. In most cases, the BED header is empty. */
  val header: Seq[String] = {
    val lines = iter.takeWhile(line => BedSource.HeaderPrefixes.exists(line.startsWith)).toIndexedSeq
    lineNumber += lines.length
    lines
  }

  /** The [[SequenceDictionary]] associated with the source. */
  val dict: Option[SequenceDictionary] = sd

  /** True if calling `next()` will yield another BED feature, false otherwise. */
  override def hasNext: Boolean = iter.hasNext

  /** Returns the next BED feature if available, or throws an exception if the feature is invalid or none is available. */
  override def next(): BEDFeature = yieldAndThen(parse(iter.next()))(lineNumber += 1)

  /** Parse a line of text and build a BED feature. */
  private def parse(line: String): BEDFeature = {
    val parsed  = codec.decode(line)
    val feature = Option(parsed).getOrElse(throw new IllegalStateException(s"No BED feature could be built from line number: $lineNumber"))
    Try(dict.foreach(_.validate(feature))) match {
      case Success(())                          => feature
      case Failure(e: NoSuchElementException)   => throw new NoSuchElementException(e.getMessage + s" Failed on line number: $lineNumber")
      case Failure(e: IllegalArgumentException) => throw new IllegalArgumentException(e.getMessage + s" Failed on line number: $lineNumber")
      case Failure(e: Throwable)                => throw new IllegalStateException(e.getMessage + s" Failed on line number: $lineNumber")
    }
  }

  /** Closes the optional underlying handle. */
  override def close(): Unit = this.source.foreach(_.close())
}

/** Companion object for [[BedSource]]. */
object BedSource {

  /** Common BED header line prefixes. */
  val HeaderPrefixes: Seq[String] = Seq("#", "browser", "track")

  /** Creates a new BED source from a sequence of lines. */
  def apply(lines: Iterable[String], dict: Option[SequenceDictionary]): BedSource = new BedSource(lines.iterator, dict)

  /** Creates a new BED source from an iterator of lines. */
  def apply(lines: Iterator[String], dict: Option[SequenceDictionary]): BedSource = new BedSource(lines, dict)

  /** Creates a new BED source from an input stream. */
  def apply(stream: InputStream, dict: Option[SequenceDictionary]): BedSource = {
    new BedSource(Source.fromInputStream(stream).getLines(), dict)
  }

  /** Creates a new BED source from a source. */
  def apply(source: Source, dict: Option[SequenceDictionary]): BedSource = {
    new BedSource(source.getLines(), dict, source = Some(source))
  }

  /** Creates a new BED source from a File. */
  def apply(file: File, dict: Option[SequenceDictionary]): BedSource = apply(path=file.toPath, dict)

  /** Creates a new BED source from a Path. */
  def apply(path: PathToIntervals, dict: Option[SequenceDictionary]): BedSource = apply(Io.readLines(path), dict)
}
