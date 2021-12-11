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

package com.fulcrumgenomics.util


import java.io.{Closeable, File, InputStream}
import com.fulcrumgenomics.FgBioDef.{PathToIntervals, yieldAndThen}
import com.fulcrumgenomics.commons.CommonsDef.BetterBufferedIteratorScalaWrapper
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.StringUtil
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.util.IntervalListSource.HeaderPrefix
import htsjdk.samtools.util.{BufferedLineReader, Interval, IntervalList}
import htsjdk.samtools.{SAMFileHeader, SAMTextHeaderCodec}

import scala.io.Source
import scala.util.{Failure, Success, Try}

/** Companion object for [[IntervalListSource]]. */
object IntervalListSource {

  /** The Interval List header line prefix. */
  val HeaderPrefix: String = "@"

  /** Creates a new interval list source from a sequence of lines. */
  def apply(lines: Iterable[String]): IntervalListSource = new IntervalListSource(lines.iterator)

  /** Creates a new interval list source from an iterator of lines. */
  def apply(lines: Iterator[String]): IntervalListSource = new IntervalListSource(lines)

  /** Creates a new interval list source from an input stream. */
  def apply(stream: InputStream): IntervalListSource = new IntervalListSource(Source.fromInputStream(stream).getLines())

  /** Creates a new interval list source from a source. */
  def apply(source: Source): IntervalListSource = new IntervalListSource(source.getLines(), Some(source))

  /** Creates a new interval list source from a File. */
  def apply(file: File): IntervalListSource = apply(path=file.toPath)

  /** Creates a new interval list source from a Path. */
  def apply(path: PathToIntervals): IntervalListSource = apply(Io.readLines(path))
}

/**
  * Reads intervals from any text based source via a reader.  The underlying reader is closed automatically when EOF is
  * reached.
  */
class IntervalListSource private(lines: Iterator[String],
                                 private[this] val source: Option[{ def close(): Unit }] = None)
  extends Iterator[Interval] with Closeable {

  private val iter: BetterBufferedIterator[String] = lines match {
    case iter: BetterBufferedIterator[String] => iter
    case iter                                 => iter.bufferBetter
  }

  private var lineNumber = 1L

  /** True if calling `next()` will yield another interval, false otherwise. */
  override def hasNext: Boolean = iter.nonEmpty

  /** Returns the next interval if available, or throws an exception if none is available. */
  override def next(): Interval = yieldAndThen(parse(iter.next())) { lineNumber += 1 }

  // Read the header
  val header: SAMFileHeader = {
    val codec = new SAMTextHeaderCodec
    val headerLines = iter.takeWhile(_.startsWith(HeaderPrefix)).toIndexedSeq
    require(headerLines.nonEmpty, "No header found")
    lineNumber += headerLines.length
    val lineReader = BufferedLineReader.fromString(headerLines.mkString("\n"))
    codec.decode(lineReader, "IntervalListSource")
  }

  require(this.dict.nonEmpty, "No reference sequences found in the header.")

  /** The [[SequenceDictionary]] associated with the source. */
  lazy val dict: SequenceDictionary = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    this.header.getSequenceDictionary.fromSam
  }

  private val parseArray = Array[String]("", "", "", "", "")

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  override def close(): Unit = this.source.foreach(_.close())

  private def parse(line: String): Interval = {
    val fieldCount = StringUtil.split(line, '\t', parseArray)
    require(fieldCount == 5, s"Expected 5 fields on line $lineNumber")
    val Array(refName: String, start: String, end: String, strand: String, name: String) = parseArray

    val negative = strand match {
      case "-" => true
      case "+" => false
      case _   => throw new IllegalArgumentException(s"Unrecognized strand '$strand' on line number $lineNumber")
    }

    val interval = new Interval(refName, start.toInt, end.toInt, negative, name)

    Try(dict.validate(interval)) match {
      case Success(())                          => interval
      case Failure(e: NoSuchElementException)   => throw new NoSuchElementException(e.getMessage + s" Failed on line number: $lineNumber")
      case Failure(e: IllegalArgumentException) => throw new IllegalArgumentException(e.getMessage + s" Failed on line number: $lineNumber")
      case Failure(e: Throwable)                => throw new IllegalStateException(e.getMessage + s" Failed on line number: $lineNumber")
    }
  }

  /** Reads in the intervals into an [[htsjdk.samtools.util.IntervalList]] */
  def toIntervalList: IntervalList = {
    import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
    val list = new IntervalList(dict.asSam)
    this.foreach { list.add }
    list
  }
}