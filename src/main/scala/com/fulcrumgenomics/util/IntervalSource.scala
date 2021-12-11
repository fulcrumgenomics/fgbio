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
import com.fulcrumgenomics.util.IntervalListSource.{HeaderPrefix => IntervalListHeaderPrefix}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.Interval
import htsjdk.tribble.annotation.Strand

import java.io.{Closeable, File, InputStream}
import scala.io.Source

/** A class for sourcing intervals from a stream of data that could either be in BED or Interval List format. */
class IntervalSource private(
  private val lines: Iterator[String],
  private val sd: Option[SequenceDictionary],
  private val source: Option[{ def close(): Unit }] = None
) extends Iterator[Interval] with Closeable {

  /** The underlying buffered iterator of interval data. */
  private val iter: BetterBufferedIterator[String] = lines match {
    case iter: BetterBufferedIterator[String] => iter
    case iter                                 => iter.bufferBetter
  }

  private val (underlying: Iterator[Interval], _dict, _header) = if (
    iter.headOption.exists(_.startsWith(IntervalListHeaderPrefix))
  ) {
    val wrapped = IntervalListSource(iter)
    require(sd.forall(wrapped.dict.sameAs), "Provided sequence dictionary does not match the input's dict header!")
    (wrapped, Some(wrapped.dict), Some(wrapped.header))
  } else {
    val wrapped = BedSource(iter, sd).map { bed =>
      new Interval(
        bed.getContig,
        bed.getStart,
        bed.getEnd,
        // BEDFeature.getStrand() can be null so wrap in an option and search for the negative enum.
        Option(bed.getStrand).contains(Strand.NEGATIVE),
        // The default name for BEDFeature is the empty string (""), but defaults to null for Interval.
        if (bed.getName == "") null else bed.getName
      )
    }
    (wrapped, sd, None)
  }

  /** The [[SAMFileHeader]] associated with the source, if it exists. */
  val header: Option[SAMFileHeader] = _header

  /** The [[SequenceDictionary]] associated with the source, if it exists. */
  val dict: Option[SequenceDictionary] = _dict

  /** True if calling `next()` will yield another interval, false otherwise. */
  override def hasNext: Boolean = underlying.hasNext

  /** Returns the next interval if available, or throws an exception if none is available. */
  override def next(): Interval = underlying.next()

  /** Closes the underlying reader. */
  override def close(): Unit = this.source.foreach(_.close())
}

/** Companion object for [[IntervalSource]]. */
object IntervalSource {

  /** Creates a new interval source from a sequence of lines. */
  def apply(lines: Iterable[String], dict: Option[SequenceDictionary]): IntervalSource = new IntervalSource(lines.iterator, dict)

  /** Creates a new interval source from an iterator of lines. */
  def apply(lines: Iterator[String], dict: Option[SequenceDictionary]): IntervalSource = new IntervalSource(lines, dict)

  /** Creates a new interval source from an input stream. */
  def apply(stream: InputStream, dict: Option[SequenceDictionary]): IntervalSource = {
    new IntervalSource(Source.fromInputStream(stream).getLines(), dict)
  }

  /** Creates a new interval source from a source. */
  def apply(source: Source, dict: Option[SequenceDictionary]): IntervalSource = {
    new IntervalSource(source.getLines(), dict, source = Some(source))
  }

  /** Creates a new interval source from a File. */
  def apply(file: File, dict: Option[SequenceDictionary]): IntervalSource = apply(path = file.toPath, dict)

  /** Creates a new interval source from a Path. */
  def apply(path: PathToIntervals, dict: Option[SequenceDictionary]): IntervalSource = apply(Io.readLines(path), dict)
}
