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

import java.io._

import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.commons.CommonsDef.{PathToFastq, yieldAndThen}

import scala.io.Source

/**
  * Provides factory methods for creating a `FastqSource` from multiple types.
  */
object FastqSource {
  /** Creates a new fastq source from a sequence of lines. */
  def apply(lines: Seq[String]): FastqSource = new FastqSource(lines.iterator)

  /** Creates a new fastq source from an iterator of lines. */
  def apply(lines: Iterator[String]): FastqSource = new FastqSource(lines)

  /** Creates a new fastq source from an input stream. */
  def apply(stream: InputStream): FastqSource = new FastqSource(Source.fromInputStream(stream).getLines(), Some(stream))

  /** Creates a new fastq source from a source. */
  def apply(source: Source): FastqSource = new FastqSource(source.getLines(), Some(source))

  /** Creates a new fastq source from a File. */
  def apply(file: File): FastqSource = apply(path=file.toPath)

  /** Creates a new fastq source from a Path. */
  def apply(path: PathToFastq): FastqSource = apply(Io.toInputStream(path))

  /** Returns an iterator over multiple fastq files that ensures:
    *   1. Either all sources or no sources have more records
    *   2. That the next records from each of the sources have the same read name
    *
    * @param sources a Seq of one or more FastqSource objects
    * @return an Iterator that returns a Seq of FastqRecord, one per source
    */
  def zipped(sources: Seq[FastqSource]): Iterator[Seq[FastqRecord]] = new Iterator[Seq[FastqRecord]] {
    require(sources.nonEmpty, "No sources provided")

    def hasNext(): Boolean = sources.exists(_.hasNext)

    def next(): Seq[FastqRecord] = {
      if (!this.hasNext) throw new NoSuchElementException("Calling next() when hasNext() is false.")
      require(sources.forall(_.hasNext) == sources.head.hasNext, "Sources are out of sync.")
      val records = sources.map(_.next)
      // Check that the FASTQ records all have the same name
      require(records.forall(_.name == records.head.name), "Fastqs are out of sync, found read names: " + records.map(_.name).mkString(", "))
      records
    }
  }
}


/**
  * Reads fastq records from any text based source via a reader. Ensures that lines come in
  * the expected groupings of four lines with the correct headers, and that bases and qualities
  * are of the same length.  The underlying reader is closed automatically when EOF is reached.
  */
class FastqSource private(val lines: Iterator[String],
                          private[this] val source: Option[{ def close(): Unit }] = None)
  extends Iterator[FastqRecord] with Closeable {
  
  private var nextRecord: Option[FastqRecord] = fetchNextRecord()

  /** True if calling `next()` will yield another record, false otherwise. */
  override def hasNext: Boolean = this.nextRecord.isDefined

  /** Returns the next record if available, or throws an exception if none is available. */
  override def next(): FastqRecord =  yieldAndThen(nextRecord.get) { this.nextRecord = fetchNextRecord() }

  /** Short hand to throw an illegal state exception. */
  private def illegal(msg: String): Nothing = throw new IllegalStateException(msg)

  /** Sets the current record to None, and then attempts to read the next record from the input. */
  private def fetchNextRecord() : Option[FastqRecord] = {
    this.lines.take(4).toList match {
      case Nil =>
        None
      case header :: seq :: qheader :: quals :: Nil =>
        if (!header.startsWith("@"))    illegal(s"Fastq sequence header must start with @: ${header}")
        if (!qheader.startsWith("+"))   illegal(s"Fastq quality header must start with +: ${qheader}")
        if (seq.length != quals.length) illegal(s"Sequence and qualities not same length for record: ${header}")

        // Destructure the header line into name, read number and comment
        val (fullName, comment) = header.indexWhere(_.isWhitespace) match {
          case -1 => (header.drop(1), None)
          case i  => (header.substring(1, i), Some(header.substring(i+1)))
        }

        val (name, readNumber) = fullName.endsWith("/1") || fullName.endsWith("/2") match {
          case true  => (fullName.dropRight(2), Some(fullName.last.asDigit))
          case false => (fullName, None)
        }

        Some(FastqRecord(
          name       = name,
          bases      = seq,
          quals      = quals,
          comment    = comment,
          readNumber = readNumber
        ))
      case header :: _ =>
        illegal(s"Fastq source terminated mid-record at ${header}")
    }
  }

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  override def close(): Unit = this.source.foreach(_.close())
}
