/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.alignment

import java.io._
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.commons.util.StringUtil

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

/** Base trait for all aligners that return a `Future[Alignment]`. */
trait AsyncAligner {

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String)(implicit e :ExecutionContext): Future[Alignment] = align(query.getBytes, target.getBytes)

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Future]] that when complete returns an [[Alignment]] object describing the optimal global alignment
    *         of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte])(implicit e :ExecutionContext): Future[Alignment]

  /** The number of completed alignments.  Returns -1 if unknown. */
  def numAligned: Long = -1
}

/** A wrapper around [[Aligner]] to return [[Alignment]]s asynchronously. */
class AsyncScalaAligner(val aligner: Aligner) extends AsyncAligner {
  private val aligned = new AtomicLong(0)

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Future]] that when complete returns an [[Alignment]] object describing the optimal global alignment
    *         of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte])(implicit e :ExecutionContext): Future[Alignment] = Future {
    val alignment = aligner.align(query, target)
    aligned.getAndIncrement()
    alignment
  }

  /** The number of completed alignments.  Returns -1 if unknown. */
  override def numAligned: Long = aligned.get()
}

private object AsyncKswAligner {

  // Developer Note: this is faster than jus split.  See [[DelimitedDataParser]].
  /** Parses an alignment result from the ksw output (a single line). */
  def toAlignment(line: String, delimiter: Char = '\t'): Alignment = {
    val tmp = new Array[String](8)
    val count = StringUtil.split(line, delimiter, tmp)
    require(count == 8, s"Could not parse line into eight: $line.")
    Alignment(
      query       = tmp(6).getBytes,
      target      = tmp(7).getBytes,
      queryStart  = tmp(1).toInt + 1, // make it one-based
      targetStart = tmp(3).toInt + 1, // make it one-based
      cigar       = Cigar(tmp(5)),
      score       = tmp(0).toInt
    )
  }
}

/** An asynchronous aligner that wraps the [ksw](https://github.com/nh13/ksw) executable for faster alignments.
  *
  * Install [the ksw executable version 0.2.0](https://github.com/nh13/ksw) manually or via conda
  * (`conda install -c bioconda ksw=0.2.0`)
  *
  * IMPORTANT: The alignments will be completed in the order added.
  *
  * @param executable the path to the `ksw` executable.
  * @param matchScore the match score (should be greater than or equal to zero).
  * @param mismatchScore the mismatch score (should be greater than or equal to zero).
  * @param gapOpen the gap opening penalty, should generally be negative or zero
  * @param gapExtend the gap extension penalty, should generally be negative or zero
  * @param mode alignment mode to use when generating alignments
  * @param buffer the size of the input and output buffer in bytes.
  */
class AsyncKswAligner(executable: FilePath,
                      matchScore: Int,
                      mismatchScore: Int,
                      gapOpen: Int,
                      gapExtend: Int,
                      mode: Mode = Mode.Glocal,
                      buffer: Int = 1024*1024) extends AsyncAligner with Closeable {
  private val aligned = new AtomicLong(0)
  private var closed = false

  /** The ksw process. */
  private val process = {
    val alignmentMode = this.mode match {
      case Mode.Local  => 0
      case Mode.Glocal => 1
      case Mode.Global => 3
    }

    val args = Seq(
      executable,
      "-a", math.abs(matchScore),
      "-b", math.abs(mismatchScore),
      "-q", math.abs(gapOpen),
      "-r", math.abs(gapExtend),
      "-M", alignmentMode,
      "-c", // add the cigar
      "-s"  // add the query and target
    ).map(_.toString)
    new ProcessBuilder(args:_*).redirectErrorStream(true).start()
  }

  /** The input stream into ksw. */
  private val kswInput  = new PrintStream(new BufferedOutputStream(this.process.getOutputStream, buffer), true)

  /** The output stream out of ksw */
  private val kswOutput = new BufferedReader(new InputStreamReader(this.process.getInputStream), buffer)

  private var previousAlignment: Future[Alignment] = Future.successful {
    new Alignment(query = Array.empty, target = Array.empty, 0, 0, Cigar.empty, 0)
  }

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Future]] that when complete returns an [[Alignment]] object describing the optimal global alignment
    *         of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte])(implicit e :ExecutionContext): Future[Alignment] = previousAlignment.synchronized {
    if (closed) return Future.failed(new IllegalStateException("ksw is closed but align() was called"))

    // add the alignment.  NB: https://github.com/nh13/ksw/issues/6
    kswInput.write(query)
    kswInput.append('\n')
    kswInput.write(target)
    kswInput.append('\n')

    // update the previous alignment with a new one
    previousAlignment = previousAlignment.map { _ =>
      val alignment = kswOutput.readLine() match {
        case null => throw new IllegalStateException("KSW error.")
        case line => AsyncKswAligner.toAlignment(line)
      }
      aligned.getAndIncrement()
      alignment
    }
    previousAlignment
  }

  /** The number of completed alignments.  Returns -1 if unknown. */
  override def numAligned: Long = aligned.get()

  /** Closes all the resource related to the running ksw instance. */
  override def close(): Unit = this.close()

  /** Closes all the resource related to the running ksw instance. */
  def closeWithDuration(duration: Duration = Duration.Inf): Unit = {
    Await.ready(this.previousAlignment, duration)
    this.closed = true
    this.kswInput.close()
    this.kswOutput.close()
    if (this.process.isAlive) this.process.destroy()
  }
}