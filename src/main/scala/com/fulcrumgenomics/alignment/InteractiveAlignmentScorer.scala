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


import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.util.Io

import scala.collection.mutable

object Alignable {
  implicit def stringToAlignable(str: String): Alignable = new Alignable {
    override def bases: Array[Byte] = str.getBytes
  }
}

/** A trait that all classes that are alignable should extend. */
trait Alignable {
  def bases: Array[Byte]
  def length: Int = bases.length
}

/** An input to the aligner.  The inputs implement the [[Alignable]] trait to allow for deferred extraction/computation
  * of the query and target respectively.  This may be useful if aligning many sub-sequences of a larger query to a
  * variety of targets, thus reducing memory. */
case class AlignmentTask[A <: Alignable, B <: Alignable](query: A, target: B) {
  def toQuery: Array[Byte]  = query.bases
  def toTarget: Array[Byte] = target.bases // rec.bases.slice(targetStart, targetEnd)
  def queryLength: Int      = query.length
  def targetLength: Int     = target.length
}

object PartialAlignment {
  def apply(alignment: Alignment): PartialAlignment = {
    PartialAlignment(
      score       = alignment.score,
      queryStart  = alignment.queryStart,
      queryEnd    = alignment.queryEnd,
      targetStart = alignment.targetStart,
      targetEnd   = alignment.targetEnd
    )
  }
}

/** The result of a the aligner.  This does not give the full-alignment (i.e. [[Alignment]]) to reduce computation. */
case class PartialAlignment(score: Int, queryStart: Int, queryEnd: Int, targetStart: Int, targetEnd: Int)

object InteractiveAlignmentScorer {
  /** Constructs a new interactive aligner. Will use the [[Ksw]] implementation if the path to the `ksw` executable is given. */
  def apply[A <: Alignable, B <: Alignable](matchScore: Int,
                                            mismatchScore: Int,
                                            gapOpen: Int,
                                            gapExtend: Int,
                                            mode: Mode = Mode.Glocal,
                                            ksw: Option[FilePath] = None): InteractiveAlignmentScorer[A, B] = {
    ksw match {
      case Some(executable) => new Ksw(executable, matchScore, mismatchScore, gapOpen, gapExtend, mode)
      case None             => Scalaligner(matchScore, mismatchScore, gapOpen, gapExtend, mode)
    }
  }
}

/** A trait that all aligners who have the following properties should implement.
  *
  * 1. They are "interactive", meaning, new alignment tasks can be give to them and the results of a previously given
  * alignment task can be retrieved.  This is useful for batching large alignment queries, or for interacting with an
  * external program.
  * 2. Produces only a limited set of information about the best alignment:
  *   a. The maximum alignment score.
  *   b. Where the end of the alignment occurred in the query and target respectively.
  *
  * This aligner is meant to be used when the full-alignment is not needed but fast-alignment or batched is needed.
  *
  * The alignment results are returned in the same order as alignment tasks.
  * */
trait InteractiveAlignmentScorer[A <: Alignable, B <: Alignable] extends Iterator[PartialAlignment] with Closeable {
  /** The number of alignments added so far. */
  private var _numAdded: Long = 0
  /** The number of alignment tasks performed so far. */
  private var _numAligned: Long = 0

  /** Adds one or more alignment tasks to this alignerk*/
  final def append(alignmentInput: AlignmentTask[A, B]*): Unit = {
    this._append(alignmentInput:_*)
    this._numAdded += alignmentInput.size
  }

  /** Gets the result of the next alignment task. */
  final def next(): PartialAlignment = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false.")
    this._numAligned += 1
    this._next()
  }

  /** True if there is an outstanding alignment task, false otherwise. */
  final def hasNext(): Boolean = this._numAligned < this._numAdded

  /** Gets an iterator over the currently enqueued alignment tasks*/
  def iterator: Iterator[PartialAlignment] = {
    val numToRetrieve = (_numAdded - _numAligned).toInt
    if (numToRetrieve == 0) Iterator.empty
    else Range.inclusive(start = 0, end = numToRetrieve - 1, step = 1).iterator.map(_ => next())
  }

  /** The number of alignments added so far. */
  def numAdded: Long = _numAdded

  /** The number of alignment tasks performed so far. */
  def numAligned: Long = _numAligned

  /** All alignment implementations should implement this to add the task to the aligner. */
  protected def _append(alignmentInput: AlignmentTask[A, B]*): Unit

  /** All alignment implementations should implement this to return the next alignment result. */
  protected def _next(): PartialAlignment
}

/** Companion to [[Scalaligner]]. */
object Scalaligner {
  def apply[A <: Alignable, B <: Alignable](matchScore: Int,
            mismatchScore: Int,
            gapOpen: Int,
            gapExtend: Int,
            mode: Mode = Mode.Glocal): Scalaligner[A, B] = {
    val aligner: Aligner = Aligner(matchScore = matchScore, mismatchScore = mismatchScore, gapOpen = gapOpen, gapExtend = gapExtend, mode = mode)
    new Scalaligner(aligner)
  }
}

/** An interactive aligner that uses the scala [[Aligner]].
  *
  * This will be very slow for large numbers of queries.
  *
  * Lazy produces alignments.
  **/
class Scalaligner[A <: Alignable, B <: Alignable](val aligner: Aligner) extends InteractiveAlignmentScorer[A, B] {


  /** The queue of alignment tasks. */
  private val queue = new mutable.Queue[AlignmentTask[A, B]]()

  /** Adds the task to the task queue. */
  protected def _append(alignmentInput: AlignmentTask[A, B]*): Unit = queue.enqueue(alignmentInput:_*)

  /** Does the next enqueued alignment task and returns its result.*/
  protected def _next(): PartialAlignment = {
    val alignmentInput = queue.dequeue()
    val alignment      = aligner.align(alignmentInput.toQuery, alignmentInput.toTarget)
    PartialAlignment(alignment)
  }

  /** Does nothing */
  def close(): Unit = Unit
}


/** Companion class to [[Ksw]]. */
object Ksw {

  // Developer Note: this is faster than jus split.  See [[DelimitedDataParser]].
  /** Parses an alignment result from the ksw output (a single line). */
  private def toAlignmentResult(line: String, delimiter: Char = '\t'): PartialAlignment = {
    val tmp = new Array[Int](5)
    val cs = line.toCharArray
    val len = cs.length
    var i = 0
    var count = 0
    var start = 0
    var end = 0

    while (i <= len && count < tmp.length) {
      if (i == len || cs(i) == delimiter) {
        tmp(count) = new String(cs, start, end-start).toInt
        count += 1
        start = i + 1
        end   = i + 1
      }
      else {
        end += 1
      }

      i += 1
    }
    require(count == 5, s"Could not parse line into five values: $line.")
    PartialAlignment(tmp(0), tmp(1), tmp(2), tmp(3), tmp(4))
  }
}


/** An interactive aligner that uses attractivechaos's `ksw` aligner.
  *
  * The aligner is SIMD vectorized implementation of alignment extension, namely find the best alignment of prefixes of
  * the query and target.  No other modes are currently supported.
  *
  * The original code is here: https://github.com/attractivechaos/klib
  *
  * The modified code to work with this class is here: https://github.com/nh13/klib/tree/fgbio_ksw
  *
  * For the latter, use `gcc -o ksw -lz ksw.c main.c` to compile.
  */
private class Ksw[A <: Alignable, B <: Alignable](executable: FilePath,
                                                  matchScore: Int,
                                                  mismatchScore: Int,
                                                  gapOpen: Int,
                                                  gapExtend: Int,
                                                  mode: Mode = Mode.Glocal,
                                                  buffer: Int = 8*Io.bufferSize) extends InteractiveAlignmentScorer[A, B]
{
  /** The ksw process. */
  private val process   = {
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
      "-M", alignmentMode
    ).map(_.toString)
    new ProcessBuilder(args:_*).redirectErrorStream(true).start()
  }

  /** The input stream into ksw. */
  private val kswInput  = new PrintStream(new BufferedOutputStream(this.process.getOutputStream, buffer), true)

  /** The output stream out of ksw */
  private val kswOutput = new BufferedReader(new InputStreamReader(this.process.getInputStream), buffer)

  /** Adds one or more alignments to Ksw's input.  This may block writing the input to ksw. */
  protected def _append(alignmentInputs: AlignmentTask[A, B]*): Unit = {
    val stringBuilder = new StringBuilder(alignmentInputs.map(a => a.targetLength + a.queryLength + 1).sum)
    alignmentInputs.foreach { alignmentInput =>
      stringBuilder.append(alignmentInput.toQuery)
      stringBuilder.append('\n')
      stringBuilder.append(alignmentInput.toTarget)
      stringBuilder.append('\n')
    }
    kswInput.print(stringBuilder)
  }

  /** Retrievs the next result from ksw.  This may block reading the output of ksw. */
  protected def _next(): PartialAlignment = {
    kswOutput.readLine() match {
      case null => throw new IllegalStateException("KSW error.")
      case line => Ksw.toAlignmentResult(line)
    }
  }

  /** Closes all the resource related to the running Primer3 instance. */
  override def close(): Unit = {
    this.kswInput.close()
    this.kswOutput.close()
    if (this.process.isAlive) this.process.destroy()
  }
}

