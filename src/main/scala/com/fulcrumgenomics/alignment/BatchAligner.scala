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
import com.fulcrumgenomics.commons.async.AsyncSink
import com.fulcrumgenomics.commons.util.StringUtil

/** The set of objects, traits, and classes are to facilitate the batching of alignment tasks.
  *
  * In cases where we wish to perform a large batch of asynchronous alignment, we can provide in one or more alignment
  * tasks (ideally many) to the batch aligner, and have alignments returned in the same order as added.
  *
  * An alignment task is provided as an [[AlignmentTask]], with a query and target respectively.  The query (or target)
  * is an [[Alignable]], which provides a wrapper around an object to provide methods to retrieve the bases to align.
  * Rather than just providing the bases to align, this allows the "alignable" query (or target) object to be returned
  * as part of the alignment.
  *
  * The batched aligners implement the [[BatchAligner]] trait.  This trait provides asychronous batching of inputs to the
  * underlying aligner.  Currently, there are two implementations:
  * 1. [[ScalaBatchAligner]] - the default implementation based on [[Aligner]] which will align immediately when an a
  *    alignment task is added.
  * 2. [[KswBatchAligner]] - an aligner using the `ksw` tool as an external executable.
  */


/** Companion to the [[Alignable]] class. */
object Alignable {
  /** Converts a [[String]] to an [[Alignable]]. */
  implicit def stringToAlignable(str: String): Alignable = new Alignable {
    override def bases: Array[Byte] = str.getBytes
  }
}

/** A trait that all classes that are alignable should extend. */
trait Alignable {
  def bases: Array[Byte]
  def length: Int = bases.length
}


object AlignmentTask {
  /** Converts two [[Alignable]]s to an [[AlignmentTask]]. */
  implicit def queryAndTargetToAlignmentTask[A <: Alignable, B <: Alignable](query: A, target: B): AlignmentTask[A, B] = new AlignmentTask(query, target)
}

/** An input to the aligner.  The inputs implement the [[Alignable]] trait to allow for deferred extraction/computation
  * of the query and target respectively.  This may be useful if aligning many sub-sequences of a larger query to a
  * variety of targets, thus reducing memory. */
case class AlignmentTask[A <: Alignable, B <: Alignable](query: A, target: B) {
  def queryBytes: Array[Byte]  = query.bases
  def targetBytes: Array[Byte] = target.bases
  def queryLength: Int         = query.length
  def targetLength: Int        = target.length
}

object BatchAligner {
  // TODO: finish docs
  /** Constructs a new interactive aligner. Will use the [[KswBatchAligner]] implementation if the path to the `ksw` executable is given.
    *
    * @param matchScore the match score (should be greater than or equal to zero).
    * @param mismatchScore the mismatch score (should be greater than or equal to zero).
    * @param gapOpen the gap opening penalty, should generally be negative or zero
    * @param gapExtend the gap extension penalty, should generally be negative or zero
    * @param mode alignment mode to use when generating alignments
    * @param ksw
    * @tparam A
    * @tparam B
    * @return
    */
  def apply[A <: Alignable, B <: Alignable](matchScore: Int,
                                            mismatchScore: Int,
                                            gapOpen: Int,
                                            gapExtend: Int,
                                            mode: Mode = Mode.Glocal,
                                            ksw: Option[FilePath] = None): BatchAligner[A, B] = {
    ksw match {
      case Some(executable) => new KswBatchAligner(executable, matchScore, mismatchScore, gapOpen, gapExtend, mode)
      case None             => ScalaBatchAligner(matchScore, mismatchScore, gapOpen, gapExtend, mode)
    }
  }
}

/** A trait that all aligners that can process alignments tasks in batches should implement.
  *
  * The aligner should be able to handle more than one input alignment task ([[AlignmentTask]]), aligning them
  * asynchronously from the task of adding or retrieving alignments.
  *
  * The alignment results are returned in the same order as alignment tasks.
  *
  * @tparam A the type of the query [[Alignable]].
  * @tparam B the type of the target [[Alignable]].
  */
trait BatchAligner[A <: Alignable, B <: Alignable] extends Iterator[Alignment] with Closeable {

  // TODO: use a cache...

  /** The number of alignments added so far. */
  private var _numAdded: Long = 0

  /** The number of alignment tasks returned so far. */
  private var _numRetrieved: Long = 0

  /** Provides an asynchronous way of adding alignment tasks to the aligner*/
  private val sink = new AsyncSink[AlignmentTask[A, B]](
    sink       = t => this._append(t),
    source     = None
  )

  /** Adds one or more alignment tasks to this aligner. */
  final def append(task: AlignmentTask[A, B]*): Unit = {
    task.foreach(this.sink.add)
    this._numAdded += 1
  }

   /** Gets the result of the next alignment task. */
  final def next(): Alignment = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false.")
    this._numRetrieved += 1
    this._next()
  }

  /** True if there is an outstanding alignment task, false otherwise. */
  final def hasNext(): Boolean = this._numRetrieved < this._numAdded

  /** Gets an iterator over the alignment tasks*/
  def iterator: Iterator[Alignment] = this

  /** The number of alignment tasks added so far. */
  def numAdded: Long = _numAdded

  /** The number of alignment tasks returned so far. */
  def numRetrieved: Long = _numRetrieved

  /** The number of alignment tasks */
  def numAvailable: Long = numAdded - numRetrieved

  /** Reset the various counts.  There must be no tasks remaining. */
  def reset(): Unit = {
    require(!this.hasNext(), "Calling rest() when hasNext() is true")
    this._numAdded = 0
    this._numRetrieved = 0
  }

  /** All alignment implementations must implement this to add the task to the aligner. */
  protected def _append(task: AlignmentTask[A, B]): Unit

  /** All alignment implementations must implement this to return the next alignment result. */
  protected def _next(): Alignment
}

/** Companion to [[ScalaBatchAligner]]. */
object ScalaBatchAligner {
  /** Creates a new [[ScalaBatchAligner]].
    *
    * @param matchScore the match score (should be greater than or equal to zero).
    * @param mismatchScore the mismatch score (should be greater than or equal to zero).
    * @param gapOpen the gap opening penalty, should generally be negative or zero
    * @param gapExtend the gap extension penalty, should generally be negative or zero
    * @param mode alignment mode to use when generating alignments
    * @tparam A the type of the query [[Alignable]].
    * @tparam B the type of the target [[Alignable]].
    */
  def apply[A <: Alignable, B <: Alignable](matchScore: Int,
            mismatchScore: Int,
            gapOpen: Int,
            gapExtend: Int,
            mode: Mode = Mode.Glocal): ScalaBatchAligner[A, B] = {
    val aligner: Aligner = Aligner(matchScore = matchScore, mismatchScore = mismatchScore, gapOpen = gapOpen, gapExtend = gapExtend, mode = mode)
    new ScalaBatchAligner(aligner)
  }
}

/** An aligner that batches alignment tasks using and aligns with [[Aligner]].
  *
  * This implementation will eagerly produces alignments, meaning it will align immediately when the alignment is added.
  **/
class ScalaBatchAligner[A <: Alignable, B <: Alignable](val aligner: Aligner) extends BatchAligner[A, B] {

  /** The queue of alignment tasks. */
  private val queue = new java.util.concurrent.LinkedBlockingQueue[Alignment]()

  /** Adds the task to the task queue. */
  protected def _append(task: AlignmentTask[A, B]): Unit = {
    queue.put(aligner.align(task.queryBytes, task.targetBytes))
  }

  /** Does the next enqueued alignment task and returns its result.*/
  protected def _next(): Alignment = queue.take()

  /** Does nothing */
  def close(): Unit = Unit
}


/** Companion class to [[KswBatchAligner]]. */
private object KswBatchAligner {

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


/** An batched aligner that wraps the [ksw](https://github.com/nh13/ksw) executable for faster alignments.
  *
  * Install [the ksw executable version 0.2.0](https://github.com/nh13/ksw) manually or via conda
  * (`conda install -c bioconda ksw=0.2.0`)
  *
  * @param executable the path to the `ksw` executable.
  * @param matchScore the match score (should be greater than or equal to zero).
  * @param mismatchScore the mismatch score (should be greater than or equal to zero).
  * @param gapOpen the gap opening penalty, should generally be negative or zero
  * @param gapExtend the gap extension penalty, should generally be negative or zero
  * @param mode alignment mode to use when generating alignments
  * @param buffer the size of the input and output buffer in bytes.
  * @tparam A the type of the query [[Alignable]].
  * @tparam B the type of the target [[Alignable]].
  */
private class KswBatchAligner[A <: Alignable, B <: Alignable](executable: FilePath,
                                                              matchScore: Int,
                                                              mismatchScore: Int,
                                                              gapOpen: Int,
                                                              gapExtend: Int,
                                                              mode: Mode = Mode.Glocal,
                                                              buffer: Int = 1024*1024) extends BatchAligner[A, B]
{
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

  /** Adds an alignment task to Ksw's input.  This may block writing the input to ksw. */
  protected def _append(task: AlignmentTask[A, B]): Unit = {
    // Align the forward and reverse: https://github.com/nh13/ksw/issues/6
    kswInput.write(task.queryBytes)
    kswInput.append('\n')
    kswInput.write(task.targetBytes)
    kswInput.append('\n')
  }

  /** Retrieves the next result from ksw.  This may block reading the output of ksw. */
  protected def _next(): Alignment = {
    kswOutput.readLine() match {
      case null => throw new IllegalStateException("KSW error.")
      case line => KswBatchAligner.toAlignment(line)
    }
  }

  /** Closes all the resource related to the running Primer3 instance. */
  override def close(): Unit = {
    this.kswInput.close()
    this.kswOutput.close()
    if (this.process.isAlive) this.process.destroy()
  }
}

