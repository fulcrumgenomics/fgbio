/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.SAMRecord

import scala.collection.mutable
import scala.collection.mutable.ListBuffer
import com.fulcrumgenomics.FgBioDef._

/**
  * An iterator that consumes from an incoming iterator of SAMRecords and generates consensus
  * read SAMRecords using the supplied caller.
  *
  * @param sourceIterator an iterator of SAMRecords
  * @param caller   the consensus caller to use to call consensus reads
  * @param progress an optional progress logger to which to log progress in input reads
  */
class ConsensusCallingIterator(sourceIterator: Iterator[SAMRecord],
                               val caller: UmiConsensusCaller[_],
                               val progress: Option[ProgressLogger] = None
                              ) extends Iterator[SAMRecord] {

  private val input = sourceIterator.bufferBetter
  private val outputQueue: mutable.Queue[SAMRecord] = mutable.Queue[SAMRecord]()

    /** True if there are more consensus reads, false otherwise. */
  def hasNext(): Boolean = this.outputQueue.nonEmpty || (this.input.nonEmpty && advance())

  /** Returns the next consensus read. */
  def next(): SAMRecord = {
    if (!this.hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false.")
    this.outputQueue.dequeue()
  }

  /**
    * Consumes the next group of records from the input iterator, based on molecule id
    * and returns them as a Seq.
    */
  private def nextGroupOfRecords(): Seq[SAMRecord] = {
    if (this.input.isEmpty) {
      Nil
    }
    else {
      val idToMatch = this.caller.sourceMoleculeId(this.input.head)
      this.input.takeWhile(this.caller.sourceMoleculeId(_) == idToMatch).toSeq
    }
  }

  /** Consumes input records until one or more consensus reads can be created, or no more input records are available.
    * Returns true if a consensus read was created and enqueued, false otherwise. */
  @annotation.tailrec
  private def advance(): Boolean = {
    // get the records to create the consensus read
    val inputs  = nextGroupOfRecords()
    val outputs = this.caller.consensusReadsFromSamRecords(inputs)
    this.outputQueue ++= outputs

    // Log progress on the _input_ reads and then return/recurse
    for (p <- progress; r <- inputs) p.record(r)
    if (outputs.nonEmpty) true
    else if (this.input.hasNext) advance()
    else false
  }
}
