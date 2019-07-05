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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.async.AsyncIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.umi.UmiConsensusCaller.SimpleRead
import com.fulcrumgenomics.util.ProgressLogger

import scala.collection.mutable.ListBuffer
import scala.concurrent.forkjoin.ForkJoinPool

/**
  * An iterator that consumes from an incoming iterator of [[SamRecord]]s and generates consensus
  * read [[SamRecord]] using the supplied consensus caller.
  *
  * @param sourceIterator the iterator over input [[SamRecord]]s.
  * @param caller the consensus caller to use to call consensus reads
  * @param progress an optional progress logger to which to log progress in input reads
  * @param threads the number of threads to use.
  * @param maxRecordsInRam the approximate maximum number of input records to store in RAM across multiple threads.
  */
class ConsensusCallingIterator[C <: SimpleRead](sourceIterator: Iterator[SamRecord],
                               caller: UmiConsensusCaller[C],
                               progress: Option[ProgressLogger] = None,
                               threads: Int = 1,
                               maxRecordsInRam: Int = 128000)
  extends Iterator[SamRecord] with LazyLogging {
  protected val groupingIterator: Iterator[Seq[SamRecord]] = {
    new SamRecordGroupedIterator(sourceIterator, caller.sourceMoleculeId)
      .map { records => for (p <- progress; r <- records) p.record(r); records }
  }
  protected val iterator: Iterator[SamRecord] = {
    if (threads <= 1) {
      groupingIterator.flatMap(caller.consensusReadsFromSamRecords)
    }
    else {
      val pool          = new ForkJoinPool(threads, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
      val bufferedIter  = groupingIterator.bufferBetter
      val callersFactor = 1000 // the number of callers per thread to create. This to reduce chance of picking the same caller
      val callers       = Array.range(start=0, threads*callersFactor).map { _ => caller.emptyClone() }

      // Read in groups of records (each from the same source molecule) until we have the maximum number of
      // individual records in RAM.  We have a few more records in RAM, if the group that pushes us over the limit is
      // large.  Then process the collected groups in parallel.
      Iterator.continually {
        var total = 0L
        bufferedIter
          .takeWhile { chunk => if (maxRecordsInRam <= total) false else {total += chunk.length; true } }
          .zipWithIndex
          .toSeq
          .parWith(pool)
          .flatMap { case (records, groupIndex) =>
            val caller = callers(groupIndex % callers.length)
            caller.synchronized { caller.consensusReadsFromSamRecords(records) }
          }
          .seq
      }.takeWhile { records =>
        if (records.nonEmpty) true else {
          // add the statistics to the original caller since there are no more reads
          require(bufferedIter.isEmpty, "Bug: input is not empty")
          callers.foreach(caller.addStatistics)
          false
        }
      }.flatten

    }
  }
  override def hasNext: Boolean = this.iterator.hasNext
  override def next(): SamRecord = this.iterator.next
}

/** Groups consecutive records based on a method to group records. */
private class SamRecordGroupedIterator[Key](sourceIterator: Iterator[SamRecord],
                                            toKey: SamRecord => Key) extends Iterator[Seq[SamRecord]] {
  private val input     = sourceIterator.bufferBetter
  private var nextChunk = IndexedSeq.empty[SamRecord]

  /** True if there are more consensus reads, false otherwise. */
  def hasNext(): Boolean = this.nextChunk.nonEmpty || (this.input.nonEmpty && advance())

  /** Returns the next consensus read. */
  def next(): Seq[SamRecord] = {
    if (!this.hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false.")
    yieldAndThen { nextChunk } { nextChunk = IndexedSeq.empty }
  }

  /** Consumes the next group of records from the input iterator, based on the vkey, and returns them as a [[IndexedSeq]]. */
  private def advance(): Boolean = {
    this.input.headOption.exists { head =>
      val idToMatch = this.toKey(head)
      this.nextChunk = this.input.takeWhile(this.toKey(_) == idToMatch).toIndexedSeq
      true
    }
  }
}