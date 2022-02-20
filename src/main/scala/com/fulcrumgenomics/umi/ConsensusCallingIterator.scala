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
import com.fulcrumgenomics.bam.api
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.async.AsyncIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.commons.util.Threads.IterableThreadLocal
import com.fulcrumgenomics.umi.UmiConsensusCaller.SimpleRead
import com.fulcrumgenomics.util.ProgressLogger

import java.util.concurrent.ForkJoinPool

/**
  * An iterator that consumes from an incoming iterator of [[SamRecord]]s and generates consensus
  * read [[SamRecord]]s using the supplied consensus caller.
  *
  * @param sourceIterator the iterator over input [[SamRecord]]s.
  * @param caller the consensus caller to use to call consensus reads
  * @param progress an optional progress logger to which to log progress in input reads
  * @param threads the number of threads to use.
  * @param maxRecordsInRam the approximate maximum number of input records to store in RAM across multiple threads.
  */
class ConsensusCallingIterator[ConsensusRead <: SimpleRead](sourceIterator: Iterator[SamRecord],
                                                            caller: UmiConsensusCaller[ConsensusRead],
                                                            progress: Option[ProgressLogger] = None,
                                                            threads: Int = 1,
                                                            maxRecordsInRam: Int = 128000)
  extends Iterator[SamRecord] with LazyLogging {

  private val progressIterator = progress match {
    case Some(p) => sourceIterator.tapEach { r => p.record(r) }
    case None    => sourceIterator
  }

  protected val iter: Iterator[SamRecord] = {
    if (threads <= 1) {
      val groupingIterator = new SamRecordGroupedIterator(progressIterator, caller.sourceMoleculeId)
      groupingIterator.flatMap(caller.consensusReadsFromSamRecords)
    }
    else {
      val halfMaxRecords   = maxRecordsInRam / 2
      val pool             = new ForkJoinPool(threads - 1, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
      val callers          = new IterableThreadLocal[UmiConsensusCaller[ConsensusRead]](() => caller.emptyClone())
      val groupingIterator = {
        val async   = AsyncIterator(progressIterator, Some(halfMaxRecords))
        val grouped = new SamRecordGroupedIterator(async, caller.sourceMoleculeId)
        grouped.bufferBetter
      }

      // Create an iterator that will pull in input records (up to the half the maximum number in memory)
      // to be consensus called.  Each chunk of records will then be called in parallel.
      new Iterator[Seq[SamRecord]] {
        private var statisticsCollected = false

        override def hasNext: Boolean = {
          if (groupingIterator.hasNext) true else {
            // If we've hit the end then aggregate statistics
            if (!statisticsCollected) {
              callers.foreach(caller.addStatistics)
              statisticsCollected = true
            }
            false
          }
        }

        override def next(): Seq[SamRecord] = {
          var total = 0L
          groupingIterator
            .takeWhile { chunk => if (halfMaxRecords <= total) false else {total += chunk.length; true } }
            .toSeq
            .parWith(pool)
            .flatMap { records =>
              val caller = callers.get()
              caller.synchronized { caller.consensusReadsFromSamRecords(records) }
            }
            .seq
        }
      }.flatten
    }
  }
  override def hasNext: Boolean = this.iter.hasNext
  override def next(): SamRecord = this.iter.next
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
