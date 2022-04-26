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
import com.fulcrumgenomics.commons.collection.ParIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.commons.util.Threads.IterableThreadLocal
import com.fulcrumgenomics.umi.UmiConsensusCaller.SimpleRead
import com.fulcrumgenomics.util.ProgressLogger

/**
  * An iterator that consumes from an incoming iterator of [[SamRecord]]s and generates consensus
  * read [[SamRecord]]s using the supplied consensus caller.
  *
  * @param sourceIterator the iterator over input [[SamRecord]]s.
  * @param caller the consensus caller to use to call consensus reads
  * @param progress an optional progress logger to which to log progress in input reads
  * @param threads the number of threads to use.
  * @param chunkSize parallel process in chunkSize units; will cause 8 * chunkSize records to be held in memory
  */
class ConsensusCallingIterator[ConsensusRead <: SimpleRead](sourceIterator: Iterator[SamRecord],
                                                            caller: UmiConsensusCaller[ConsensusRead],
                                                            progress: Option[ProgressLogger] = None,
                                                            threads: Int = 1,
                                                            chunkSize: Int = ParIterator.DefaultChunkSize)
  extends Iterator[SamRecord] with LazyLogging {

  private val callers = new IterableThreadLocal[UmiConsensusCaller[ConsensusRead]](() => caller.emptyClone())
  private var collectedStats: Boolean = false

  protected val iter: Iterator[SamRecord] = {
    // Wrap our input iterator in a progress logging iterator if we have a progress logger
    val progressIterator = progress match {
      case Some(p) => sourceIterator.tapEach { r => p.record(r) }
      case None    => sourceIterator
    }

    // Then turn it into a grouping iterator
    val groupingIterator = new SamRecordGroupedIterator(progressIterator, caller.sourceMoleculeId)

    // Then call consensus either single-threaded or multi-threaded
    if (threads <= 1) {
      groupingIterator.flatMap(caller.consensusReadsFromSamRecords)
    }
    else {
      ParIterator(groupingIterator, threads=threads).flatMap { rs =>
        val caller = callers.get()
        caller.synchronized { caller.consensusReadsFromSamRecords(rs) }
      }.toAsync(chunkSize * 8)
    }
  }

  // Responsible for adding statistics to the main caller once we hit the end of the iterator.
  override def hasNext: Boolean = {
    if (this.iter.hasNext) true else {
      if (!collectedStats) {
        callers.foreach(c => caller.addStatistics(c))
        collectedStats = true
      }

      false
    }
  }

  override def next(): SamRecord = this.iter.next()
}

/** Groups consecutive records based on a method to group records. */
private class SamRecordGroupedIterator[Key](sourceIterator: Iterator[SamRecord],
                                            toKey: SamRecord => Key) extends Iterator[Seq[SamRecord]] {
  private val input     = sourceIterator.bufferBetter
  private var nextChunk = IndexedSeq.empty[SamRecord]

  /** True if there are more consensus reads, false otherwise. */
  def hasNext: Boolean = this.nextChunk.nonEmpty || (this.input.nonEmpty && advance())

  /** Returns the next consensus read. */
  def next(): Seq[SamRecord] = {
    if (!this.hasNext) throw new NoSuchElementException("Calling next() when hasNext() is false.")
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
