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
import com.fulcrumgenomics.bam.api.{SamRecord, TransientAttrs}
import com.fulcrumgenomics.umi.UmiConsensusCaller.SimpleRead
import com.fulcrumgenomics.util.ProgressLogger

import scala.concurrent.forkjoin.ForkJoinPool

object ConsensusCallingIterator {

  type ConsensusCallingIterator = Iterator[SamRecord]

  /**
    * An iterator that consumes from an incoming iterator of [[SamRecord]]s and generates consensus
    * read [[SamRecord]] using the supplied consensus caller.
    *
    * @param sourceIterator the iterator over input [[SamRecord]]s.
    * @param caller the consensus caller to use to call consensus reads
    * @param progress an optional progress logger to which to log progress in input reads
    */
  def apply(sourceIterator: Iterator[SamRecord],
            caller: UmiConsensusCaller[_],
            progress: Option[ProgressLogger] = None
           ): ConsensusCallingIterator = {
    val groupingIterator = new SamRecordGroupedIterator(sourceIterator, caller.sourceMoleculeId)
    groupingIterator.flatMap { records =>
      for (p <- progress; r <- records) p.record(r)
      caller.consensusReadsFromSamRecords(records)
    }
  }

  /**
    * An iterator that consumes from an incoming iterator of [[SamRecord]]s and generates consensus
    * read [[SamRecord]] using the supplied consensus caller.
    *
    * @param sourceIterator the iterator over input [[SamRecord]]s.
    * @param toCaller method to build a consensus caller to use to call consensus reads
    * @param progress an optional progress logger to which to log progress in input reads
    * @param threads the number of threads to use.
    * @param maxRecordsInRam the maximum number of input records to store in RAM.  May be larger if a single consensus
    *                        read contains more records.
    */
  def parWith[T<:SimpleRead](sourceIterator: Iterator[SamRecord],
              toCaller: () => UmiConsensusCaller[T],
              progress: Option[ProgressLogger] = None,
              threads: Int = 1,
              maxRecordsInRam: Int = 128000
             ): (ConsensusCallingIterator, UmiConsensusCaller[T]) = {
    if (threads <= 1) {
      val caller = toCaller()
      (this.apply(sourceIterator, caller, progress), caller)
    }
    else {
      val callers          = IndexedSeq.range(start=0, end=threads).map(_ => toCaller())
      val groupingIterator = new SamRecordGroupedIterator(sourceIterator, callers.head.sourceMoleculeId)
      val pool             = new ForkJoinPool(threads, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
      val iter             = groupingIterator.bufferBetter
      val outIter          = Iterator.continually {
        // Collect sets of input reads, each set to call a consensus read, until we have consumed the specified number
        // of input records.  Then we can process that batch.
        var total = 0L
        iter
          .takeWhile { chunk => if (maxRecordsInRam <= total) false else { total += chunk.length; true } }
          .map { records => yieldAndThen(records) { for (p <- progress; r <- records) p.record(r) } }
          .zipWithIndex
          .toSeq
          .parWith(pool)
          .flatMap { case (records, idx) => callers(idx % threads).consensusReadsFromSamRecords(records) }
          .seq
      }.takeWhile { records =>
        if (records.nonEmpty) true else {
          // add the statistics
          callers.tail.foreach { caller => callers.head.addStatistics(caller) }
          false
        }
      }.flatten
      (outIter, callers.head)
    }
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
}