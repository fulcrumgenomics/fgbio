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
 *
 */

package com.fulcrumgenomics.vcf

import java.util.NoSuchElementException

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util._
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.annotation.tailrec

object ByIntervalListVariantContextIterator {

  /**
    * Creates an iterator over variant contexts that overlap any interval in an interval list.
    *
    * All variants will be read in and compared to the intervals.
    */
  def apply(iterator: Iterator[VariantContext],
            intervalList: IntervalList,
            dict: SequenceDictionary): Iterator[VariantContext] = {
    new OverlapDetectionVariantContextIterator(iterator, intervalList, dict)
  }

  /**
    * Creates an iterator over variant contexts that overlap any interval in an interval list.
    *
    * The interval lists should be sorted and uniqued.
    *
    * The VCF will be queried when moving to the next variant context, and so may be quite slow
    * if we jump around the VCF a lot.
    */
  def apply(reader: VCFFileReader,
            intervalList: IntervalList): Iterator[VariantContext] = {
    new IndexQueryVariantContextIterator(reader, intervalList)
  }
}

private class OverlapDetectionVariantContextIterator(val iter: Iterator[VariantContext],
                                                     val intervalList: IntervalList,
                                                     val dict: SequenceDictionary)
  extends Iterator[VariantContext] {

  require(dict != null)
  private val intervals = intervalList.uniqued(false).iterator().buffered
  private var nextVariantContext: Option[VariantContext] = None

  this.advance()

  def hasNext(): Boolean = this.nextVariantContext.isDefined

  def next(): VariantContext = {
    this.nextVariantContext match {
      case None      => throw new NoSuchElementException("Called next when hasNext is false")
      case Some(ctx) => yieldAndThen(ctx) { this.nextVariantContext = None; this.advance() }
    }
  }

  def remove(): Unit = throw new UnsupportedOperationException

  private def contigIdx(locatable: Locatable): Int = dict(locatable.getContig).index

  private def coordLessThan(left: Locatable, right: Locatable): Boolean = {
    val leftContigIdx = contigIdx(left)
    val rightContigIdx = contigIdx(right)
    leftContigIdx < rightContigIdx || (leftContigIdx == rightContigIdx && left.getEnd < right.getStart)
  }

  private def overlaps(left: Locatable, right: Locatable): Boolean = {
    contigIdx(left) == contigIdx(right) && CoordMath.overlaps(left.getStart, left.getEnd, right.getStart, right.getEnd)
  }

  @tailrec
  private def advance(): Unit = {
    if (this.iter.isEmpty) return

    val ctx = iter.next()

    // Move to the interval that overlaps or is past this context...
    while (intervals.hasNext && coordLessThan(intervals.head, ctx)) {
      intervals.next()
    }

    if (intervals.isEmpty) { } // no more intervals
    else if (overlaps(ctx, intervals.head)) { nextVariantContext = Some(ctx) } // overlaps!!!
    else if (iter.isEmpty) { } // no more variants
    else { this.advance() } // move to the next context
  }
}

/** NB: if a variant overlaps multiple intervals, only returns it once. */
private class IndexQueryVariantContextIterator(private val reader: VCFFileReader, intervalList: IntervalList)
  extends Iterator[VariantContext] {
  private var iter: Option[Iterator[VariantContext]] = None
  private val intervals = intervalList.iterator()
  private var previousInterval: Option[Interval] = None

  this.advance()

  def hasNext(): Boolean = {
    this.iter.exists(_.hasNext)
  }

  def next(): VariantContext = {
    this.iter match {
      case None => throw new NoSuchElementException("Called next when hasNext is false")
      case Some(i) => yieldAndThen(i.next()) { advance() }
    }
  }

  def remove(): Unit = throw new UnsupportedOperationException

  private def overlapsInterval(ctx: VariantContext, interval: Interval): Boolean = {
    if (!ctx.getContig.equals(interval.getContig)) false // different contig
    else if (interval.getStart <= ctx.getStart && ctx.getStart <= interval.getEnd) true // start falls within this interval, count it
    else if (ctx.getStart < interval.getStart && interval.getEnd <= ctx.getEnd) true // the variant encloses the interval
    else false
  }

  private def advance(): Unit = {
    while (!this.iter.exists(_.hasNext) && this.intervals.hasNext) {
      val interval = this.intervals.next()

      // Validate sort order of the intervals
      previousInterval.foreach { lastInterval =>
        val lastIntervalIdx = intervalList.getHeader.getSequenceIndex(lastInterval.getContig)
        val intervalIdx     = intervalList.getHeader.getSequenceIndex(interval.getContig)
        if (intervalIdx < lastIntervalIdx || (lastIntervalIdx == intervalIdx && interval.getStart <= lastInterval.getEnd)) {
          throw new IllegalStateException(s"Intervals are out of order: '$lastInterval' and '$lastInterval'")
        }
      }

      val lastInterval = previousInterval
      previousInterval = Some(interval)

      // NB: for variants that span an indel, make sure it was not output in the previous interval
      val iter = this.reader.query(interval.getContig, interval.getStart, interval.getEnd)
        .filter { ctx => overlapsInterval(ctx, interval) }
        .filter { ctx => lastInterval.forall { interval => !overlapsInterval(ctx, interval) } }

      this.iter = Some(iter)
    }
  }
}
