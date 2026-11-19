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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util._
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}

import scala.annotation.tailrec

object ByIntervalListVariantContextIterator {

  /**
    * Creates an iterator over variants that overlap any interval in an interval list.
    *
    * All variants will be read in and compared to the intervals.
    */
  def apply(iterator: Iterator[Variant],
            intervalList: IntervalList,
            dict: SequenceDictionary): Iterator[Variant] = {
    new OverlapDetectionVariantContextIterator(iterator, intervalList, dict)
  }

  /**
    * Creates an iterator over variants that start in, or enclose, any interval in an interval list.
    *
    * The interval lists should be sorted and uniqued.  The VCF must be indexed.
    *
    * The VCF will be queried when moving to the next interval, and so may be quite slow
    * if we jump around the VCF a lot.
    */
  def apply(source: VcfSource,
            intervalList: IntervalList): Iterator[Variant] = {
    new IndexQueryVariantContextIterator(source, intervalList)
  }
}

private class OverlapDetectionVariantContextIterator(val iter: Iterator[Variant],
                                                     val intervalList: IntervalList,
                                                     val dict: SequenceDictionary)
  extends Iterator[Variant] {

  require(dict != null)
  private val intervals = intervalList.uniqued(false).iterator().buffered
  private var nextVariant: Option[Variant] = None

  this.advance()

  def hasNext: Boolean = this.nextVariant.isDefined

  def next(): Variant = {
    this.nextVariant match {
      case None          => throw new NoSuchElementException("Called next when hasNext is false")
      case Some(variant) => yieldAndThen(variant) { this.nextVariant = None; this.advance() }
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

    val variant   = iter.next()
    val locatable = variant.toLocatable

    // Move to the interval that overlaps or is past this variant...
    while (intervals.hasNext && coordLessThan(intervals.head, locatable)) {
      intervals.next()
    }

    if (intervals.isEmpty) { } // no more intervals
    else if (overlaps(locatable, intervals.head)) { nextVariant = Some(variant) } // overlaps!!!
    else if (iter.isEmpty) { } // no more variants
    else { this.advance() } // move to the next variant
  }
}

/** NB: if a variant overlaps multiple intervals, only returns it once. */
private class IndexQueryVariantContextIterator(private val source: VcfSource, intervalList: IntervalList)
  extends Iterator[Variant] {
  private var iter: Option[Iterator[Variant]] = None
  private val intervals = intervalList.iterator()
  private var previousInterval: Option[Interval] = None

  this.advance()

  def hasNext: Boolean = {
    this.iter.exists(_.hasNext)
  }

  def next(): Variant = {
    this.iter match {
      case None => throw new NoSuchElementException("Called next when hasNext is false")
      case Some(i) => yieldAndThen(i.next()) { advance() }
    }
  }

  def remove(): Unit = throw new UnsupportedOperationException

  private def overlapsInterval(variant: Variant, interval: Interval): Boolean = {
    if (variant.chrom != interval.getContig) false // different contig
    else if (interval.getStart <= variant.pos && variant.pos <= interval.getEnd) true // start falls within this interval, count it
    else if (variant.pos < interval.getStart && interval.getEnd <= variant.end) true // the variant encloses the interval
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
      val iter = this.source.query(interval.getContig, interval.getStart, interval.getEnd)
        .filter { variant => overlapsInterval(variant, interval) }
        .filter { variant => lastInterval.forall { interval => !overlapsInterval(variant, interval) } }

      this.iter = Some(iter)
    }
  }
}
