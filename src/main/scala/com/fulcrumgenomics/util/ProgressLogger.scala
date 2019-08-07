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
 */
package com.fulcrumgenomics.util

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.vcf.api.Variant
import htsjdk.samtools.util.{AbstractProgressLogger, Locatable}

import scala.reflect.runtime.universe._

/**
  * A subclass of HTSJDK's progress logger that uses fgbio's logging system.
  */
case class ProgressLogger(logger: Logger,
                          noun: String = "records",
                          verb: String = "processed",
                          unit: Int = 1000 * 1000) extends AbstractProgressLogger(noun, verb, unit) {

  override def log(message: String*): Unit = {
    this.logger.info(message:_*)
  }

  /** Method to use to record progress when genome location isn't available or relevant. */
  def record(): Boolean = record(null, -1)

  def record(rec: SamRecord): Boolean = super.record(rec.asSam)

  def record(variant: Variant): Boolean = record(variant.chrom, variant.pos)

  def record(locus: Locatable): Boolean = record(locus.getContig, locus.getStart)

  /** Logs the last record if it wasn't already logged. */
  def logLast(): Boolean = super.log() // Calls the super's log() method, not log(message: String*)
}

object ProgressLogger {

  /** Provides a method to record items of the given type.  If [[ItemType]] is not one of [[SamRecord]],
    * `(String, Int)`, [[Variant]], or [[Locatable], [[ProgressLogger.record()]] will be used. */
  sealed abstract class ProgressLoggingHelper[ItemType: TypeTag] {
    protected val record: (ProgressLogger, ItemType) => Boolean = {
      typeOf[ItemType] match {
        case t if t <:< typeOf[SamRecord]     => (p: ProgressLogger, item: ItemType) => p.record(item.asInstanceOf[SamRecord])
        case t if t <:< typeOf[(String, Int)] => (p: ProgressLogger, item: ItemType) => { val (contig, pos) = item.asInstanceOf[(String, Int)]; p.record(contig, pos) }
        case t if t <:< typeOf[Variant]       => (p: ProgressLogger, item: ItemType) => p.record(item.asInstanceOf[Variant])
        case t if t <:< typeOf[Locatable]     => (p: ProgressLogger, item: ItemType) => p.record(item.asInstanceOf[Locatable])
        case _                                => (p: ProgressLogger, _: ItemType)    => p.record()
      }
    }
  }

  /** Wraps an [[Iterator]] and provides a method to record progress as items are consumed. */
  implicit class ProgressLoggingIterator[ItemType: TypeTag](iterator: Iterator[ItemType]) extends ProgressLoggingHelper[ItemType] {
    def progress(progressLogger: ProgressLogger): Iterator[ItemType] = {
      new SelfClosingIterator(iterator.map { item => this.record(progressLogger, item); item }, () => progressLogger.logLast())
    }
  }

  /** Wraps an [[Iterator]] over [[ItemType]]s and provides a method to transform them into the [[ResultType]] to record progress as items are consumed. */
  implicit class TransformedProgressLoggingIterator[ItemType, ResultType: TypeTag](iterator: Iterator[ItemType]) extends ProgressLoggingHelper[ResultType] {
    def progress(progressLogger: ProgressLogger, transform: ItemType => ResultType): Iterator[ItemType] = {
      new SelfClosingIterator(iterator.map { item => this.record(progressLogger, transform(item)); item }, () => progressLogger.logLast())
    }
  }
}
