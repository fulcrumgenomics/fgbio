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
import com.fulcrumgenomics.commons.util.Logger
import htsjdk.samtools.util.AbstractProgressLogger

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

  def record(rec: SamRecord): Boolean = record(rec.asSam)

  /** Logs the last record if it wasn't already logged. */
  def logLast(): Boolean = super.log() // Calls the super's log() method, not log(message: String*)

}
