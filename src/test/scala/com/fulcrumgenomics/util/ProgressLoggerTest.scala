/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import java.io.{ByteArrayOutputStream, PrintStream}
import java.nio.charset.StandardCharsets

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.vcf.api.Variant
import org.scalatest.concurrent.PatienceConfiguration.Interval

class ProgressLoggerTest extends UnitSpec {

  private class LoggerHelper extends Logger(this.getClass) {
    private val baos = new ByteArrayOutputStream()
    out = Some(new PrintStream(baos, true, "UTF-8"))
    def lines: IndexedSeq[String] = new String(baos.toByteArray, StandardCharsets.UTF_8).split('\n').toIndexedSeq
  }
  
  // For Scala 2.12 compatibility
  private def emptyIterator[T]: Iterator[T] = Iterator.empty

  private val progressLogger = ProgressLogger(new LoggerHelper())

  "ProgressLoggingIterator" should "wrap a SamRecord, (String, Int), Variant, and Interval" in {
    import com.fulcrumgenomics.util.ProgressLogger.ProgressLoggingIterator

    // Check typing
    emptyIterator[SamRecord].progress(progressLogger)
    emptyIterator[(String, Int)].progress(progressLogger)
    emptyIterator[Variant].progress(progressLogger)
    emptyIterator[Interval].progress(progressLogger)

    // Do an actual test
    val logger = new LoggerHelper()
    val progress = ProgressLogger(logger, unit=2)
    Iterator(("chr1", 1), ("chr2", 2), ("chr3", 3)).progress(progress).foreach(_ => ())
    val lines = logger.lines
    lines.length shouldBe 2
    lines(0) should include("chr2:2")
    lines(1) should include("chr3:3")
  }

  it should "wrap unsupported types" in {
    import com.fulcrumgenomics.util.ProgressLogger.ProgressLoggingIterator

    // Check typing
    emptyIterator[Double].progress(progressLogger)
    emptyIterator[String].progress(progressLogger)

    // Do an actual test
    val logger = new LoggerHelper()
    val progress = ProgressLogger(logger, unit=2)
    Iterator("foo", "bar", "car").progress(progress).foreach(_ => ())
    val lines = logger.lines
    lines.length shouldBe 2
    lines(0) should include("*/*")
    lines(1) should include("*/*")
  }

  "TransformedProgressLoggingIterator" should "convert items to a supported type" in {
    import com.fulcrumgenomics.util.ProgressLogger.TransformedProgressLoggingIterator
    emptyIterator[(String, String)].progress(progressLogger, (item: (String, String)) => (item._1, item._2.toInt))

    // Do an actual test
    val logger = new LoggerHelper()
    val progress = ProgressLogger(logger, unit=2)
    Iterator(("chr1", "1"), ("chr2", "2"), ("chr3", "3")).progress(progress, (item: (String, String)) => (item._1, item._2.toInt)).foreach(_ => ())
    val lines = logger.lines
    lines.length shouldBe 2
    lines(0) should include("chr2:2")
    lines(1) should include("chr3:3")
  }
}
