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

package com.fulcrumgenomics.bam

import java.nio.file.{Path, Paths}

import com.fulcrumgenomics.testing.SamRecordSetBuilder.{Minus, Plus, Strand}
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric
import dagr.commons.io.PathUtil
import htsjdk.samtools.SAMFileHeader.SortOrder
import org.scalatest.ParallelTestExecution

/**
  * Tests for ErrorRateByReadPosition.
  */
class ErrorRateByReadPositionTest extends UnitSpec with ParallelTestExecution {
  val dir = Paths.get("src/test/resources/com/fulcrumgenomics/bam")
  private val referenceFasta = dir.resolve("error_rate_by_read_position.fasta")

  def outputAndPrefix: (Path, Path) = {
    val out = makeTempFile("ErrorRateByReadPositionTest", ErrorRateByReadPositionMetric.FileExtension)
    val pre = PathUtil.removeExtension(PathUtil.removeExtension(out))
    (out, pre)
  }

  "ErrorRateByReadPosition" should "compute the error rate for all unmapped reads" in {
    Seq(None, Some(referenceFasta)).foreach { case maybeRef =>
      val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
      builder.addFrag(unmapped = true)
      val input = builder.toTempFile()
      val (out, pre) = outputAndPrefix
      new ErrorRateByReadPosition(input = input, output = Some(pre), ref=maybeRef).execute()
      val metrics = Metric.read[ErrorRateByReadPositionMetric](out)
      metrics.length shouldBe 0
    }
  }

  it should "compute the error rate for paired reads" in {
    Seq(None, Some(referenceFasta)).foreach { case maybeRef =>
      val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
      builder.addPair(start1=1, start2=1).foreach { rec =>
        if (rec.getFirstOfPairFlag) {
          rec.setAttribute("MD", "100")
          rec.setReadString("A" * 100)
        }
        else {
          rec.setAttribute("MD", "A" * 100)
          rec.setReadString("C" * 100)
        }
      }
      builder.addPair(start1=1, start2=1).foreach { rec => rec.setAttribute("MD", "100"); rec.setReadString("A" * 100) }
      val input = builder.toTempFile()
      val (out, pre) = outputAndPrefix
      new ErrorRateByReadPosition(input=input, output=Some(pre), ref=maybeRef).execute()
      val metrics = Metric.read[ErrorRateByReadPositionMetric](out)
      metrics.length shouldBe 200
      metrics.foreach { case metric =>
        (metric.read_number == 1 || metric.read_number == 2) shouldBe true
        metric.read_number match {
          case 1 => metric shouldBe new ErrorRateByReadPositionMetric(1, metric.position, 0.0d, 2l)
          case 2 => metric shouldBe new ErrorRateByReadPositionMetric(2, metric.position, 0.5d, 2l)
        }
      }
    }
  }

  def runTwoFragments(cigar: String,
                      md: String,
                      readString: String,
                      strand: Strand,
                      includeInsertions: Boolean = false): Seq[Seq[ErrorRateByReadPositionMetric]] = {
    Seq(None, Some(referenceFasta)).map { case maybeRef =>
      val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
      builder.addFrag(start = 100, strand=strand).foreach { rec => rec.setAttribute("MD", "100"); rec.setReadString("A" * 100) }
      builder.addFrag(start = 100, cigar=cigar, strand=strand).foreach { rec => rec.setAttribute("MD", md); rec.setReadString(readString) }
      val input = builder.toTempFile()
      val (out, pre) = outputAndPrefix
      new ErrorRateByReadPosition(input = input, output = Some(pre), ref=maybeRef, includeInsertions=includeInsertions).execute()
      val metrics = Metric.read[ErrorRateByReadPositionMetric](out)
      metrics.length shouldBe 100
      metrics
    }
  }

  it should "compute the error rate for fragment reads with no mismatches" in {
    val metrics = runTwoFragments(cigar="100M", md="100", readString="A"*100, strand=Plus) ++
      runTwoFragments(cigar="100M", md="100", readString="A"*100, strand=Minus)
    metrics.flatten.foreach { metric =>
      metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
    }
  }

  it should "compute the error rate for fragment reads with a mismatch at the end" in {
    val metrics = runTwoFragments(cigar="100M", md="99A", readString=("A"*99)+"C", strand=Plus) ++
      runTwoFragments(cigar="100M", md="A99", readString="C"+("A"*99), strand=Minus)
    metrics.flatten.foreach { metric =>
      if (metric.position == 100) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "compute the error rate for fragment reads with a mismatch at the start" in {
    val metrics = runTwoFragments(cigar="100M", md="A99", readString="C"+("A"*99), strand=Plus) ++
      runTwoFragments(cigar="100M", md="99A", readString=("A"*99)+"C", strand=Minus)
    metrics.flatten.foreach { metric =>
      if (metric.position == 1) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "compute the error rate for fragment reads with a mismatch in the middle" in {
    val metrics = runTwoFragments(cigar="100M", md="50A49", readString=("A"*50)+"C"+("A"*49), strand=Plus) ++
      runTwoFragments(cigar="100M", md="49A50", readString=("A"*49)+"C"+("A"*50), strand=Minus)
    metrics.flatten.foreach { metric =>
      if (metric.position == 51) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "compute the error rate for fragment reads with all mismatches" in {
    val metrics = runTwoFragments(cigar="100M", md="A"*100, readString="C"*100, strand=Plus) ++
      runTwoFragments(cigar="100M", md="A"*100, readString="C"*100, strand=Minus)
    metrics.flatten.foreach { metric =>
       metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
    }
  }

  it should "compute the error rate for fragment reads with a leading insertion" in {
    val metrics = runTwoFragments(cigar="1I99M", md="99", readString="A"*100, strand=Plus, includeInsertions=true)
      runTwoFragments(cigar="99M1I", md="99", readString="A"*100, strand=Minus, includeInsertions=true)
    metrics.flatten.foreach { metric =>
      if (metric.position == 1) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "compute the error rate for fragment reads with a trailing insertion" in {
    val metrics = runTwoFragments(cigar="99M1I", md="99", readString="A"*100, strand=Plus, includeInsertions=true) ++
      runTwoFragments(cigar="1I99M", md="99", readString="A"*100, strand=Minus, includeInsertions=true)
    metrics.flatten.foreach { metric =>
      if (metric.position == 100) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "not count insertions as errors if --include-insertions is false" in {
    val metrics = runTwoFragments(cigar="99M1I", md="99", readString="A"*100, strand=Plus, includeInsertions=false) ++
      runTwoFragments(cigar="1I99M", md="99", readString="A"*100, strand=Minus, includeInsertions=false)
    metrics.flatten.foreach { metric =>
      val numReads = if (metric.position == 100) 1 else 2
      metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, numReads)
    }
  }

  it should "ignore soft-clipped bases at the start of the read" in {
    val metrics = runTwoFragments(cigar="10S90M", md="90", readString="A"*100, strand=Plus) ++
      runTwoFragments(cigar="90M10S", md="90", readString="A"*100, strand=Minus)
    metrics.flatten.foreach { metric =>
      if (metric.position <= 10) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 1l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "ignore soft-clipped bases at the end of the read" in {
    val metrics = runTwoFragments(cigar="90M10S", md="90", readString="A"*100, strand=Plus) ++
      runTwoFragments(cigar="10S90M", md="90", readString="A"*100, strand=Minus)
    metrics.flatten.foreach { metric =>
      if (90 < metric.position) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 1l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "ignore soft-clipped bases and not count insertions at the start if --include-insertions is false" in {
    val metrics = runTwoFragments(cigar="10S10I80M", md="80", readString="A"*100, strand=Plus) ++
      runTwoFragments(cigar="80M10I10S", md="80", readString="A"*100, strand=Minus)
     metrics.flatten.foreach { metric =>
      if (metric.position <= 20) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 1l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }

  it should "ignore soft-clipped bases and count insertions at the start if --include-insertions is true" in {
    val metrics = runTwoFragments(cigar="10S10I80M", md="80", readString="A"*100, strand=Plus, includeInsertions=true) ++
      runTwoFragments(cigar="80M10I10S", md="80", readString="A"*100, strand=Minus, includeInsertions=true)
    metrics.flatten.foreach { metric =>
      if (metric.position <= 10) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 1l)
      }
      else if (metric.position <= 20) {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.5d, 2l)
      }
      else {
        metric shouldBe new ErrorRateByReadPositionMetric(0, metric.position, 0.0d, 2l)
      }
    }
  }
}
