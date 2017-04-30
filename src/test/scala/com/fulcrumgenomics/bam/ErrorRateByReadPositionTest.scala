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

import java.nio.file.{Files, Path}

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.fasta.ReferenceSetBuilder
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec, VariantContextSetBuilder}
import com.fulcrumgenomics.util.{Metric, Rscript}
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

/**
  * Tests for ErrorRateByReadPosition.
  */
class ErrorRateByReadPositionTest extends UnitSpec {
  /////////////////////////////////////////////////////////////////////////////
  // Make a reference and a set of variants for use in the tests
  /////////////////////////////////////////////////////////////////////////////
  private val ref = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr0").add("RYSWKMBDHV", 50).add("N", 500) // Sooo ambiguous
    builder.add("chr1").add("A", 1000)
    builder.add("chr2").add("C", 1000)
    builder.add("chr3").add("G", 1000)
    builder.add("chr4").add("T", 1000)
    builder.toTempFile()
  }

  private val dict = SAMSequenceDictionaryExtractor.extractDictionary(ref.toFile)

  private val vcf = {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(1, 500, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"))
    builder.addVariant(2, 500, variantAlleles=List("C", "T"), genotypeAlleles=List("C", "T"))
    builder.addVariant(3, 500, variantAlleles=List("G", "A"), genotypeAlleles=List("G", "A"))
    builder.addVariant(4, 500, variantAlleles=List("T", "C"), genotypeAlleles=List("T", "C"))
    builder.toTempFile()
  }

  private def outputAndPrefix: (Path, Path) = {
    val out = makeTempFile("ErrorRateByReadPositionTest.", ErrorRateByReadPositionMetric.FileExtension)
    val pre = PathUtil.pathTo(out.toString.replace(ErrorRateByReadPositionMetric.FileExtension, ""))
    (out, pre)
  }

  private def newSamBuilder = {
    val builder = new SamBuilder(readLength=20, sort=Some(SamOrder.Coordinate))
    builder.header.setSequenceDictionary(dict)
    builder
  }

  "ErrorRateByReadPosition" should "work on an empty BAM" in {
    val builder = newSamBuilder
    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref).computeMetrics
    metrics.size shouldBe 0
  }

  it should "work on a file with all unmapped reads" in {
    val builder = newSamBuilder
    Range.inclusive(1, 20).foreach { i => builder.addFrag(unmapped=true) }
    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref).computeMetrics
    metrics.size shouldBe 0
  }

  it should "compute the error rate for some simple paired end reads" in {
    val builder = newSamBuilder
    1 to 9 foreach {i => builder.addPair(name="$i", contig=1, start1=100, start2=200, bases1="A"*20, bases2="A"*20) }
    builder.addPair("err", contig=1, start1=300, start2=400).foreach { rec =>
      if (rec.firstOfPair) rec.bases = "AAAAAAAAAAAAAAAAACGT"
      else rec.bases = "AAAAAAAAAAAAAAAAACGT"
    }

    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref).computeMetrics
    metrics.size shouldBe 40
    metrics.foreach { m =>
      m.bases_total shouldBe 10
      m.error_rate shouldBe m.errors / m.bases_total.toDouble
      if (m.read_number == 1 && m.position >= 18) {
        m.error_rate shouldBe 0.1
        m.position match {
          case 18 => m.a_to_c_error_rate shouldBe 0.1
          case 19 => m.a_to_g_error_rate shouldBe 0.1
          case 20 => m.a_to_t_error_rate shouldBe 0.1
        }
      }
      else if (m.read_number == 2 && m.position <= 3) {
        m.position match {
          case 1 => m.a_to_t_error_rate shouldBe 0.1
          case 2 => m.a_to_g_error_rate shouldBe 0.1
          case 3 => m.a_to_c_error_rate shouldBe 0.1
        }
      }
      else {
        m.error_rate shouldBe 0
      }
    }
  }

  it should "not count as errors either Ns in reads or ambiguous bases in the reference" in {
    val builder = newSamBuilder
    builder.addPair(contig=0, start1=100, start2=100, bases1="A"*20  , bases2="A"*20   )
    builder.addPair(contig=1, start1=100, start2=200, bases1="AANA"*5, bases2="AANA"*5 )
    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref).computeMetrics
    metrics.size shouldBe 40
    metrics.map(_.bases_total).sum shouldBe 30 // nothing on contig 0, and 3/4 on contig 1
    metrics.forall(_.error_rate == 0.0) shouldBe true
  }

  it should "restrict to a set of intervals if provided" in {
    val builder = newSamBuilder
    builder.addPair(contig=1, start1=100, start2=150, bases1="AAAA"*5, bases2="AAAA"*5)
    builder.addPair(contig=1, start1=200, start2=250, bases1="AAAA"*5, bases2="AAAA"*5)
    builder.addPair(contig=1, start1=300, start2=350, bases1="ACGA"*5, bases2="ACGA"*5)
    builder.addPair(contig=1, start1=400, start2=450, bases1="AAAA"*5, bases2="AAAA"*5)

    val intervals = new IntervalList(dict)
    intervals.add(new Interval("chr1", 100, 275))
    intervals.add(new Interval("chr1", 400, 500))
    val intervalPath = makeTempFile("regions.", ".interval_list")
    intervals.write(intervalPath.toFile)

    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref, intervals=Some(intervalPath)).computeMetrics
    metrics.forall(_.error_rate == 0.0) shouldBe true
    metrics.map(_.bases_total).sum shouldBe (3 * 2 * 20)
  }

  it should "not count errors that occur at sites with variants" in {
    val builder = newSamBuilder
    builder.addFrag(contig=1, start=490, bases="GAAAAAAAAAGAAAAAAAAA")
    builder.addFrag(contig=2, start=490, bases="TCCCCCCCCCTCCCCCCCCC")
    builder.addFrag(contig=3, start=490, bases="AGGGGGGGGGAGGGGGGGGG")
    builder.addFrag(contig=4, start=490, bases="CTTTTTTTTTTTTTTTTTTT")

    val (out, pre) = outputAndPrefix
    val metrics = new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref, variants=Some(vcf)).computeMetrics
    metrics.find(_.position == 1).get.error_rate shouldBe 1.0
    metrics.find(_.position == 11).get.error_rate shouldBe 0.0
    metrics.filter(m => m.position != 1 && m.position != 11).foreach { m =>
      m.bases_total shouldBe 4
      m.error_rate  shouldBe 0.0
    }
  }

  it should "run end to end and maybe generate a PDF" in {
    val builder = newSamBuilder
    1 to 99 foreach {i => builder.addPair(contig=1, start1=i, start2=i+50, bases1="A"*20, bases2="A"*20) }
    builder.addPair(contig=1, start1=100, start2=200, bases1="AAAAAAAAATTAAAAAAAAA", bases2="AAAAAAAAATTAAAAAAAAA")

    val (out, pre) = outputAndPrefix
    new ErrorRateByReadPosition(input=builder.toTempFile(), output=Some(pre), ref=ref).execute()
    val metrics = Metric.read[ErrorRateByReadPositionMetric](out)
    metrics.size shouldBe 40
    metrics.map(_.bases_total).sum shouldBe 4000
    metrics.foreach { m =>
      if (m.position == 10 || m.position == 11) m.error_rate shouldBe 0.01
      else m.error_rate shouldBe 0.0
    }

    if (Rscript.Available) {
      val plot = PathUtil.pathTo(pre + ErrorRateByReadPositionMetric.PlotExtension)
      Files.exists(plot) shouldBe true
    }
  }

  it should "handle interval lists that are unsorted and/or contain duplicate entries" in {
    val builder = newSamBuilder
    1 to 99 foreach {i => builder.addPair(contig=1, start1=i, start2=i+50, bases1="A"*20, bases2="A"*20) }
    builder.addPair(contig=1, start1=100, start2=200, bases1="AAAAAAAAATTAAAAAAAAA", bases2="AAAAAAAAATTAAAAAAAAA")

    val dict = builder.dict
    val intervals = new IntervalList(dict)
    intervals.add(new Interval(dict.getSequence(1).getSequenceName, 200, 300))
    intervals.add(new Interval(dict.getSequence(1).getSequenceName, 100, 150))
    intervals.add(new Interval(dict.getSequence(1).getSequenceName, 100, 200))
    intervals.add(new Interval(dict.getSequence(1).getSequenceName, 1,  1000))

    val intervalsPath = makeTempFile("regions.", ".interval_list")
    intervals.write(intervalsPath.toFile)

    val (out, pre) = outputAndPrefix
    new ErrorRateByReadPosition(input=builder.toTempFile(), intervals=Some(intervalsPath), output=Some(pre), ref=ref, variants=Some(vcf)).execute()
    val metrics = Metric.read[ErrorRateByReadPositionMetric](out)
    metrics.size shouldBe 40
    metrics.map(_.bases_total).sum shouldBe 4000
  }
}
