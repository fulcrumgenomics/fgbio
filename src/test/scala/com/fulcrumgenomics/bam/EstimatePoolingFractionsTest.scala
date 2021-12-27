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

package com.fulcrumgenomics.bam

import java.nio.file.Paths
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Metric
import com.fulcrumgenomics.vcf.api
import com.fulcrumgenomics.vcf.api.{Genotype, VcfCount, VcfFieldType, VcfFormatHeader, VcfSource, VcfWriter}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.{MergingSamRecordIterator, SamFileHeaderMerger}
import org.scalatest.ParallelTestExecution

class EstimatePoolingFractionsTest extends UnitSpec with ParallelTestExecution {
  private val Samples = Seq("HG01879", "HG01112", "HG01583", "HG01500", "HG03742", "HG03052")
  private val DataDir = Paths.get("src/test/resources/com/fulcrumgenomics/bam/estimate_pooling_fractions")
  private val Bams    = Samples.map(s => DataDir.resolve(s + ".bam"))
  private val Vcf     = DataDir.resolve("variants.vcf.gz")
  private val Regions = DataDir.resolve("regions.interval_list")

  /** Merges one or more BAMs and returns the path to the merged BAM. */
  def merge(bams: Seq[PathToBam]): PathToBam = {
    val readers = bams.map(bam => SamSource(bam))

    // Mangle the library names in the header so that the merger sees duplicate RGs as different RGs.
    readers.zipWithIndex.foreach { case (reader, index) =>
        reader.header.getReadGroups.foreach(rg => rg.setLibrary(rg.getLibrary + ":" + index))
    }
    val headerMerger = new SamFileHeaderMerger(SortOrder.coordinate, readers.iterator.map(_.header).toJavaList, false)
    val iterator     = new MergingSamRecordIterator(headerMerger, readers.iterator.map(_.toSamReader).toJavaList, true)

    val output = makeTempFile("merged.", ".bam")
    val out    = SamWriter(output, headerMerger.getMergedHeader, compression = 0)
    iterator.map(_.asInstanceOf[SamRecord]).foreach { r =>
       // Add the RG ID to the read name so we don't have identical read names when merging the same BAM 2+ times
      r.name = r.readGroup.getReadGroupId + ":" + r.name
      out += r
    }
    out.close()
    readers.foreach(_.safelyClose())
    output
  }

  "EstimatePoolingFractions" should "estimate approximately 50/50 for two samples mixed 50/50" in {
    val bam = merge(Bams.take(2))
    val out = makeTempFile("pooling_metrics.", ".txt")
    new EstimatePoolingFractions(vcf=Vcf, bam=bam, output=out, samples=Samples.take(2)).execute()
    val metrics = Metric.read[PoolingFractionMetric](out)
    metrics should have size 2
    metrics.foreach(m => 0.5 should (be >= m.ci99_low and be <= m.ci99_high))
  }

  Range.inclusive(3, Samples.size-1).foreach { n =>
    it should s"accurately estimate a mixof $n samples" in {
        val bam = merge(Bams.take(n))
        val out = makeTempFile("pooling_metrics.", ".txt")
        new EstimatePoolingFractions(vcf=Vcf, bam=bam, output=out, samples=Samples.take(n)).execute()
        val metrics = Metric.read[PoolingFractionMetric](out)
        metrics should have size n
        metrics.foreach(m => (1/n.toDouble) should (be >= m.ci99_low and be <= m.ci99_high))
    }
  }

  it should "work with an interval list, and also use all samples if no samples are provided" in {
    val bam = merge(Bams)
    val out = makeTempFile("pooling_metrics.", ".txt")
    new EstimatePoolingFractions(vcf=Vcf, bam=bam, output=out, intervals=Seq(Regions)).execute()
    val metrics = Metric.read[PoolingFractionMetric](out)
    metrics should have size Samples.size
    metrics.foreach(m => (1/Samples.size.toDouble) should (be >= m.ci99_low and be <= m.ci99_high))
  }

  it should "accurately estimate unequal mixes of two samples" in {
    val samples         = Samples.take(2)
    val Seq(bam1, bam2) = Bams.take(2)
    val bam = merge(Seq(bam1, bam1, bam1, bam2))
    val out = makeTempFile("pooling_metrics.", ".txt")
    new EstimatePoolingFractions(vcf=Vcf, bam=bam, output=out, samples=samples).execute()
    val metrics = Metric.read[PoolingFractionMetric](out)
    metrics should have size 2
    metrics.foreach {m =>
      val expected = if (m.sample == samples.head) 0.75 else 0.25
      expected should (be >= m.ci99_low and be <= m.ci99_high)
    }
  }

  it should "accurately estimate a three sample mixture using the AF genotype field" in {
    val samples         = Samples.take(3)
    val Seq(s1, s2, s3) = samples
    val bams            = Bams.take(3)
    val bam             = merge(bams)

    val vcf = {
      val vcf = makeTempFile("mixture.", ".vcf.gz")
      val in  = api.VcfSource(Vcf)
      val hd  = in.header.copy(
        samples = IndexedSeq(s1, "two_sample_mixture"),
        formats = VcfFormatHeader("AF", VcfCount.OnePerAltAllele, kind=VcfFieldType.Float, description="Allele Frequency") +: in.header.formats
      )
      val out = VcfWriter(vcf, hd)

      in.filter(_.alleles.size == 2).foreach { v =>
        val gts = samples.map(v.gt)

        // Only bother with sites where all samples have called genotypes and there is variation
        if (gts.forall((_.isFullyCalled)) && gts.flatMap(_.calls).toSet.size > 1) {
          // Make a mixture of the 2nd and 3rd samples
          val (mixCalls, mixAf) = {
            val input = gts.drop(1)
            if      (input.forall(_.isHomRef)) (IndexedSeq(v.alleles.ref, v.alleles.ref), 0.0)
            else if (input.forall(_.isHomVar)) (IndexedSeq(v.alleles.alts.head, v.alleles.alts.head), 1.0)
            else {
              val calls = input.flatMap(_.calls)
              (IndexedSeq(v.alleles.ref, v.alleles.alts.head), calls.count(_ != v.alleles.ref) / calls.size.toDouble)
            }
          }

          val mixtureGt = Genotype(
            alleles = v.alleles,
            sample  = "two_sample_mixture",
            calls   = mixCalls,
            attrs   = Map("AF" -> IndexedSeq[Float](mixAf.toFloat))
          )

          out += v.copy(genotypes=Map(s1 -> gts.head, mixtureGt.sample -> mixtureGt))
        }
      }

      in.safelyClose()
      out.close()
      vcf
    }

    // Run the estimator and test the outputs
    val out = makeTempFile("pooling_metrics.", ".txt")
    new EstimatePoolingFractions(vcf=vcf, bam=bam, output=out, minGenotypeQuality = -1).execute()
    val metrics = Metric.read[PoolingFractionMetric](out)

    metrics should have size 2
    metrics.foreach {m =>
      val expected = if (m.sample == samples.head) 1/3.0 else 2/3.0
      expected should (be >= m.ci99_low and be <= m.ci99_high)
    }
  }
}
