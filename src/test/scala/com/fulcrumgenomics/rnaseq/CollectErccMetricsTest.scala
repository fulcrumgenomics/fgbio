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
 *
 */

package com.fulcrumgenomics.rnaseq

import java.nio.file.Files

import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathPrefix}
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Metric, Rscript}
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}
import org.scalatest.OptionValues

import scala.io.Source

class CollectErccMetricsTest extends UnitSpec with OptionValues {

  private val dict = {
    val sd = new SAMSequenceDictionary()
    sd.addSequence(new SAMSequenceRecord("chr1", 10000))
    sd.addSequence(new SAMSequenceRecord("ERCC-00001", 1000))
    sd.addSequence(new SAMSequenceRecord("ERCC-00002", 1000))
    sd.addSequence(new SAMSequenceRecord("ERCC-00003", 1000))
    sd.addSequence(new SAMSequenceRecord("ERCC-00004", 1000))
    sd.addSequence(new SAMSequenceRecord("ERCC-00005", 1000))
    sd
  }

  private val builder = new SamBuilder(sd=Some(dict))
  builder.addFrag(contig=0, start=1) // non-ERCC frag
  builder.addFrag(contig=1, start=1) // ERCC frag
  builder.addPair(contig=0, start1=1, start2=2) // non-ERCC pair
  builder.addPair(contig=1, start1=1, start2=2).foreach { r => r.secondary = true } // ERCC pair, secondary
  builder.addPair(contig=1, start1=1, start2=2).foreach { r => r.supplementary = true } // ERCC pair, supplementary
  builder.addPair(contig=1, start1=1, start2=2, mapq1 = 1, mapq2 = 1) // ERCC pair, low mapq
  builder.addPair(contig=0, start1=1, start2=2).foreach { r => if (r.secondOfPair) r.refIndex = 1 } // second end is ERCC
  builder.addPair(contig=1, start1=1, start2=2) // ERCC-00001 pair
  Seq.range(0, 2).foreach { _ => builder.addPair(contig=2, start1=1, start2=2) } // ERCC-00002 pair
  Seq.range(0, 4).foreach { _ => builder.addPair(contig=3, start1=1, start2=2) } // ERCC-00002 pair
  Seq.range(0, 8).foreach { _ => builder.addPair(contig=4, start1=1, start2=2) } // ERCC-00002 pair
  Seq.range(0, 16).foreach { _ => builder.addPair(contig=5, start1=1, start2=2) } // ERCC-00002 pair

  private val bam = builder.toTempFile()

  private val metadata = {
    val path = makeTempFile("CollectErccMetricsTest.", ".tab")
    val lines = Seq(
      Seq("Id", "Concentration"),
      Seq("ERCC-00001", 1),
      Seq("ERCC-00002", 2),
      Seq("ERCC-00003", 4),
      Seq("ERCC-00004", 8),
      Seq("ERCC-00005", 16)
    ).map(_.mkString("\t"))
    Io.writeLines(path, lines)
    path
  }

  private def toOutput: FilePath = {
    val output = Files.createTempDirectory("CollectErccMetricsTest.")
    output.toFile.deleteOnExit()
    output.resolve("output")
  }

  object Outputs {
    def apply(prefix: PathPrefix): Outputs = {
      def f(ext: String): FilePath = PathUtil.pathTo(prefix + ext)
      Outputs(
        summaryPath  = f(".ercc_summary_metrics.txt"),
        detailedPath = f(".ercc_detailed_metrics.txt"),
        plotPath     = f(".ercc_plot.pdf")
      )
    }
  }

  case class Outputs(summaryPath: FilePath, detailedPath: FilePath, plotPath: FilePath)

  "CollectErccMetrics" should "fail if the metadata file did not contain any ERCC transcripts" in {
    val metadata = {
      val path = makeTempFile("CollectErccMetricsTest.", ".tab")
      Io.writeLines(path, Seq("Id\tConcentration"))
      path
    }

    val output = makeTempFile("CollectErccMetricsTest.", ".tab")

    an[IllegalArgumentException] should be thrownBy new CollectErccMetrics(input=bam, customMixture=Some(metadata), output=output).execute()
  }

  it should "fail if not all the ERCC IDs in the metadata file are found in the SAM/BAM header" in {
    val metadata = {
      val path = makeTempFile("CollectErccMetricsTest.", ".tab")
      val lines = Seq(
        Seq("Id", "Concentration"),
        Seq("ERCC-00001", 1),
        Seq("ERCC-00002", 2),
        Seq("ERCC-00003", 4),
        Seq("ERCC-00004", 8),
        Seq("ERCC-00005", 16),
        Seq("ERCC-00006", 32) // an extra ERCC
      ).map(_.mkString("\t"))
      Io.writeLines(path, lines)
      path
    }

    an[IllegalArgumentException] should be thrownBy new CollectErccMetrics(input=bam, customMixture=Some(metadata), output=toOutput).execute()
  }

  it should "filter secondary alignments" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=1, start1=1, start2=2).foreach { r => r.secondary = true } // ERCC pair, secondary

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 0
  }

  it should "filter supplementary alignments" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=1, start1=1, start2=2).foreach { r => r.supplementary = true } // ERCC pair, supplementary

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 0
  }

  it should "require minimum mapq for ERCC records" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=1, start1=1, start2=2, mapq1 = 1, mapq2 = 1) // ERCC pair, low mapq

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 2
    metrics.head.ercc_reads shouldBe 0
  }

  it should "require both ends to map to the same ERCC for paired reads/fragments" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=0, start1=1, start2=2).foreach { r => if (r.secondOfPair) r.refIndex = 1 } // second end is ERCC

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 2
    metrics.head.ercc_reads shouldBe 1
    metrics.head.ercc_templates shouldBe 0
  }

  it should "require at least two observations to compute correlation coefficients and other metrics" in {

    // one transcript
    {
      val builder = new SamBuilder(sd=Some(dict))
      builder.addPair(contig=1, start1=1, start2=2) // ERCC-00001 pair

      val output = toOutput
      new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output, minTranscriptCount=1).execute()
      val outputs = Outputs(output)
      val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
      metrics.head.total_reads shouldBe 2
      metrics.head.ercc_reads shouldBe 2
      metrics.head.ercc_templates shouldBe 1
      metrics.head.pearsons_correlation shouldBe 'empty
      metrics.head.slope shouldBe 'empty
      Files.exists(outputs.plotPath) shouldBe false
    }

    // two different transcripts with counts
    {
      val builder = new SamBuilder(sd=Some(dict))
      builder.addPair(contig=1, start1=1, start2=2) // ERCC-00001 pair
      builder.addPair(contig=2, start1=1, start2=2) // ERCC-00002 pair
      builder.addPair(contig=2, start1=1, start2=2) // ERCC-00002 pair

      val output = toOutput
      new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output, minTranscriptCount=1).execute()
      val outputs = Outputs(output)
      val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
      metrics.head.total_reads shouldBe 6
      metrics.head.ercc_reads shouldBe 6
      metrics.head.ercc_templates shouldBe 3
      metrics.head.pearsons_correlation.value shouldBe 1
      metrics.head.slope.value shouldBe 1
      Files.exists(outputs.plotPath) shouldBe Rscript.Available
    }
  }

  it should "require a minimum transcript count" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=1, start1=1, start2=2) // ERCC-00001 pair
    builder.addPair(contig=2, start1=1, start2=2) // ERCC-00002 pair
    builder.addPair(contig=2, start1=1, start2=2) // ERCC-00002 pair

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output, minTranscriptCount=2).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 6
    metrics.head.ercc_reads shouldBe 6
    metrics.head.ercc_templates shouldBe 3
    metrics.head.total_transcripts shouldBe 2
    metrics.head.passing_filter_transcripts shouldBe 1
  }

  it should "support single end reads" in {
    val builder = new SamBuilder(sd=Some(dict))
    builder.addFrag(contig=0, start=1) // non-ERCC frag
    builder.addFrag(contig=1, start=1) // ERCC-00001 frag
    builder.addFrag(contig=1, start=1) // ERCC-00001 frag
    builder.addFrag(contig=2, start=1) // ERCC-00002 frag
    builder.addFrag(contig=2, start=1) // ERCC-00002 frag
    builder.addFrag(contig=2, start=1) // ERCC-00002 frag

    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output, minTranscriptCount=2).execute()
    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    metrics.head.total_reads shouldBe 6
    metrics.head.ercc_reads shouldBe 5
    metrics.head.ercc_templates shouldBe 5
    metrics.head.total_transcripts shouldBe 2
    metrics.head.passing_filter_transcripts shouldBe 2
  }

  it should "support using standard ERCC mixtures" in {
    val dict = {
      // Read in the standard file
      val stream = getClass.getResourceAsStream(s"${CollectErccMetrics.StandardMixtureMetadataPath}")
      val lines  = Source.fromInputStream(stream).withClose(() => stream.close()).getLines.filterNot(_.startsWith("#")).toSeq
      val names  = new DelimitedDataParser(lines, '\t').map { row => row[String](1)  }.toSeq

      // Create a dummy dictionary
      val sd = new SAMSequenceDictionary()
      names.foreach { name => sd.addSequence(new SAMSequenceRecord(name, 1000)) }
      sd.addSequence(new SAMSequenceRecord("chr1", 10000))
      sd
    }

    val builder = new SamBuilder(sd=Some(dict))
    builder.addPair(contig=63, start1=1, start2=2) // ERCC-00120
    builder.addPair(contig=62, start1=1, start2=2) // ERCC-00058
    builder.addPair(contig=62, start1=1, start2=2) // ERCC-00058

    // Mix #1
    {
      val output = toOutput
      new CollectErccMetrics(input=builder.toTempFile(), mixtureName=Some(ErccMixture.Mix1), output=output, minTranscriptCount=2).execute()
      val outputs = Outputs(output)
      val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
      metrics.head.total_reads shouldBe 6
      metrics.head.ercc_reads shouldBe 6
      metrics.head.ercc_templates shouldBe 3
      metrics.head.total_transcripts shouldBe 2
      metrics.head.passing_filter_transcripts shouldBe 1
    }

    // Mix #2
    {
      val output = toOutput
      new CollectErccMetrics(input=builder.toTempFile(), mixtureName=Some(ErccMixture.Mix2), output=output, minTranscriptCount=2).execute()
      val outputs = Outputs(output)
      val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
      metrics.head.total_reads shouldBe 6
      metrics.head.ercc_reads shouldBe 6
      metrics.head.ercc_templates shouldBe 3
      metrics.head.total_transcripts shouldBe 2
      metrics.head.passing_filter_transcripts shouldBe 1
    }
  }

  it should "run end-to-end" in {
    val output = toOutput
    new CollectErccMetrics(input=builder.toTempFile(), customMixture=Some(metadata), output=output).execute()

    val outputs = Outputs(output)
    val metrics = Metric.read[ErccSummaryMetrics](outputs.summaryPath)
    val metric = metrics.head
    metric.total_reads shouldBe 70 // 2 frags, and 34 pairs (= 2 + 34*2 = 70 reads)
    metric.ercc_reads shouldBe 64 // 2 frags, 2 two ends of pairs, are not ERCC ( = 70 - 2 - 2 = 66)
    metric.fraction_ercc_reads shouldBe 0.914286
    metric.ercc_templates shouldBe 32 // one frag not mapping to an ERCC, one pair did not have both ends map to an ERCC, one pair with low mapq, and one pair with ends mapping to different ERCC transcripts ( = 36 - 4 = 32)
    metric.total_transcripts shouldBe 5
    metric.passing_filter_transcripts shouldBe 3 // ERCC-00001 and ERCC-00002 do not meet the threshold
    metric.pearsons_correlation.value shouldBe 1
    metric.spearmans_correlation.value shouldBe 1
    metric.slope.value shouldBe 1
    metric.intercept.value shouldBe 0
    metric.r_squared.value shouldBe 1
  }
}
