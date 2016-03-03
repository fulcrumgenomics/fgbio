/*
 * The MIT License
 *
 * Copyright (c) $year Fulcrum Genomics
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

package com.fulcrumgenomics.fastq

import java.nio.file.{Files, Path, Paths}

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.miseq.SampleSheet
import com.fulcrumgenomics.util.{Metric, ReadStructure, SampleBarcodeMetric}
import dagr.commons.CommonsDef.{FilePath, PathToFastq}
import dagr.commons.io.PathUtil
import dagr.commons.util.{LogLevel, Logger}
import org.scalatest.BeforeAndAfterAll

class DemuxFastqsTest extends UnitSpec with BeforeAndAfterAll {

  private val testDir = Paths.get("src/test/resources/com/fulcrumgenomics/util/DemuxFastqsTest/")
  private val testDirPairedDualIndexed = testDir.resolve("paired")
  private val testDirFragmentInlineIndexed = testDir.resolve("fragment_inline_indexed")

  private val logLevel = Logger.level

  override def beforeAll(): Unit = {
    htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)
    Logger.level = LogLevel.Error
  }

  override def afterAll(): Unit = {
    Logger.level = logLevel
    htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.INFO)
  }

  private def testPairedDemuxFastqs(fastq1: PathToFastq,
                                               fastq2: PathToFastq,
                                               fastq7: PathToFastq,
                                               fastq5: PathToFastq,
                                               rs1: String,
                                               rs2: String = "",
                                               rs7: String = "",
                                               rs5: String = "",
                                               sampleSheet: FilePath,
                                               maxMismatches: Int,
                                               barcodes: Array[String],
                                               numReadsPerBarcode: Array[Int],
                                               multiUmiTags: Seq[String] = Seq.empty): Unit = {

    def nonEmptyOrNone(str: String): Option[String] = if (str.isEmpty) None else Some(str)
    testDemuxFastqs(
      fastq1 = fastq1,
      fastq2 = Some(fastq2),
      fastq7 = Some(fastq7),
      fastq5 = Some(fastq5),
      rs1    = rs1,
      rs2    = nonEmptyOrNone(rs2),
      rs7    = nonEmptyOrNone(rs7),
      rs5    = nonEmptyOrNone(rs5),
      sampleSheet = sampleSheet,
      maxMismatches = maxMismatches,
      barcodes = barcodes,
      numReadsPerBarcode = numReadsPerBarcode,
      multiUmiTags = multiUmiTags
    )
  }

  private def testDemuxFastqs(fastq1: PathToFastq,
                              fastq2: Option[PathToFastq] = None,
                              fastq7: Option[PathToFastq] = None,
                              fastq5: Option[PathToFastq] = None,
                              rs1: String,
                              rs2: Option[String] = None,
                              rs7: Option[String] = None,
                              rs5: Option[String] = None,
                              sampleSheet: FilePath,
                              maxMismatches: Int,
                              barcodes: Array[String],
                              numReadsPerBarcode: Array[Int],
                              multiUmiTags: Seq[String] = Seq.empty,
                              matchReverseComplement: Boolean = false): Unit = {

    val output = Files.createTempDirectory("DemuxFastqs")
    val metricsFile: FilePath = output.resolve("DemuxFastqsTest.metrics.txt")

    val program = new DemuxFastqs(
      fastq1=List(fastq1),
      fastq2=fastq2.map(List(_)).getOrElse(List.empty),
      fastq7=fastq7.map(List(_)).getOrElse(List.empty),
      fastq5=fastq5.map(List(_)).getOrElse(List.empty),
      rs1=rs1,
      rs2=rs2,
      rs7=rs7,
      rs5=rs5,
      sampleSheet=sampleSheet,
      output=output,
      metrics=metricsFile,
      maxMismatches=maxMismatches,
      minMismatchDelta=1, // for testing
      multiUmiTags=multiUmiTags,
      matchReverseComplement=matchReverseComplement
    )

    program.execute()

    val readStructures = Seq(Some(rs1), rs2, rs7, rs5).flatten.map(ReadStructure(_))
    val samples        = DemuxFastqs.sampleSheet(sampleSheet).sampleBarcodes(readStructures:_*).toList

    // check that all samples are represented in the metrics file
    val metrics = Metric.read[SampleBarcodeMetric](metricsFile)
    metrics.size shouldBe barcodes.length
    for (i <- metrics.indices) {
      val metric = metrics(i)
      val barcode = barcodes(i)
      val numReads = numReadsPerBarcode(i)
      metric.barcode shouldBe barcode
      val sampleName = if (i == (barcodes.length-1)) DemuxFastqs.UnmatchedSampleId else s"Sample_Name_${i+1}"
      metric.barcode_name shouldBe sampleName
      metric.library_name shouldBe sampleName
      metric.reads shouldBe numReads

      // verify: (1) the templates are the correct length in the output
      //         (2) the sample barcodes are of the correct length
      //         (3) the molecular barcode bases are the correct length
      // NB: we must also count the delimiters for sample and molecular barcodes
      val outputBam                = if (i < samples.length) DemuxFastqs.sampleOutputBam(output, samples(i)) else PathUtil.pathTo(output.toString, DemuxFastqs.UnmatchedSampleId + ".bam")
      val numTemplateBasesReadOne  = ReadStructure(rs1).template.map(_.length).sum
      val numTemplateBasesReadTwo  = rs2.map(ReadStructure(_).template.map(_.length).sum).getOrElse(-1)
      val numSampleBarcodeBases    = readStructures.map(_.sampleBarcode.map(_.length).sum).sum + readStructures.count(_.sampleBarcode.nonEmpty) - 1
      val numMolecularBarcodeBases = readStructures.map(_.molecularBarcode.map(_.length).sum).sum + readStructures.count(_.molecularBarcode.nonEmpty) - 1
      readBamRecs(outputBam).foreach { rec =>
        val numTemplateBases = if (!rec.getReadPairedFlag || rec.getFirstOfPairFlag) numTemplateBasesReadOne else numTemplateBasesReadTwo
        rec.getReadLength shouldBe numTemplateBases
        rec.getReadString.length shouldBe numTemplateBases
        rec.getBaseQualityString.length shouldBe numTemplateBases
        rec.getStringAttribute("BC") match {
          case null => numSampleBarcodeBases shouldBe -1
          case mb   => mb.length shouldBe numSampleBarcodeBases
        }
        rec.getStringAttribute("RX") match {
          case null => numMolecularBarcodeBases shouldBe -1
          case mb   => mb.length shouldBe numMolecularBarcodeBases
        }
      }
    }
  }

  private def pairedDualIndexed(index: Int): Path = {
    if      (index == 0) testDirPairedDualIndexed.resolve("whole_run_S0_L001_R1_001.fastq")
    else if (index == 1) testDirPairedDualIndexed.resolve("whole_run_S0_L001_R2_001.fastq")
    else if (index == 2) testDirPairedDualIndexed.resolve("whole_run_S0_L001_I1_001.fastq")
    else if (index == 3) testDirPairedDualIndexed.resolve("whole_run_S0_L001_I2_001.fastq")
    else throw new IllegalStateException("Index out of range: " + index)
  }

  "DemixFastqs" should "demux with molecular id on one index read, sample barcode on the other" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7B", "7M", sampleSheet, 1,
      Array[String]("GATTACA", "GATTACC", "NNNNNNN"), Array[Int](1, 1, 1))
  }

  it should "demux both molecular id and sample barcode present on both index reads" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.dual_indexed.1B.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "6M1B", "1B6M", sampleSheet, 1,
      Array[String]("A-A", "C-C", "N-N"), Array[Int](1, 1, 1))
  }

  it should "demux molecular id on one index read, sample barcode on the other, tests matching with at most one mismatch (third read is two mismatches away from sample two, so no match)" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index2.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7M", "7B", sampleSheet, 1,
      Array[String]("ACATTAG", "CGGGGGG", "NNNNNNN"), Array[Int](2, 0, 1))
  }

  it should "demux molecular id on one index read, sample barcode on the other, tests matching with two mismatches (now the third read matches sample two)" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index2.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7M", "7B", sampleSheet, 2,
      Array[String]("ACATTAG", "CGGGGGG", "NNNNNNN"), Array[Int](2, 1, 0))
  }

  it should "demux with molecular id on one index read, sample barcode on the other, with a single molecular id SAM tag" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7B", "7M", sampleSheet, 1,
      Array[String]("GATTACA", "GATTACC", "NNNNNNN"), Array[Int](1, 1, 1), Seq("MI"))
  }

  it should "demux with molecular id on both index reads, with a single molecular id SAM tag" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.3B.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "4M3B", "7M", sampleSheet, 1,
      Array[String]("ACA", "ACC", "NNN"), Array[Int](1, 1, 1), Seq("MI"))
  }

  it should "demux with molecular id on both index reads, with a two molecular id SAM tags" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.3B.csv")
    testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "4M3B", "7M", sampleSheet, 1,
      Array[String]("ACA", "ACC", "NNN"), Array[Int](1, 1, 1), Seq("M1", "M2"))
  }

  it should "demux inline sample barcodes on both ends of the read for fragment reads" in {
    val sampleSheet = testDirFragmentInlineIndexed.resolve("SampleSheet.dual_indexed.csv")
    val fastq1 = testDirFragmentInlineIndexed.resolve("whole_run_S0_L001_R1_001.fastq")
    testDemuxFastqs(fastq1=fastq1, fastq2=None, fastq7=None, fastq5=None, "7B137T7B", None, None, None, sampleSheet, 1,
      Array[String]("GATTACAGATTACA", "CTAATGTCTAATGT", "NNNNNNNNNNNNNN"), Array[Int](1, 1, 1))
  }

  it should "demux inline a sample barcode on the start of the read for fragment" in {
    val sampleSheet = testDirFragmentInlineIndexed.resolve("SampleSheet.single_indexed.csv")
    val fastq1 = testDirFragmentInlineIndexed.resolve("whole_run_S0_L001_R1_001.fastq")
    testDemuxFastqs(fastq1=fastq1, fastq2=None, fastq7=None, fastq5=None, "7B144T", None, None, None, sampleSheet, 1,
      Array[String]("GATTACA", "CTAATGT", "NNNNNNN"), Array[Int](1, 1, 1))
  }

  it should "demux inline a sample barcode on the start of the read for fragment and consider the reverse complement" in {
    val sampleSheet = testDirFragmentInlineIndexed.resolve("SampleSheet.single_indexed.csv")
    val fastq1 = testDirFragmentInlineIndexed.resolve("whole_run_S0_L001_R1_001.rc.fastq")
    testDemuxFastqs(fastq1=fastq1, fastq2=None, fastq7=None, fastq5=None, "7B144T", None, None, None, sampleSheet, 1,
      Array[String]("GATTACA", "CTAATGT", "NNNNNNN"), Array[Int](1, 1, 1), matchReverseComplement=true)
  }

  it should "throw an exception when molecular barcode SAM tags are given but no molecular barcodes are found" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.csv")
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7B", "7S", sampleSheet, 1,
      Array.empty[String], Array.empty[Int], Seq("M1", "M2"))
  }

  it should "throw an exception when not enough molecular barcode SAM tags are given" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.csv")
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "1M1S1M1S3M", "7S", sampleSheet, 1,
      Array.empty[String], Array.empty[Int], Seq("M1", "M2"))
  }

  it should "throw an exception with a sample barcode in the read structure but not the sample sheet" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.csv")
    // index 2 has a sample barcode in the read structure, but the sample sheet does not.
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7B", "7B", sampleSheet, 1,
      Array.empty[String], Array.empty[Int])
  }

  it should "throw an exception with a sample barcode the sample sheet but not the read structure" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.csv")
    // index 1 has a sample barcode in the sample sheet, but the read structure does not.
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7M", "7M", sampleSheet, 1,
      Array.empty[String], Array.empty[Int])
  }

  it should "throw an exception when not enough read structures are given" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.index1.csv")
    // index 2 has a FASTQ but no read structure
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "7M", "", sampleSheet, 1,
      Array.empty[String], Array.empty[Int])
  }

  it should "throw an exception when the sample barcode # of bases in the sample sheet does not match the # of sample barcodes in the read structures" in {
    val sampleSheet = testDirPairedDualIndexed.resolve("SampleSheet.dual_indexed.1B.csv")
    // second index read has too many sample barcode bases in the index read
    an[Exception] should be thrownBy testPairedDemuxFastqs(pairedDualIndexed(0), pairedDualIndexed(1), pairedDualIndexed(2), pairedDualIndexed(3), "151T", "151T", "6M1B", "2B5M", sampleSheet, 1,
      Array.empty[String], Array.empty[Int])
  }

  "DemuxFastqs.countMismatches" should "find no mismatches" in {
    DemuxFastqs.countMismatches("GATTACA".getBytes, "GATTACA".getBytes) shouldBe 0
  }

  it should "find two mismatches" in {
    DemuxFastqs.countMismatches("GATTACA".getBytes, "GACCACA".getBytes) shouldBe 2
  }

  it should "not count no calls" in {
    DemuxFastqs.countMismatches("GATTACA".getBytes, "GANNACA".getBytes) shouldBe 0
  }

  it should "find compare two sequences that have all mismatches" in {
    DemuxFastqs.countMismatches("GATTACA".getBytes, "CTAATGT".getBytes) shouldBe 7
  }
}
