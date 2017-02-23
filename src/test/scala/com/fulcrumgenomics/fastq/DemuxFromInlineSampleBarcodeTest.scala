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

package com.fulcrumgenomics.fastq

import java.nio.file.Files

import com.fulcrumgenomics.testing.{ErrorLogLevel, SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Io, Metric, ReadStructure, SampleBarcodeMetric}
import com.fulcrumgenomics.FgBioDef.FilePath
import com.fulcrumgenomics.fastq.DemuxFromInlineSampleBarcode.UnmatchedSampleId
import com.fulcrumgenomics.fastq.FastqDemultiplexer.DemuxRecord
import org.scalatest.OptionValues
import com.fulcrumgenomics.util.miseq.{Sample, SampleBarcode, SampleSheet}
import dagr.commons.io.PathUtil
import htsjdk.samtools.SAMUtils

import scala.collection.mutable.ListBuffer

class DemuxFromInlineSampleBarcodeTest extends UnitSpec with OptionValues with ErrorLogLevel {
  import DemuxFromInlineSampleBarcode._

  private val sampleSheetPath: FilePath = {
    val path = makeTempFile("SampleSheet", ".csv")
    val lines = """
      |[Header],,,,,,,,,
      |IEMFileVersion,4,,,,,,,,
      |Investigator Name,Joe,,,,,,,,
      |Experiment Name,EXPID,,,,,,,,
      |Date,01/01/2000,,,,,,,,
      |Workflow,GenerateFASTQ,,,,,,,,
      |Application,FASTQ Only,,,,,,,,
      |Assay,Assay Name,,,,,,,,
      |Description,The Description,,,,,,,,
      |Chemistry,Amplicon,,,,,,,,
      |,,,,,,,,,
      |[Reads],,,,,,,,,
      |151,,,,,,,,,
      |151,,,,,,,,,
      |,,,,,,,,,
      |[Settings],,,,,,,,,
      |ReverseComplement,0,,,,,,,,
      |Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
      |AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
      |,,,,,,,,,
      |[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,R1_Barcode_Bases,R2_Barcode_Bases,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,AAAAAAAA,GATTACAGA,I7_1,GATTACAACGT,I5_1,AAAAAAA,Sample_Project_1,Description_1
      |20000101-EXPID-2,Sample_Name_2,,,CCCCCCCC,GATTACAGA,I7_2,GATTACAACGT,I5_2,CCCCCCC,Sample_Project_2,Description_2
      |20000101-EXPID-3,Sample_Name_3,,,GGGGGGGG,GATTACAGA,I7_3,GATTACAACGT,I5_3,GGGGGGG,Sample_Project_3,Description_3
      |20000101-EXPID-4,Sample_Name_4,,,GGGGGGTT,GATTACAGA,I7_4,GATTACAACGT,I5_4,GGGGGTT,Sample_Project_4,Description_4
    """.stripMargin.trim
    Io.writeLines(path=path, lines=lines.split("\n"))
    path
  }

  private val sampleSheet      = CustomSampleSheet.maybeCustomSampleSheet(sampleSheet=this.sampleSheetPath, useIndexColumns=false).setSampleBarcodes()
  private val sampleHeaders    = this.sampleSheet.map(createSamFileHeader).toSeq
  private val unmatchedSample  = {
    val sample           = new Sample(this.sampleHeaders.length, UnmatchedSampleId, UnmatchedSampleId, Some(UnmatchedSampleId), None, None)
    sample.sampleBarcode = Some(new SampleBarcode(Seq("N"*17)))
    sample
  }
  private val samples          = this.sampleSheet.setSampleBarcodes().toSeq ++ Seq(this.unmatchedSample)
  private val headers           = {
    val unmatchedHeader = createSamFileHeader(this.unmatchedSample)
    this.sampleHeaders ++ Seq(unmatchedHeader)
  }

  "CustomSampleSheet.maybeCustomSampleSheet" should "use the index/index2 columns for sample barcodes when useIndexColumns=true" in {
    val samples = CustomSampleSheet.maybeCustomSampleSheet(sampleSheet=this.sampleSheetPath, useIndexColumns=true)
    samples.foreach { sample =>
      val barcode = sample.sampleBarcode.value
      barcode.concatenatedBarcode.startsWith("GATTACAACGT") shouldBe true
    }
  }

  it should "use the R1_Barcode_Bases/R2_Barcode_Bases columns for sample barcodes when useIndexColumns=false" in {
    val samples = CustomSampleSheet.maybeCustomSampleSheet(sampleSheet=this.sampleSheetPath, useIndexColumns=false)
    samples.foreach { sample =>
      val barcode = sample.sampleBarcode.value
      barcode.concatenatedBarcode.endsWith("GATTACAGA") shouldBe true
    }
  }

  "SampleBarcodeMetricCollector" should "fail if no sample barcodes are present in the read struture" in {
    an[Exception] should be thrownBy new SampleBarcodeMetricCollector(this.sampleSheet.toSeq, Seq(ReadStructure("50T"), ReadStructure("50T")))
    an[Exception] should be thrownBy new SampleBarcodeMetricCollector(this.sampleSheet.toSeq, Seq(ReadStructure("50T")))
    an[Exception] should be thrownBy new SampleBarcodeMetricCollector(this.sampleSheet.toSeq, Seq.empty)
  }

  it should "collect matching statistic on sample barcodes" in {
    val collector = new SampleBarcodeMetricCollector(this.sampleSheet.toSeq, Seq(ReadStructure("8B50T")))

    // Sample 1
    collector.increment(0, 0)
    collector.increment(0, 1)
    collector.increment(0, 1)
    collector.increment(0, 2)

    // Unmatched sample
    collector.increment(-1, 0)
    collector.increment(-1, 1)
    collector.increment(this.sampleSheet.size, 0) // Since there's an extra sample for the unmatched sample

    val metrics = collector.barcodeMetrics

    // Zero counts all around
    metrics.zip(sampleSheet.toSeq).slice(1, this.sampleSheet.size-1).foreach { case (metric, sample) =>
      metric shouldBe SampleBarcodeMetric(barcodeName=sample.sampleName, libraryName=sample.libraryId.value, barcode=sample.sampleBarcode.value.concatenatedBarcode)
    }

    // First sample
    val firstSample = this.sampleSheet.toSeq.headOption.value
    val firstMetric = metrics.headOption.value
    firstMetric.barcode_name            shouldBe firstSample.sampleName
    firstMetric.library_name            shouldBe firstSample.libraryId.value
    firstMetric.barcode                 shouldBe firstSample.sampleBarcode.value.concatenatedBarcode
    firstMetric.reads                   shouldBe 4
    firstMetric.pf_reads                shouldBe 4
    firstMetric.perfect_matches         shouldBe 1
    firstMetric.pf_perfect_matches      shouldBe 1
    firstMetric.one_mismatch_matches    shouldBe 2
    firstMetric.pf_one_mismatch_matches shouldBe 2

    // Unmatched sample
    val unmatchedMetric = metrics.lastOption.value
    unmatchedMetric.barcode_name shouldBe UnmatchedSampleId
    unmatchedMetric.library_name shouldBe UnmatchedSampleId
    unmatchedMetric.barcode      shouldBe "N"*8
    unmatchedMetric.reads        shouldBe 3
    unmatchedMetric.pf_reads     shouldBe 3
  }

    /** Helper method to create a [[FastqDemultiplexer]] */
  private def dx(structures: Seq[ReadStructure], minQ: Int = 0, maxQ: Int = SAMUtils.MAX_PHRED_SCORE, mm: Int = 2, md: Int = 1, mn: Int = 1): FastqDemultiplexer = {
    new FastqDemultiplexer(headers=this.headers, samples=this.sampleSheet.toSeq, readStructures=structures,
      minQ=minQ, maxQ=maxQ, maxMismatches=mm, minMismatchDelta=md, maxNoCalls=mn)
  }

  private def fq(name: String, bases: String, readNumber: Option[Int]=None): FastqRecord = FastqRecord(name=name, bases=bases, quals="I"*bases.length, comment=None, readNumber=readNumber)

  private def verifyFragUnmatchedSample(demuxRecord: DemuxRecord): Unit = {
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.getReadName shouldBe "frag"
    record.getReadString shouldBe "A"*100
    record.getBaseQualityString shouldBe "I"*100
    record.getReadPairedFlag shouldBe false
    record.getReadUnmappedFlag shouldBe true
    record.getReadGroup.getId shouldBe UnmatchedSampleId
    record.getAttribute("RX") shouldBe null
  }

  "FastqDemultiplexer" should "demux fragment reads" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")))
    val fastqRecord = fq(name="frag", bases="AAAAAAAAGATTACAGA" + "A"*100)

    val demuxRecord = demuxer.demultiplex(fastqRecord)
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.getReadName shouldBe "frag"
    record.getReadString shouldBe "A"*100
    record.getBaseQualityString shouldBe "I"*100
    record.getReadPairedFlag shouldBe false
    record.getReadUnmappedFlag shouldBe true
    record.getReadGroup.getId shouldBe "20000101-EXPID-1"
    record.getAttribute("RX") shouldBe null
  }

  it should "demux paired reads" in {
    val demuxer = dx(structures=Seq(ReadStructure("8B100T"),ReadStructure("9B100T")))
    val fq1 = fq(name="pair", bases="AAAAAAAA" + "A"*100, readNumber=Some(1))
    val fq2 = fq(name="pair", bases="GATTACAGA" + "T"*100, readNumber=Some(2))

    val demuxRecord = demuxer.demultiplex(fq1, fq2)
    demuxRecord.records.length shouldBe 2
    val r1 = demuxRecord.records.headOption.value
    val r2 = demuxRecord.records.lastOption.value

    r1.getReadName shouldBe "pair"
    r1.getReadString shouldBe "A"*100
    r1.getBaseQualityString shouldBe "I"*100
    r1.getReadPairedFlag shouldBe true
    r1.getFirstOfPairFlag shouldBe true
    r1.getReadUnmappedFlag shouldBe true
    r1.getMateUnmappedFlag shouldBe true
    r1.getReadGroup.getId shouldBe "20000101-EXPID-1"
    r1.getAttribute("RX") shouldBe null

    r2.getReadName shouldBe "pair"
    r2.getReadString shouldBe "T"*100
    r2.getBaseQualityString shouldBe "I"*100
    r2.getReadPairedFlag shouldBe true
    r2.getSecondOfPairFlag shouldBe true
    r2.getReadUnmappedFlag shouldBe true
    r2.getMateUnmappedFlag shouldBe true
    r2.getReadGroup.getId shouldBe "20000101-EXPID-1"
    r2.getAttribute("RX") shouldBe null
  }

  it should "set molecular barcodes" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B5M100T")))
    val fastqRecord = fq(name="frag", bases="AAAAAAAAGATTACAGA" + "NNNNN" + "A"*100)

    val demuxRecord = demuxer.demultiplex(fastqRecord)
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.getReadName shouldBe "frag"
    record.getReadString shouldBe "A"*100
    record.getBaseQualityString shouldBe "I"*100
    record.getReadPairedFlag shouldBe false
    record.getReadUnmappedFlag shouldBe true
    record.getReadGroup.getId shouldBe "20000101-EXPID-1"
    record.getAttribute("RX") shouldBe "NNNNN"
  }

  it should "assign to the 'unmatched' sample if there are too many mismatches" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=0) // no mismatches allowed
    val fastqRecord = fq(name="frag", bases="AAAAAAAAGATTACAGT" + "A"*100) // last sample barcode base mismatches
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  it should "assign to the 'unmatched' sample if the read matched two sample barcodes within the mismatch delta" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=10, md=3) // so many mismatches allowed!
    val fastqRecord = fq(name="frag", bases="GGGGGGTTGATTACAGA" + "A"*100) // matches the 4th barcode perfectly and the 3rd barcode with two mismatches
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  it should "assign to the 'unmatched' sample if the read's sample barcode has too many Ns" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=10, md=1, mn=0) // so many mismatches allowed!
    val fastqRecord = fq(name="frag", bases="GGGGGGTTGATTACAGN" + "A"*100) // one N
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  "FastqDemultiplexer.countMismatches" should "find no mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GATTACA".getBytes) shouldBe 0
  }

  it should "find two mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GACCACA".getBytes) shouldBe 2
  }

  it should "not count no calls" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GANNACA".getBytes) shouldBe 0
  }

  it should "find compare two sequences that have all mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "CTAATGT".getBytes) shouldBe 7
  }

  "DemuxFromInlineSampleBarcode" should "run end-to-end" in {
    val fastqs = new ListBuffer[FastqRecord]()
    fastqs += fq(name="frag1", bases="AAAAAAAAGATTACAGA" + "A"*100) // matches the first sample -> first sample
    fastqs += fq(name="frag2", bases="AAAAAAAAGATTACAGT" + "A"*100) // matches the first sample, one mismatch -> first sample
    fastqs += fq(name="frag3", bases="AAAAAAAAGATTACTTT" + "A"*100) // matches the first sample, three mismatches -> unmatched
    fastqs += fq(name="frag4", bases="GGGGGGTTGATTACAGA" + "A"*100) // matches the 4th barcode perfectly and the 3rd barcode with two mismatches, delta too small -> unmatched
    fastqs += fq(name="frag5", bases="AAAAAAAAGANNNNNNN" + "A"*100) // matches the first sample, too many Ns -> unmatched


    val fastq = makeTempFile("test", ".fastq")
    Io.writeLines(fastq, fastqs.map(_.toString))

    val output = {
      val dir = Files.createTempDirectory("DemuxFromInlineSampleBarcodeTest")
      dir.toFile.deleteOnExit()
      dir
    }

    val metrics = makeTempFile("metrics", ".txt")

    new DemuxFromInlineSampleBarcode(inputReadOne=Seq(fastq), output=output, sampleSheet=sampleSheetPath, readStructureReadOne=ReadStructure("17B100T"), metrics=metrics, maxMismatches=2, minMismatchDelta=3).execute()

    // Check the BAMs
    this.sampleSheet.foreach { sample =>
      val bam = sampleOutputBam(output=output, sample=sample)
      val records = readBamRecs(bam)

      if (sample.sampleOrdinal == 1) {
        records.length shouldBe 2
        records.map(_.getReadName) should contain theSameElementsInOrderAs Seq("frag1", "frag2")
      }
      else if (sample.sampleId == UnmatchedSampleId) {
        records.length shouldBe 3
        records.map(_.getReadName) should contain theSameElementsInOrderAs Seq("frag2", "frag4", "frag5")
      }
      else {
        records shouldBe 'empty
      }
    }

    // Check the metrics
    val sampleBarcodMetrics = Metric.read[SampleBarcodeMetric](metrics).toSeq

    sampleBarcodMetrics.length shouldBe this.samples.length

    sampleBarcodMetrics.zip(samples).foreach { case (metric, sample) =>
      metric.barcode_name shouldBe sample.sampleName
      metric.barcode      shouldBe sample.sampleBarcode.value.concatenatedBarcode
      metric.library_name shouldBe sample.libraryId.value
      if (sample.sampleOrdinal == 1) {
        metric.reads                   shouldBe 2
        metric.pf_reads                shouldBe 2
        metric.perfect_matches         shouldBe 1
        metric.one_mismatch_matches    shouldBe 1
        metric.pf_perfect_matches      shouldBe 1
        metric.pf_one_mismatch_matches shouldBe 1
      }
      else if (sample.sampleId == UnmatchedSampleId) {
        metric.reads                   shouldBe 3
        metric.pf_reads                shouldBe 3
        metric.perfect_matches         shouldBe 0
        metric.one_mismatch_matches    shouldBe 0
        metric.pf_perfect_matches      shouldBe 0
        metric.pf_one_mismatch_matches shouldBe 0
      }
      else {
        metric.reads                   shouldBe 0
        metric.pf_reads                shouldBe 0
        metric.perfect_matches         shouldBe 0
        metric.one_mismatch_matches    shouldBe 0
        metric.pf_perfect_matches      shouldBe 0
        metric.pf_one_mismatch_matches shouldBe 0
      }
    }
  }
}
