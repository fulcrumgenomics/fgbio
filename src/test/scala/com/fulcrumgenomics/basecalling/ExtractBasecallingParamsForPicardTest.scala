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

package com.fulcrumgenomics.basecalling

import java.nio.file.Files

import com.fulcrumgenomics.illumina.{Sample, SampleSheet}
import com.fulcrumgenomics.testing.{ErrorLogLevel, UnitSpec}
import com.fulcrumgenomics.util.Io

class ExtractBasecallingParamsForPicardTest extends UnitSpec with ErrorLogLevel {

  private val outputDir = Files.createTempDirectory("ExtractBasecallingParamsTest")

  "BasecallingParams.bamFileFrom" should "create a path to a BAM file with a single-indexed sample" in {
    val bam = BasecallingParams.bamFileFrom(
      output = outputDir,
      sample = new Sample(sampleOrdinal=0, sampleId="sampleId", sampleName="sampleName", libraryId="libraryId", i7IndexBases=Some("GATTACA")),
      lane   = 1
    )
    bam.toString shouldBe outputDir.resolve("sampleName.GATTACA.1.bam").toString
  }

  "BasecallingParams.unmatchedBamFileFrom" should "create a path to an unmatched BAM" in {
    val bam = BasecallingParams.unmatchedBamFileFrom(output = outputDir, lane = 1)
    bam.toString shouldBe outputDir.resolve("unmatched.1.bam").toString
  }

  it should "create a path to a BAM file without a library identifier" in {
    val bam = BasecallingParams.bamFileFrom(
      output = outputDir,
      sample = new Sample(sampleOrdinal=0, sampleId="sampleId", sampleName="sampleName", libraryId="libraryId", i7IndexBases=Some("GATTACA"), i5IndexBases=Some("ACATTAG")),
      lane   = 2
    )
    bam.toString shouldBe outputDir.resolve("sampleName.GATTACAACATTAG.2.bam").toString
  }

  private val singleIndexSampleSheet =
    """[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,R1_Barcode_Bases,R2_Barcode_Bases,I7_Index_ID,index,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,GATTACAG,GATTACAGA,I7_1,GATTACAACGT,Sample_Project_1,Description_1
      |20000101-EXPID-2,Sample_Name_2,,,GATTACAG,GATTACAGA,I7_2,GATTACAACGT,Sample_Project_2,Description_2
      |20000101-EXPID-3,Sample_Name_3,,,GATTACAG,GATTACAGA,I7_3,GATTACAACGT,Sample_Project_3,Description_3
      |20000101-EXPID-4,Sample_Name_4,,,GATTACAG,GATTACAGA,I7_4,GATTACAACGT,Sample_Project_4,Description_4
      |20000101-EXPID-5,Sample_Name_5,,,GATTACAG,GATTACAGA,I7_5,GATTACAACGT,Sample_Project_5,Description_5
      |20000101-EXPID-6,Sample_Name_6,,,GATTACAG,GATTACAGA,I7_6,GATTACAACGT,Sample_Project_6,Description_6
      |20000101-EXPID-7,Sample_Name_7,,,GATTACAG,GATTACAGA,I7_7,GATTACAACGT,Sample_Project_7,Description_7
      |20000101-EXPID-8,Sample_Name_8,,,GATTACAG,GATTACAGA,I7_8,GATTACAACGT,Sample_Project_8,Description_8
      |20000101-EXPID-9,Sample_Name_9,,,GATTACAG,GATTACAGA,I7_9,GATTACAACGT,Sample_Project_9,Description_9
      |20000101-EXPID-10,Sample_Name_10,,,GATTACAG,GATTACAGA,I7_10,GATTACAACGT,Sample_Project_10,Description_10
      |20000101-EXPID-11,Sample_Name_11,,,GATTACAG,GATTACAGA,I7_11,GATTACAACGT,Sample_Project_11,Description_11
      |20000101-EXPID-12,Sample_Name_12,,,GATTACAG,GATTACAGA,I7_12,GATTACAACGT,Sample_Project_12,Description_12"""
      .stripMargin.split("\n").toIndexedSeq

  private val dualIndexedSampleSheet =
    """[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,R1_Barcode_Bases,R2_Barcode_Bases,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,GATTACAG,GATTACAGA,I7_1,GATTACAACGT,I5_1,GATTACA,Sample_Project_1,Description_1
      |20000101-EXPID-2,Sample_Name_2,,,GATTACAG,GATTACAGA,I7_2,GATTACAACGT,I5_2,GATTACA,Sample_Project_2,Description_2
      |20000101-EXPID-3,Sample_Name_3,,,GATTACAG,GATTACAGA,I7_3,GATTACAACGT,I5_3,GATTACA,Sample_Project_3,Description_3
      |20000101-EXPID-4,Sample_Name_4,,,GATTACAG,GATTACAGA,I7_4,GATTACAACGT,I5_4,GATTACA,Sample_Project_4,Description_4
      |20000101-EXPID-5,Sample_Name_5,,,GATTACAG,GATTACAGA,I7_5,GATTACAACGT,I5_5,GATTACA,Sample_Project_5,Description_5
      |20000101-EXPID-6,Sample_Name_6,,,GATTACAG,GATTACAGA,I7_6,GATTACAACGT,I5_6,GATTACA,Sample_Project_6,Description_6
      |20000101-EXPID-7,Sample_Name_7,,,GATTACAG,GATTACAGA,I7_7,GATTACAACGT,I5_7,GATTACA,Sample_Project_7,Description_7
      |20000101-EXPID-8,Sample_Name_8,,,GATTACAG,GATTACAGA,I7_8,GATTACAACGT,I5_8,GATTACA,Sample_Project_8,Description_8
      |20000101-EXPID-9,Sample_Name_9,,,GATTACAG,GATTACAGA,I7_9,GATTACAACGT,I5_9,GATTACA,Sample_Project_9,Description_9
      |20000101-EXPID-10,Sample_Name_10,,,GATTACAG,GATTACAGA,I7_10,GATTACAACGT,I5_10,GATTACA,Sample_Project_10,Description_10
      |20000101-EXPID-11,Sample_Name_11,,,GATTACAG,GATTACAGA,I7_11,GATTACAACGT,I5_11,GATTACA,Sample_Project_11,Description_11
      |20000101-EXPID-12,Sample_Name_12,,,GATTACAG,GATTACAGA,I7_12,GATTACAACGT,I5_12,GATTACA,Sample_Project_12,Description_12"""
  .stripMargin.split("\n").toIndexedSeq

  private val dualIndexedSampleSheetNoProjectOrDescription =
    """[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,R1_Barcode_Bases,R2_Barcode_Bases,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,GATTACAG,GATTACAGA,I7_1,GATTACAACGT,I5_1,GATTACA,,"""
  .stripMargin.split("\n").toIndexedSeq

  private val duplicateNamesSampleSheet =
    """[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,R1_Barcode_Bases,R2_Barcode_Bases,I7_Index_ID,index,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,GATTACAG,GATTACAGA,I7_1,GATTACAACGT,Sample_Project_1,Description_1
      |20000101-EXPID-2,Sample_Name_1,,,GATTACAG,GATTACAGA,I7_2,GATTACAACGT,Sample_Project_2,Description_2
      |20000101-EXPID-3,Sample_Name_3,,,GATTACAG,GATTACAGA,I7_3,GATTACAACGT,Sample_Project_3,Description_3
      |20000101-EXPID-4,Sample_Name_4,,,GATTACAG,GATTACAGA,I7_4,GATTACAACGT,Sample_Project_4,Description_4
      |20000101-EXPID-5,Sample_Name_5,,,GATTACAG,GATTACAGA,I7_5,GATTACAACGT,Sample_Project_5,Description_5
      |20000101-EXPID-6,Sample_Name_6,,,GATTACAG,GATTACAGA,I7_6,GATTACAACGT,Sample_Project_6,Description_6
      |20000101-EXPID-7,Sample_Name_7,,,GATTACAG,GATTACAGA,I7_7,GATTACAACGT,Sample_Project_7,Description_7
      |20000101-EXPID-8,Sample_Name_8,,,GATTACAG,GATTACAGA,I7_8,GATTACAACGT,Sample_Project_8,Description_8
      |20000101-EXPID-9,Sample_Name_9,,,GATTACAG,GATTACAGA,I7_9,GATTACAACGT,Sample_Project_9,Description_9
      |20000101-EXPID-10,Sample_Name_10,,,GATTACAG,GATTACAGA,I7_10,GATTACAACGT,Sample_Project_10,Description_10
      |20000101-EXPID-11,Sample_Name_11,,,GATTACAG,GATTACAGA,I7_11,GATTACAACGT,Sample_Project_11,Description_11
      |20000101-EXPID-12,Sample_Name_12,,,GATTACAG,GATTACAGA,I7_12,GATTACAACGT,Sample_Project_12,Description_12"""
      .stripMargin.split("\n").toIndexedSeq


  "BasecallingParams.from" should "extract params from a single-index sequencing run" in {
    val sampleSheet = SampleSheet(singleIndexSampleSheet.toIterator, lane=None)
    val params = BasecallingParams.from(sampleSheet=sampleSheet, lanes=Seq(1), output=outputDir)
    params should have size 1
    val param = params.head

    param.barcodeFile shouldBe BasecallingParams.barcodeFileFrom(output=outputDir, lane=1)
    param.libraryParamsFile shouldBe BasecallingParams.libraryParamsFileFrom(output=outputDir, lane=1)
    param.bams.head shouldBe BasecallingParams.bamFileFrom(output=outputDir, sample=sampleSheet.head, lane=1)

    // Check the header, first sample, and last line (last sample)
    val barcodeParams = Io.readLines(param.barcodeFile).toSeq
    barcodeParams.head shouldBe "barcode_sequence_1\tbarcode_name\tlibrary_name"
    barcodeParams.drop(1).head shouldBe "GATTACAACGT\tGATTACAACGT\t20000101-EXPID-1"
    barcodeParams.last shouldBe "GATTACAACGT\tGATTACAACGT\t20000101-EXPID-12"

    // Check the header, first sample, and last line (unmatched sample)
    val libraryParams = Io.readLines(param.libraryParamsFile).toSeq
    libraryParams.head shouldBe "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tDS"
    libraryParams.drop(1).head shouldBe s"GATTACAACGT\t${outputDir.resolve("Sample_Name_1.GATTACAACGT.1.bam")}\tSample_Name_1\t20000101-EXPID-1\tDescription_1"
    libraryParams.last shouldBe s"N\t${outputDir.resolve("unmatched.1.bam")}\tunmatched\tunmatched\tunmatched"
  }

  it should "extract params for multiple lanes from a sequencing run" in {
    val sampleSheet = SampleSheet(singleIndexSampleSheet.toIterator, lane=None)
    val params = BasecallingParams.from(sampleSheet=sampleSheet, lanes=Seq(1, 2, 4), output=outputDir)
    params should have size 3
    val param = params.last

    param.barcodeFile shouldBe BasecallingParams.barcodeFileFrom(output=outputDir, lane=4)
    param.libraryParamsFile shouldBe BasecallingParams.libraryParamsFileFrom(output=outputDir, lane=4)
    param.bams.head shouldBe BasecallingParams.bamFileFrom(output=outputDir, sample=sampleSheet.head, lane=4)

    // Check the header, first sample, and last line (last sample)
    val barcodeParams = Io.readLines(param.barcodeFile).toSeq
    barcodeParams.head shouldBe "barcode_sequence_1\tbarcode_name\tlibrary_name"
    barcodeParams.drop(1).head shouldBe "GATTACAACGT\tGATTACAACGT\t20000101-EXPID-1"
    barcodeParams.last shouldBe "GATTACAACGT\tGATTACAACGT\t20000101-EXPID-12"

    // Check the header, first sample, and last line (unmatched sample)
    val libraryParams = Io.readLines(param.libraryParamsFile).toSeq
    libraryParams.head shouldBe "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tDS"
    libraryParams.drop(1).head shouldBe s"GATTACAACGT\t${outputDir.resolve("Sample_Name_1.GATTACAACGT.4.bam")}\tSample_Name_1\t20000101-EXPID-1\tDescription_1"
    libraryParams.last shouldBe s"N\t${outputDir.resolve("unmatched.4.bam")}\tunmatched\tunmatched\tunmatched"
  }

  it should "extract params from a dual-indexed sequencing run" in {
    val sampleSheet = SampleSheet(dualIndexedSampleSheet.toIterator, lane=None)
    val params = BasecallingParams.from(sampleSheet=sampleSheet, lanes=Seq(1), output=outputDir)
    params should have size 1
    val param = params.head

    param.barcodeFile shouldBe BasecallingParams.barcodeFileFrom(output=outputDir, lane=1)
    param.libraryParamsFile shouldBe BasecallingParams.libraryParamsFileFrom(output=outputDir, lane=1)
    param.bams.head shouldBe BasecallingParams.bamFileFrom(output=outputDir, sample=sampleSheet.head, lane=1)

    // Check the header, first sample, and last line (last sample)
    val barcodeParams = Io.readLines(param.barcodeFile).toSeq
    barcodeParams.head shouldBe "barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name"
    barcodeParams.drop(1).head shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-1"
    barcodeParams.last shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-12"

    // Check the header, first sample, and last line (unmatched sample)
    val libraryParams = Io.readLines(param.libraryParamsFile).toSeq
    libraryParams.head shouldBe "BARCODE_1\tBARCODE_2\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tDS"
    libraryParams.drop(1).head shouldBe s"GATTACAACGT\tGATTACA\t${outputDir.resolve("Sample_Name_1.GATTACAACGTGATTACA.1.bam")}\tSample_Name_1\t20000101-EXPID-1\tDescription_1"
    libraryParams.last shouldBe s"N\tN\t${outputDir.resolve("unmatched.1.bam")}\tunmatched\tunmatched\tunmatched"
  }

  it should "exclude the description if the project and description are not defined on all samples" in {
    val sampleSheet = SampleSheet(dualIndexedSampleSheetNoProjectOrDescription.toIterator, lane=None)
    val params = BasecallingParams.from(sampleSheet=sampleSheet, lanes=Seq(1), output=outputDir)
    params should have size 1
    val param = params.head

    param.barcodeFile shouldBe BasecallingParams.barcodeFileFrom(output=outputDir, lane=1)
    param.libraryParamsFile shouldBe BasecallingParams.libraryParamsFileFrom(output=outputDir, lane=1)
    param.bams.head shouldBe BasecallingParams.bamFileFrom(output=outputDir, sample=sampleSheet.head, lane=1)

    // Check the header, first sample, and last line (last sample)
    val barcodeParams = Io.readLines(param.barcodeFile).toSeq
    barcodeParams.head shouldBe "barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name"
    barcodeParams.drop(1).head shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-1"
    barcodeParams.last shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-1"

    // Check the header, first sample, and last line (unmatched sample)
    val libraryParams = Io.readLines(param.libraryParamsFile).toSeq
    libraryParams.head shouldBe "BARCODE_1\tBARCODE_2\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME"
    libraryParams.drop(1).head shouldBe s"GATTACAACGT\tGATTACA\t${outputDir.resolve("Sample_Name_1.GATTACAACGTGATTACA.1.bam")}\tSample_Name_1\t20000101-EXPID-1"
    libraryParams.last shouldBe s"N\tN\t${outputDir.resolve("unmatched.1.bam")}\tunmatched\tunmatched"
  }

  "ExtractBasecallingParams" should "run end-to-end" in {
    val sampleSheet = makeTempFile("SampleSheet", ".csv")
    Io.writeLines(sampleSheet, dualIndexedSampleSheet.toSeq)
    new ExtractBasecallingParamsForPicard(input=sampleSheet, output=outputDir, lanes=Seq(1)).execute()

    val barcodeFile       = BasecallingParams.barcodeFileFrom(output=outputDir, lane=1)
    val libraryParamsFile = BasecallingParams.libraryParamsFileFrom(output=outputDir, lane=1)

    // Check the header, first sample, and last line (last sample)
    val barcodeParams = Io.readLines(barcodeFile).toSeq
    barcodeParams.head shouldBe "barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name"
    barcodeParams.drop(1).head shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-1"
    barcodeParams.last shouldBe "GATTACAACGT\tGATTACA\tGATTACAACGTGATTACA\t20000101-EXPID-12"

    // Check the header, first sample, and last line (unmatched sample)
    val libraryParams = Io.readLines(libraryParamsFile).toSeq
    libraryParams.head shouldBe "BARCODE_1\tBARCODE_2\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tDS"
    libraryParams.drop(1).head shouldBe s"GATTACAACGT\tGATTACA\t${outputDir.resolve("Sample_Name_1.GATTACAACGTGATTACA.1.bam")}\tSample_Name_1\t20000101-EXPID-1\tDescription_1"
    libraryParams.last shouldBe s"N\tN\t${outputDir.resolve("unmatched.1.bam")}\tunmatched\tunmatched\tunmatched"
  }
}
