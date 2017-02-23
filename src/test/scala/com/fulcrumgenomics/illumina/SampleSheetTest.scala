/**
  * Copyright (c) 2016, Fulcrum Genomics LLC
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *
  * 1. Redistributions of source code must retain the above copyright notice,
  * this list of conditions and the following disclaimer.
  *
  * 2. Redistributions in binary form must reproduce the above copyright notice,
  * this list of conditions and the following disclaimer in the documentation
  * and/or other materials provided with the distribution.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */

package com.fulcrumgenomics.illumina

import java.nio.file.Paths

import org.scalatest.{FlatSpec, Matchers, OptionValues}

class SampleSheetTest extends FlatSpec with Matchers with OptionValues {
  private val testDir = Paths.get("src/test/resources/com/fulcrumgenomics/util/miseq/")

  "SampleSheet" should "load a simple sample sheet" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.csv"))
    sampleSheet.size shouldBe 12
    sampleSheet.foreach { sample =>
      sample.sampleId                    shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.sampleName                  shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.libraryId                   shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.project.value               shouldBe "Sample_Project_" + sample.sampleOrdinal
      sample.description.value           shouldBe "Description_" + sample.sampleOrdinal
      sample.i7IndexBases.value          shouldBe "GATTACAACGT"
      sample.i5IndexBases.value          shouldBe "GATTACA"
      sample.extendedAttributes.keySet should contain ("R1_Barcode_Bases".toUpperCase)
      sample.extendedAttributes.keySet should contain ("R2_Barcode_Bases".toUpperCase)
    }
  }

  it should "load a sample sheet with only required values" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheetOnlyRequired.csv"))
    sampleSheet.size shouldBe 12
    sampleSheet.foreach { sample =>
      sample.sampleId             shouldBe "20000101-EXPID-" + sample.sampleOrdinal
      sample.sampleName           shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.libraryId            shouldBe "Sample_Name_" + sample.sampleOrdinal
      sample.project.value        shouldBe "Sample_Project_" + sample.sampleOrdinal
      sample.description.value    shouldBe "Description_" + sample.sampleOrdinal
      sample.i7IndexBases         shouldBe None
      sample.i5IndexBases         shouldBe None
      sample.extendedAttributes.keySet should not contain "R1_Barcode_Bases"
      sample.extendedAttributes.keySet should not contain "R2_Barcode_Bases"
    }
  }

  it should "load a sample sheet and a subset of samples by lane" in {
    Seq(1, 2, 3).foreach { lane =>
      val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=Some(lane))
      sampleSheet.foreach { sample => sample.lane.value shouldBe lane }
      sampleSheet.size shouldBe 4
    }
    SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=Some(4)) shouldBe 'empty
    SampleSheet(testDir.resolve("SampleSheet.lanes.csv"), lane=None).size shouldBe 12
  }

  it should "throw an exception when no sample section is found" in {
    an[Exception] should be thrownBy SampleSheet(testDir.resolve("SampleSheetNoSampleSection.csv"))
  }

  it should "throw an exception when a row in the sample section has fewer entries than expected" in {
    an[Exception] should be thrownBy SampleSheet(testDir.resolve("SampleSheetSampleMissingColumns.csv"))
  }

  "SampleSheet.get" should "get the ith sample" in {
    val sampleSheet = SampleSheet(testDir.resolve("SampleSheet.csv"))
    val samples = sampleSheet.toList
    samples.zipWithIndex.foreach { case (sample, i) =>
        sampleSheet.get(i) shouldBe samples(i)
    }
  }

  "SampleSheet.cleanMiseqSampleId" should "clean [#] to ' ' and [_+ ] to '-'" in {
    SampleSheet.cleanMiseqSampleId("ABCDEFGHIJKLMNOP") shouldBe "ABCDEFGHIJKLMNOP"
    SampleSheet.cleanMiseqSampleId("#") shouldBe " "
    SampleSheet.cleanMiseqSampleId("_+ ") shouldBe "---"
    SampleSheet.cleanMiseqSampleId("ABCDE#_+ ") shouldBe "ABCDE ---"
  }

  "SampleSheet.SplitRegex" should "split by commas but commas in quotes" in {
    "A,B,C".split(",") should have size 3
    "A,B,C".split(SampleSheet.SplitRegex) should have size 3
    "A,\"B,C\",D".split(SampleSheet.SplitRegex) should have size 3
    "\"B,C\",D,E".split(SampleSheet.SplitRegex) should have size 3
    "A,B,\"C,D\"".split(SampleSheet.SplitRegex) should have size 3
  }
}

