/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.ConsensusTags.PerRead.{AbRawReadCount, BaRawReadCount, RawReadCount}
import org.scalatest.OptionValues

class UmisTest extends UnitSpec with OptionValues {
  import Umis.copyUmiFromReadName

  private val builder = new SamBuilder()
  private def rec(name: String): SamRecord = builder.addFrag(name=name, unmapped=true).value

  implicit private class SamRecordUmi(rec: SamRecord) {
    def nameAndUmi: (String, String) = (rec.name, rec[String](ConsensusTags.UmiBases))
  }

  "Umis.extractUmisFromReadName" should "extract a UMI from a well-formed read name" in {
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:ACGTACGT", strict=true).value shouldBe "ACGTACGT"
  }

  it should "return None if the read name has only 7 parts in strict mode" in {
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7", strict=true) shouldBe None
  }

  it should "throw an exception in strict mode if the read has too many or too few segments" in {
    an[Exception] shouldBe thrownBy { Umis.extractUmisFromReadName("1:2:3:4:5:6", strict=true) }
    an[Exception] shouldBe thrownBy { Umis.extractUmisFromReadName("1:2:3:4:5:6:7:8:ACGT", strict=true) }
  }

  it should "translate pluses to hyphens when multiple UMIs are present" in {
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:ACGTACGT+TTGCGGCT", strict=true).value shouldBe "ACGTACGT-TTGCGGCT"
  }

  it should "extract a UMI from the last segment in non-strict mode, if it looks like a UMI" in {
    Umis.extractUmisFromReadName("1:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:5:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:5:6:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:8:ACGT", strict=false).value shouldBe "ACGT"
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:8:9:ACGT", strict=false).value shouldBe "ACGT"
  }

  it should "return None in non-strict mode if the last segment doesn't look like a UMI" in {
    Umis.extractUmisFromReadName("1:2", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3:4", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3:4:5", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3:4:5:6", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7", strict=false) shouldBe None
    Umis.extractUmisFromReadName("1:2:3:4:5:6:7:8", strict=false) shouldBe None
  }


  "Umis.copyUmiFromReadName" should "copy the UMI from the read name" in {
    copyUmiFromReadName(rec=rec("UMI:A")).nameAndUmi shouldBe ("UMI:A", "A")
    copyUmiFromReadName(rec=rec("UMI:C:A")).nameAndUmi shouldBe ("UMI:C:A", "A")
    copyUmiFromReadName(rec=rec("UMI:C:ACC-GGT")).nameAndUmi shouldBe ("UMI:C:ACC-GGT", "ACC-GGT")
  }

  it should "remove the UMI if specified" in {
    copyUmiFromReadName(rec=rec("UMI:A"), removeUmi=true).nameAndUmi shouldBe ("UMI", "A")
    copyUmiFromReadName(rec=rec("UMI:C:A"), removeUmi=true).nameAndUmi shouldBe ("UMI:C", "A")
    copyUmiFromReadName(rec=rec("UMI:C:ACC+GGT"), removeUmi=true).nameAndUmi shouldBe ("UMI:C", "ACC-GGT")
  }
  
  it should "split on a different name delimiter if specified" in {
    copyUmiFromReadName(rec=rec("UMI-A"), fieldDelimiter='-').nameAndUmi shouldBe ("UMI-A", "A")
    copyUmiFromReadName(rec=rec("UMI-C-A"), fieldDelimiter='-').nameAndUmi shouldBe ("UMI-C-A", "A")
    copyUmiFromReadName(rec=rec("UMI-C-ACC+GGT"), fieldDelimiter='-').nameAndUmi shouldBe ("UMI-C-ACC+GGT", "ACC-GGT")
  }

  it should "change the UMI delimiter from + to -" in {
    copyUmiFromReadName(rec=rec("UMI:C:ACC+GGT")).nameAndUmi shouldBe ("UMI:C:ACC+GGT", "ACC-GGT")
  }

  it should "fail if the read name has only one field" in {
    an[Exception] should be thrownBy copyUmiFromReadName(rec=rec("NAME"))
    an[Exception] should be thrownBy copyUmiFromReadName(rec=rec("1-2"))
  }

  it should "fail if UMI contains illegal characters" in {
    an[Exception] should be thrownBy copyUmiFromReadName(rec=rec("NAME:XYZ"))
    an[Exception] should be thrownBy copyUmiFromReadName(rec=rec("NAME:ACGT-CCKC"))
    an[Exception] should be thrownBy copyUmiFromReadName(rec=rec("NAME:CCKC-ACGT"))
  }

  "Umis.isFgbioStyleConsensus" should "return false for reads without consensus tags" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    Umis.isFgbioStyleConsensus(builder.addFrag(start=100).value) shouldBe false
    val pair = builder.addPair(start1=100, start2=100, unmapped2=true)
    pair.length shouldBe 2
    pair.forall(rec => Umis.isFgbioStyleConsensus(rec)) shouldBe false
  }

  it should "return true for reads with consensus tags" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    Umis.isFgbioStyleConsensus(builder.addFrag(start=10, attrs=Map(RawReadCount -> 10)).value) shouldBe true
    Umis.isFgbioStyleConsensus(builder.addFrag(start=10, attrs=Map(AbRawReadCount -> 10, BaRawReadCount -> 10)).value) shouldBe true
  }

  "Umis.isFgbioSimplexConsensus" should "return true for reads with simplex only consensus tags" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    Umis.isFgbioSimplexConsensus(builder.addFrag(start=100).value) shouldBe false
    Umis.isFgbioSimplexConsensus(builder.addFrag(start=10, attrs=Map(RawReadCount -> 10)).value) shouldBe true
    Umis.isFgbioSimplexConsensus(builder.addFrag(start=10, attrs=Map(AbRawReadCount -> 10, BaRawReadCount -> 10)).value) shouldBe false
    Umis.isFgbioSimplexConsensus(builder.addFrag(start=10, attrs=Map(RawReadCount -> 20, AbRawReadCount -> 10, BaRawReadCount -> 10)).value) shouldBe false

    val pair = builder.addPair(start1=100, start2=100, unmapped2=true)
    pair.length shouldBe 2
    pair.forall(rec => Umis.isFgbioSimplexConsensus(rec)) shouldBe false
  }

  "Umis.isFgbioDuplexConsensus" should "return true for reads with duplex only consensus tags" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    Umis.isFgbioDuplexConsensus(builder.addFrag(start=100).value) shouldBe false
    Umis.isFgbioDuplexConsensus(builder.addFrag(start=10, attrs=Map(RawReadCount -> 10)).value) shouldBe false
    Umis.isFgbioDuplexConsensus(builder.addFrag(start=10, attrs=Map(AbRawReadCount -> 10, BaRawReadCount -> 10)).value) shouldBe true
    Umis.isFgbioDuplexConsensus(builder.addFrag(start=10, attrs=Map(RawReadCount -> 20, AbRawReadCount -> 10, BaRawReadCount -> 10)).value) shouldBe true

    val pair = builder.addPair(start1=100, start2=100, unmapped2=true)
    pair.length shouldBe 2
    pair.forall(rec => Umis.isFgbioDuplexConsensus(rec)) shouldBe false
  }
}
