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

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class UmisTest extends UnitSpec with OptionValues {
  import Umis.copyUmiFromReadName

  private val builder = new SamBuilder()
  private def rec(name: String): SamRecord = builder.addFrag(name=name, unmapped=true).value

  implicit private class SamRecordUmi(rec: SamRecord) {
    def nameAndUmi: (String, String) = (rec.name, rec[String](ConsensusTags.UmiBases))
  }

  "copyUmiFromReadName" should "copy the UMI from the read name" in {
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
    copyUmiFromReadName(rec=rec("UMI-A"), nameDelimiter='-').nameAndUmi shouldBe ("UMI-A", "A")
    copyUmiFromReadName(rec=rec("UMI-C-A"), nameDelimiter='-').nameAndUmi shouldBe ("UMI-C-A", "A")
    copyUmiFromReadName(rec=rec("UMI-C-ACC+GGT"), nameDelimiter='-').nameAndUmi shouldBe ("UMI-C-ACC+GGT", "ACC-GGT")
  }

  it should "change the UMI delimiter if specified" in {
    copyUmiFromReadName(rec=rec("UMI:C:ACC+GGT")).nameAndUmi shouldBe ("UMI:C:ACC+GGT", "ACC-GGT")
    copyUmiFromReadName(rec=rec("UMI:C:ACCXGGT"), umiDelimiter=Some('X')).nameAndUmi shouldBe ("UMI:C:ACCXGGT", "ACC-GGT")
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
}
