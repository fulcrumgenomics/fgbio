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

package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec

class QualityEncodingTest extends UnitSpec {
  private val ExpectedQuals = Array(2,3,4,5,10,15,20,25,30,35,40,45,50).map(_.toByte)
  private val ExpectedAscii = ExpectedQuals.map(b => (b + 33).toChar).mkString

  "Solexa quality encoding" should "decode some qualities correctly" in {
    QualityEncoding.Solexa.toStandardNumeric(">@BDJOTY^chmr") shouldBe ExpectedQuals
  }

  it should "translate some qualities into standard ascii correctly" in {
    QualityEncoding.Solexa.toStandardAscii(">@BDJOTY^chmr") shouldBe ExpectedAscii
  }

  "Illumina quality encoding" should "decode some qualities correctly" in {
    QualityEncoding.Illumina.toStandardNumeric("BCDEJOTY^chmr") shouldBe ExpectedQuals
  }

  it should "translate some qualities into standard ascii correctly" in {
    QualityEncoding.Illumina.toStandardAscii("BCDEJOTY^chmr") shouldBe ExpectedAscii
  }

  "Standard quality encoding" should "decode some qualities correctly" in {
    QualityEncoding.Standard.toStandardNumeric("#$%&+05:?DINS") shouldBe ExpectedQuals
  }

  it should "return an identical string when translateing to standard ascii" in {
    QualityEncoding.Standard.toStandardAscii("#$%&+05:?DINS") shouldBe ExpectedAscii
  }

  "QualityEncodingDetector" should "find every encoding compatible with no evidence" in {
    val detector = new QualityEncodingDetector()
    detector.compatibleEncodings should have length 3
    detector.rankedCompatibleEncodings(q=30) should have length 3
    detector.isAmbiguous shouldBe true
    QualityEncoding.all.forall(detector.isCompatible) shouldBe true
  }

  it should "determine Standard format from qualities in range 2-40" in {
    val detector = new QualityEncodingDetector()
    Range.inclusive(2, 40).foreach(q => detector.add((q+33).toChar))
    detector.isAmbiguous shouldBe false
    detector.compatibleEncodings shouldBe Seq(QualityEncoding.Standard)
    detector.rankedCompatibleEncodings(q=30) shouldBe Seq(QualityEncoding.Standard)
  }

  it should "match multiple encodings and resolve to Solexa by ranking" in {
    val detector = new QualityEncodingDetector()
    Range.inclusive(-5, 40).foreach(q => detector.add((q+64).toChar))
    detector.isAmbiguous shouldBe true
    detector.compatibleEncodings should contain theSameElementsAs Seq(QualityEncoding.Solexa, QualityEncoding.Standard)
    detector.rankedCompatibleEncodings(q=30) shouldBe Seq(QualityEncoding.Solexa, QualityEncoding.Standard)
  }

  it should "match multiple encodings and resolve to Standard by ranking" in {
    val detector = new QualityEncodingDetector()
    Range.inclusive(35, 40).foreach(q => detector.add((q+33).toChar))
    detector.isAmbiguous shouldBe true
    detector.rankedCompatibleEncodings(q=30).head shouldBe QualityEncoding.Standard
  }

  it should "sample reads in order to determine quality encoding" in {
    val qualStrings =
    """
      |AA9-AEF9,EE9FF9C,EFF999,C,EF9,@,EFF99F,,,C,C99E,,6CFF9CEEFG9,@F,,C,,,,,C,C,
      |-ACCCFFGCFGGG9FFFGG,EFC9F,CFFF9<F9,,,;,,,CFEFCF9,CFGGGG9@F,,,CFEF,,,;CE<EEF
      |A-6ABEF9FF:C7@BEF8,,DFC88CC,,C,CF@FGFEEE7CCC:CF7FE,B,CEE9FF,C,CEEEEEE,CE8,C
      |A<@@-<-C<;FFFG,++6:,@,,,,CCDEF,,,,,,6,,,,,C9CEE9C@@,C<,,,,CC,,,,+,,:,:,:+CC
      |@C-BC9FE9<,CE9<,,6CFC9,,CEEDF,@<<C+++@CFG9,,9,,<,CEF<FED,@E9E9,,CF9FE,,,,,6
      |ACC9<D@;8@FEEF>FEG@,8E8,,,,,,;,CEFF88,,CF8,C,,+CF89,88E9,C9FF,,CE8,,CEEGG9C
      |-<BCCGGGFCAC<<<,CE<FF8@E<<,,,B88E@F<9;9E<FFFDEC8E8FFFGG8E,,,@@@88CE8E,,,C9<
      |-A<9@F@EGGFGF<FF9EC,++F7:+CFFFGGFFGCFFFCFGGGA8F<7F,6C++@F7F89C@@BFG9=C+=FE,
    """.stripMargin.trim.lines

    val detector = new QualityEncodingDetector
    detector.sample(qualStrings)
    detector.rankedCompatibleEncodings(30) shouldBe Seq(QualityEncoding.Standard)
  }
}
