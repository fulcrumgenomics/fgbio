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
 */

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.ReadStructure.SubReadWithQuals
import com.fulcrumgenomics.util.SegmentType._
import org.scalatest.OptionValues

import scala.util.Try

class ReadStructureTest extends UnitSpec with OptionValues {

  private def compareReadStructures(actual: ReadStructure, expected: Seq[ReadSegment]): Unit = {
    // make sure the segments match
    actual shouldBe expected
    // make sure the string representations are the same
    actual.toString shouldBe ReadStructure(expected).toString
  }

  private def compareReadStructuresResetOffset(actual: Seq[ReadSegment], expected: ReadStructure): Unit = {
    val actualReadStructure = ReadStructure(actual, resetOffsets=true)
    compareReadStructures(actualReadStructure, expected)
  }

  private def T(off: Int, len: Int) = ReadSegment(offset=off, length=Some(len), kind=SegmentType.Template)
  private def B(off: Int, len: Int) = ReadSegment(offset=off, length=Some(len), kind=SegmentType.SampleBarcode)
  private def M(off: Int, len: Int) = ReadSegment(offset=off, length=Some(len), kind=SegmentType.MolecularBarcode)
  private def S(off: Int, len: Int) = ReadSegment(offset=off, length=Some(len), kind=SegmentType.Skip)

  "ReadStructure" should "be built from a string" in {
    compareReadStructures(ReadStructure("1T"), Seq(T(0, 1)))
    compareReadStructures(ReadStructure("1B"), Seq(B(0, 1)))
    compareReadStructures(ReadStructure("1M"), Seq(M(0, 1)))
    compareReadStructures(ReadStructure("1S"), Seq(S(0, 1)))
    compareReadStructures(ReadStructure("101T"), Seq(T(0, 101)))
    compareReadStructures(ReadStructure("5B101T"), Seq(B(0, 5), T(5, 101)))
    compareReadStructures(ReadStructure("123456789T"), Seq(T(0, 123456789)))
    compareReadStructures(ReadStructure("10T10B10B10S10M"), Seq(T(0, 10), B(10, 10), B(20, 10), S(30, 10), M(40, 10)))
  }

  it should "allow + only once and only for the last segment of the read" in {
    import SegmentType._
    ReadStructure("5M+T") shouldBe Seq(ReadSegment(0, Some(5), MolecularBarcode), ReadSegment(5, None, Template))
    ReadStructure("+M")   shouldBe Seq(ReadSegment(0, None, MolecularBarcode))
    an[Exception] shouldBe thrownBy { ReadStructure("++M") }
    an[Exception] shouldBe thrownBy { ReadStructure("5M++T") }
    an[Exception] shouldBe thrownBy { ReadStructure("5M70+T") }
    an[Exception] shouldBe thrownBy { ReadStructure("+M+T") }
    an[Exception] shouldBe thrownBy { ReadStructure("+M70T") }
  }
  
  it should "be built from segments while resetting their offset" in {
    compareReadStructuresResetOffset(Seq(T(Int.MaxValue, 1)), ReadStructure("1T"))
    compareReadStructuresResetOffset(Seq(B(Int.MaxValue, 1)), ReadStructure("1B"))
    compareReadStructuresResetOffset(Seq(M(Int.MaxValue, 1)), ReadStructure("1M"))
    compareReadStructuresResetOffset(Seq(S(Int.MaxValue, 1)), ReadStructure("1S"))
    compareReadStructuresResetOffset(Seq(T(Int.MaxValue, 101)), ReadStructure("101T"))
    compareReadStructuresResetOffset(Seq(B(Int.MaxValue, 5), T(5, 101)), ReadStructure("5B101T"))
    compareReadStructuresResetOffset(Seq(T(Int.MaxValue, 123456789)), ReadStructure("123456789T"))
    compareReadStructuresResetOffset(Seq(T(Int.MaxValue, 10), B(Int.MaxValue, 10), B(Int.MaxValue, 10), S(Int.MaxValue, 10), M(Int.MaxValue, 10)), ReadStructure("10T10B10B10S10M"))
  }

  it should "not be built from invalid structures" in {
    an[Exception] shouldBe thrownBy { ReadStructure("0T") }
    Try { ReadStructure("9R")        }.failed.get.getMessage should include ("[9R]")
    Try { ReadStructure("T")         }.failed.get.getMessage should include ("[T]")
    Try { ReadStructure("23TT")      }.failed.get.getMessage should include ("23T[T]")
    Try { ReadStructure("23T2")      }.failed.get.getMessage should include ("23T[2]")
    Try { ReadStructure("23T2TT23T") }.failed.get.getMessage should include ("23T2T[T]23T")
  }

  it should "collect segments of a single type" in {
    val rs = ReadStructure("10M9T8B7S10M9T8B7S")
    rs.templateSegments         should contain theSameElementsInOrderAs Seq(T(10, 9), T(44, 9))
    rs.molecularBarcodeSegments should contain theSameElementsInOrderAs Seq(M(0, 10), M(34, 10))
    rs.sampleBarcodeSegments    should contain theSameElementsInOrderAs Seq(B(19, 8), B(53, 8))
    rs.skipSegments             should contain theSameElementsInOrderAs Seq(S(27, 7), S(61, 7))
  }

  "ReadStructure.withVariableLastSegment" should "convert the last segment to a variant length segment" in {
    ReadStructure("75T").withVariableLastSegment.toString   shouldBe "+T"
    ReadStructure("5M70T").withVariableLastSegment.toString shouldBe "5M+T"
    ReadStructure("+B").withVariableLastSegment.toString    shouldBe "+B"
    ReadStructure("5B+T").withVariableLastSegment.toString  shouldBe "5B+T"
  }

  "ReadStructure.extract" should "get extract the bases for each segment" in {
    val rs = ReadStructure("2T2B2M2S")
    rs.extract("AACCGGTT").foreach { r =>
        r.kind match {
          case Template         => r.bases shouldBe "AA"
          case SampleBarcode    => r.bases shouldBe "CC"
          case MolecularBarcode => r.bases shouldBe "GG"
          case Skip             => r.bases shouldBe "TT"
        }
    }

    an[Exception] should be thrownBy rs.extract("AAAAAAA")

    // the last segment is truncated
    rs.withVariableLastSegment.extract("AACCGGT").foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"
        case SampleBarcode    => r.bases shouldBe "CC"
        case MolecularBarcode => r.bases shouldBe "GG"
        case Skip             => r.bases shouldBe "T"
      }
    }
    // the last segment is skipped
    rs.withVariableLastSegment.extract("AACCGG").foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"
        case SampleBarcode    => r.bases shouldBe "CC"
        case MolecularBarcode => r.bases shouldBe "GG"
        case Skip             => r.bases shouldBe ""
      }
    }
  }


  "ReadStructure.extract(bases, quals)" should "get extract the bases and qualities for each segment" in {
    val rs = ReadStructure("2T2B2M2S")
    rs.extract("AACCGGTT", "11223344").foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case Skip             => r.bases shouldBe "TT"; r.quals shouldBe "44"
      }
    }
    an[Exception] should be thrownBy rs.extract("AAAAAAA", "AAAAAAA")

    // the last segment is truncated
    rs.withVariableLastSegment.extract("AACCGGT", "1122334").foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case Skip             => r.bases shouldBe "T";  r.quals shouldBe "4"
      }
    }
    // the last segment is skipped
    rs.withVariableLastSegment.extract("AACCGG", "112233").foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case Skip             => r.bases shouldBe ""  ; r.quals shouldBe ""
      }
    }
  }

  "ReadStructure.length" should "return the number of segments" in {
    ReadStructure("1T").length shouldBe 1
    ReadStructure("1B").length shouldBe 1
    ReadStructure("1M").length shouldBe 1
    ReadStructure("1S").length shouldBe 1
    ReadStructure("101T").length shouldBe 1
    ReadStructure("5B101T").length shouldBe 2
    ReadStructure("123456789T").length shouldBe 1
    ReadStructure("10T10B10B10S10M").length shouldBe 5
  }

  "ReadStructure.apply(idx: Int)" should "return the segment at the 0-based index" in {
    ReadStructure("1T")(0) shouldBe T(0, 1)
    ReadStructure("1B")(0) shouldBe B(0, 1)
    ReadStructure("1M")(0) shouldBe M(0, 1)
    ReadStructure("1S")(0) shouldBe S(0, 1)
    ReadStructure("101T")(0) shouldBe T(0, 101)

    ReadStructure("5B101T")(0) shouldBe B(0, 5)
    ReadStructure("5B101T")(1) shouldBe T(5, 101)

    ReadStructure("123456789T")(0) shouldBe T(0, 123456789)

    ReadStructure("10T10B10B10S10M")(0) shouldBe T(0, 10)
    ReadStructure("10T10B10B10S10M")(1) shouldBe B(10, 10)
    ReadStructure("10T10B10B10S10M")(2) shouldBe B(20, 10)
    ReadStructure("10T10B10B10S10M")(3) shouldBe S(30, 10)
    ReadStructure("10T10B10B10S10M")(4) shouldBe M(40, 10)

    an[Exception] should be thrownBy ReadStructure("101T")(1)
  }

  "ReadSegment.apply(offset: Int, length: Int, ch: Char)" should "create a new ReadSegment" in {
    ReadSegment(1, 2, 'T') shouldBe T(1, 2)
    ReadSegment(1, 2, 'M') shouldBe M(1, 2)
    ReadSegment(1, 2, 'S') shouldBe S(1, 2)
    ReadSegment(1, 2, 'B') shouldBe B(1, 2)
    an[Exception] should be thrownBy ReadSegment(1, 2, 'G')
  }

  "ReadSegment.extract(bases)" should "get extract the bases and qualities for a segment" in {
    val molecularBarcodeSegment = ReadStructure("2T2B2M2S").segments(MolecularBarcode).head
    molecularBarcodeSegment.extract("AACCGGTT").bases shouldBe "GG"
  }

  "ReadSegment.extract(bases, quals)" should "get extract the bases and qualities for a segment" in {
    val molecularBarcodeSegment = ReadStructure("2T2B2M2S").segments(MolecularBarcode).head
    molecularBarcodeSegment.extract("AACCGGTT", "11223344") match {
      case SubReadWithQuals(bases, quals, seg) =>  bases shouldBe "GG"; quals shouldBe "33"
    }
  }

  "SegmentType.toString" should "return the String value of the segment code letter" in {
    SegmentType.values.foreach { kind => kind.toString shouldBe String.valueOf(kind.code) }
  }
}
