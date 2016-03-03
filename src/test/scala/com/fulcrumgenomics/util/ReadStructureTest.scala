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
import org.scalatest.OptionValues

import scala.util.Try

class ReadStructureTest extends UnitSpec with OptionValues {

  private def compareReadStructures(actual: ReadStructure, expected: Seq[ReadSegment]): Unit = {
    // make sure the segments match
    actual shouldBe expected
    // make sure the string representations are the same
    actual.toString shouldBe new ReadStructure(expected).toString
  }

  private def compareReadStructuresResetOffset(actual: Seq[ReadSegment], expected: ReadStructure): Unit = {
    val actualReadStructure = new ReadStructure(actual, resetOffsets=true)
    compareReadStructures(actualReadStructure, expected)
  }

  "ReadStructure" should "be built from a string" in {
    compareReadStructures(ReadStructure("1T"), Seq(Template(0, 1)))
    compareReadStructures(ReadStructure("1B"), Seq(SampleBarcode(0, 1)))
    compareReadStructures(ReadStructure("1M"), Seq(MolecularBarcode(0, 1)))
    compareReadStructures(ReadStructure("1S"), Seq(Skip(0, 1)))
    compareReadStructures(ReadStructure("101T"), Seq(Template(0, 101)))
    compareReadStructures(ReadStructure("5B101T"), Seq(SampleBarcode(0, 5), Template(5, 101)))
    compareReadStructures(ReadStructure("123456789T"), Seq(Template(0, 123456789)))
    compareReadStructures(ReadStructure("10T10B10B10S10M"), Seq(Template(0, 10), SampleBarcode(10, 10), SampleBarcode(20, 10), Skip(30, 10), MolecularBarcode(40, 10)))
  }
  
  it should "be built from segments while resetting their offset" in {
    compareReadStructuresResetOffset(Seq(Template(Int.MaxValue, 1)), ReadStructure("1T"))
    compareReadStructuresResetOffset(Seq(SampleBarcode(Int.MaxValue, 1)), ReadStructure("1B"))
    compareReadStructuresResetOffset(Seq(MolecularBarcode(Int.MaxValue, 1)), ReadStructure("1M"))
    compareReadStructuresResetOffset(Seq(Skip(Int.MaxValue, 1)), ReadStructure("1S"))
    compareReadStructuresResetOffset(Seq(Template(Int.MaxValue, 101)), ReadStructure("101T"))
    compareReadStructuresResetOffset(Seq(SampleBarcode(Int.MaxValue, 5), Template(5, 101)), ReadStructure("5B101T"))
    compareReadStructuresResetOffset(Seq(Template(Int.MaxValue, 123456789)), ReadStructure("123456789T"))
    compareReadStructuresResetOffset(Seq(Template(Int.MaxValue, 10), SampleBarcode(Int.MaxValue, 10), SampleBarcode(Int.MaxValue, 10), Skip(Int.MaxValue, 10), MolecularBarcode(Int.MaxValue, 10)), ReadStructure("10T10B10B10S10M"))
  }

  it should "not be built from invalid structures" in {
    Try { ReadStructure("0T")        }.failed.get.getMessage should include ("[0]T")
    Try { ReadStructure("9R")        }.failed.get.getMessage should include ("9[R]")
    Try { ReadStructure("T")         }.failed.get.getMessage should include ("[T]")
    Try { ReadStructure("23TT")      }.failed.get.getMessage should include ("23T[T]")
    Try { ReadStructure("23T2")      }.failed.get.getMessage should include ("23T[2]")
    Try { ReadStructure("23T2TT23T") }.failed.get.getMessage should include ("23T2T[T]23T")
  }

  it should "collect segments of a single type" in {
    val rs = ReadStructure("10M9T8B7S10M9T8B7S")
    rs.template should contain theSameElementsInOrderAs Seq(Template(10, 9), Template(44, 9))
    rs.molecularBarcode should contain theSameElementsInOrderAs Seq(MolecularBarcode(0, 10), MolecularBarcode(34, 10))
    rs.sampleBarcode should contain theSameElementsInOrderAs Seq(SampleBarcode(19, 8), SampleBarcode(53, 8))
    rs.skip should contain theSameElementsInOrderAs Seq(Skip(27, 7), Skip(61, 7))
  }


  "ReadStructure.totalBases" should "return the total bases described by the read structure" in {
    ReadStructure("1T").totalBases shouldBe 1
    ReadStructure("1B").totalBases shouldBe 1
    ReadStructure("1M").totalBases shouldBe 1
    ReadStructure("1S").totalBases shouldBe 1
    ReadStructure("101T").totalBases shouldBe 101
    ReadStructure("123456789T").totalBases shouldBe 123456789
    ReadStructure("10T10B10B10S10M").totalBases shouldBe 50
  }

  "ReadStructure.structureRead" should "get extract the bases for each segment" in {
    val rs = ReadStructure("2T2B2M2S")
    rs.structureRead("AACCGGTT").foreach { case r =>
        r.segment match {
          case t: Template         => r.bases shouldBe "AA"
          case s: SampleBarcode    => r.bases shouldBe "CC"
          case m: MolecularBarcode => r.bases shouldBe "GG"
          case x: Skip             => r.bases shouldBe "TT"
        }
    }
    an[IllegalStateException] should be thrownBy rs.structureRead("AAAAAAA")
    // the last segment is truncated
    rs.structureRead("AACCGGT", strict=false).foreach { case r =>
      r.segment match {
        case t: Template         => r.bases shouldBe "AA"
        case s: SampleBarcode    => r.bases shouldBe "CC"
        case m: MolecularBarcode => r.bases shouldBe "GG"
        case x: Skip             => r.bases shouldBe "T"
      }
    }
    // the last segment is skipped
    rs.structureRead("AACCGG", strict=false).foreach { case r =>
      r.segment match {
        case t: Template         => r.bases shouldBe "AA"
        case s: SampleBarcode    => r.bases shouldBe "CC"
        case m: MolecularBarcode => r.bases shouldBe "GG"
        case x: Skip             => throw new IllegalStateException("Skip should not be found")
      }
    }
  }


  "ReadStructure.structureReadWithQualities" should "get extract the bases and qualities for each segment" in {
    val rs = ReadStructure("2T2B2M2S")
    rs.structureReadWithQualities("AACCGGTT", "11223344").foreach { case r =>
      r.segment match {
        case t: Template         => r.bases shouldBe "AA"; r.quals.value shouldBe "11"
        case s: SampleBarcode    => r.bases shouldBe "CC"; r.quals.value shouldBe "22"
        case m: MolecularBarcode => r.bases shouldBe "GG"; r.quals.value shouldBe "33"
        case x: Skip             => r.bases shouldBe "TT"; r.quals.value shouldBe "44"
      }
    }
    an[IllegalStateException] should be thrownBy rs.structureReadWithQualities("AAAAAAA", "AAAAAAA")
    // the last segment is truncated
    rs.structureReadWithQualities("AACCGGT", "1122334", strict=false).foreach { case r =>
      r.segment match {
        case t: Template         => r.bases shouldBe "AA"; r.quals.value shouldBe "11"
        case s: SampleBarcode    => r.bases shouldBe "CC"; r.quals.value shouldBe "22"
        case m: MolecularBarcode => r.bases shouldBe "GG"; r.quals.value shouldBe "33"
        case x: Skip             => r.bases shouldBe "T";  r.quals.value shouldBe "4"
      }
    }
    // the last segment is skipped
    rs.structureReadWithQualities("AACCGG", "112233", strict=false).foreach { case r =>
      r.segment match {
        case t: Template         => r.bases shouldBe "AA"; r.quals.value shouldBe "11"
        case s: SampleBarcode    => r.bases shouldBe "CC"; r.quals.value shouldBe "22"
        case m: MolecularBarcode => r.bases shouldBe "GG"; r.quals.value shouldBe "33"
        case x: Skip             => throw new IllegalStateException("Skip should not be found")
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
    ReadStructure("1T")(0) shouldBe Template(0, 1)
    ReadStructure("1B")(0) shouldBe SampleBarcode(0, 1)
    ReadStructure("1M")(0) shouldBe MolecularBarcode(0, 1)
    ReadStructure("1S")(0) shouldBe Skip(0, 1)
    ReadStructure("101T")(0) shouldBe Template(0, 101)

    ReadStructure("5B101T")(0) shouldBe SampleBarcode(0, 5)
    ReadStructure("5B101T")(1) shouldBe Template(5, 101)

    ReadStructure("123456789T")(0) shouldBe Template(0, 123456789)

    ReadStructure("10T10B10B10S10M")(0) shouldBe Template(0, 10)
    ReadStructure("10T10B10B10S10M")(1) shouldBe SampleBarcode(10, 10)
    ReadStructure("10T10B10B10S10M")(2) shouldBe SampleBarcode(20, 10)
    ReadStructure("10T10B10B10S10M")(3) shouldBe Skip(30, 10)
    ReadStructure("10T10B10B10S10M")(4) shouldBe MolecularBarcode(40, 10)

    an[Exception] should be thrownBy ReadStructure("101T")(1)
  }

  "ReadStructure.sampleBarcodeStructure and friends" should "return a new read structure with the correct offsets" in {
    val readStructure = ReadStructure("1B1M1S1T1B1M1S1T")

    readStructure.templateStructure(resetOffsets=false) should contain theSameElementsInOrderAs readStructure.template
    readStructure.sampleBarcodeStructure(resetOffsets=false) should contain theSameElementsInOrderAs readStructure.sampleBarcode
    readStructure.molecularBarcodeStructure(resetOffsets=false) should contain theSameElementsInOrderAs readStructure.molecularBarcode
    readStructure.skipStructure(resetOffsets=false) should contain theSameElementsInOrderAs readStructure.skip

    readStructure.templateStructure(resetOffsets=true) should contain theSameElementsInOrderAs ReadStructure("1T1T")
    readStructure.sampleBarcodeStructure(resetOffsets=true) should contain theSameElementsInOrderAs ReadStructure("1B1B")
    readStructure.molecularBarcodeStructure(resetOffsets=true) should contain theSameElementsInOrderAs ReadStructure("1M1M")
    readStructure.skipStructure(resetOffsets=true) should contain theSameElementsInOrderAs ReadStructure("1S1S")
  }

  "ReadSegment.apply(offset: Int, length: Int, s: String)" should "create a new ReadSegment" in {
    ReadSegment(1, 2, "T") shouldBe Template(1, 2)
    ReadSegment(1, 2, "M") shouldBe MolecularBarcode(1, 2)
    ReadSegment(1, 2, "S") shouldBe Skip(1, 2)
    ReadSegment(1, 2, "B") shouldBe SampleBarcode(1, 2)
    ReadSegment(1, 2, "TGARBAGE") shouldBe Template(1, 2)
    an[Exception] should be thrownBy ReadSegment(1, 2, "GARBAGE")
  }

  "ReadSegment.apply(segment: ReadSegment, length: Int)" should "create a new ReadSegment" in {
    ReadSegment(ReadSegment(1, 2, "T"), 3) shouldBe Template(1, 3)
    ReadSegment(ReadSegment(1, 2, "M"), 3) shouldBe MolecularBarcode(1, 3)
    ReadSegment(ReadSegment(1, 2, "S"), 3) shouldBe Skip(1, 3)
    ReadSegment(ReadSegment(1, 2, "B"), 3) shouldBe SampleBarcode(1, 3)
  }

  "ReadSegment.structureRead" should "get extract the bases and qualities for a segment" in {
    val molecularBarcodeSegment = ReadStructure("2T2B2M2S").filter { case x: MolecularBarcode => true; case _ => false }.head
    molecularBarcodeSegment.extractBases("AACCGGTT") shouldBe "GG"
  }

  "ReadSegment.structureReadWithQualities" should "get extract the bases and qualities for a segment" in {
    val molecularBarcodeSegment = ReadStructure("2T2B2M2S").filter { case x: MolecularBarcode => true; case _ => false }.head
    molecularBarcodeSegment.extractBasesAndQualities("AACCGGTT", "11223344") match {
      case (bases, quals) =>  bases shouldBe "GG"; quals shouldBe "33"
    }
  }
}
