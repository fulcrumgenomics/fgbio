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
import com.fulcrumgenomics.util.SegmentType._
import org.scalatest.OptionValues

import scala.util.Try

class ReadStructureTest extends UnitSpec with OptionValues {

  private def compareReadStructures(actual: ReadStructure, expected: Seq[ReadSegment]) = {
    actual shouldBe expected
    actual.toString shouldBe ReadStructure(expected).toString
  }

  private def T(len: Int) = ReadSegment(length = Some(len), kind = Template)
  private def B(len: Int) = ReadSegment(length = Some(len), kind = SampleBarcode)
  private def M(len: Int) = ReadSegment(length = Some(len), kind = MolecularBarcode)
  private def S(len: Int) = ReadSegment(length = Some(len), kind = Skip)
  private val AnyT = ReadSegment(length = None, kind = Template)
  private val AnyM = ReadSegment(length = None, kind = MolecularBarcode)

  "ReadStructure" should "be built from a string" in {
    compareReadStructures(ReadStructure("1T"), Seq(T(1)))
    compareReadStructures(ReadStructure("1B"), Seq(B(1)))
    compareReadStructures(ReadStructure("1M"), Seq(M(1)))
    compareReadStructures(ReadStructure("1S"), Seq(S(1)))
    compareReadStructures(ReadStructure("101T"), Seq(T(101)))
    compareReadStructures(ReadStructure("5B101T"), Seq(B(5), T(101)))
    compareReadStructures(ReadStructure("123456789T"), Seq(T(123456789)))
    compareReadStructures(ReadStructure("10T10B10B10S10M"), Seq(T(10), B(10), B(10), S(10), M(10)))
  }

  it should "allow whitespace within the read structure and strip it" in {
    ReadStructure("75T 8B 8B 75T").toString shouldBe "75T8B8B75T"
    ReadStructure(" 75T  8B   8B     75T  ").toString shouldBe "75T8B8B75T"
  }

  it should "accept the + segment at any position, at most once" in {
    ReadStructure("5M+T") shouldBe Seq(M(5), AnyT)
    ReadStructure("+M")   shouldBe Seq(AnyM)
    ReadStructure("+M70T") shouldBe Seq(AnyM, T(70))
    ReadStructure("8B+M10T") shouldBe Seq(B(8), AnyM, T(10))
    ReadStructure("10T8B+M10T") shouldBe Seq(T(10), B(8), AnyM, T(10))

    an[Exception] shouldBe thrownBy { ReadStructure("++M") }
    an[Exception] shouldBe thrownBy { ReadStructure("5M++T") }
    an[Exception] shouldBe thrownBy { ReadStructure("5M70+T") }
    an[Exception] shouldBe thrownBy { ReadStructure("+M+T") }
    an[Exception] shouldBe thrownBy { ReadStructure("5M+T+B") }
  }

  it should "round-trip non-terminal + structures through toString" in {
    Seq("8B+M10T", "+B10T", "10T8B+M10T", "5M+T", "+M", "+T").foreach { s =>
      ReadStructure(s).toString shouldBe s
    }
  }

  it should "be built from a sequence of segments" in {
    ReadStructure(Seq(T(1))) shouldBe Seq(T(1))
    ReadStructure(Seq(B(5), T(101))) shouldBe Seq(B(5), T(101))
    ReadStructure(Seq(T(10), B(10), B(10), S(10), M(10))) shouldBe Seq(T(10), B(10), B(10), S(10), M(10))
    ReadStructure(Seq(B(8), AnyM, T(10))) shouldBe Seq(B(8), AnyM, T(10))
  }

  it should "not be built from invalid structures" in {
    an[Exception] shouldBe thrownBy { ReadStructure("0T") }
    Try { ReadStructure("9R")        }.failed.get.getMessage should include ("[9R]")
    Try { ReadStructure("T")         }.failed.get.getMessage should include ("[T]")
    Try { ReadStructure("23TT")      }.failed.get.getMessage should include ("23T[T]")
    Try { ReadStructure("23T2")      }.failed.get.getMessage should include ("23T[2]")
    Try { ReadStructure("23T2TT23T") }.failed.get.getMessage should include ("23T2T[T]23T")
  }

  it should "reject segments with zero or negative length" in {
    an[Exception] shouldBe thrownBy { ReadStructure(Seq(ReadSegment(length = Some(0), kind = Template))) }
    an[Exception] shouldBe thrownBy { ReadStructure(Seq(ReadSegment(length = Some(-1), kind = Template))) }
  }

  it should "reject construction with multiple + segments via the segments constructor" in {
    an[Exception] shouldBe thrownBy { ReadStructure(Seq(AnyM, AnyT)) }
  }

  it should "collect segments of a single type" in {
    val rs = ReadStructure("10M9T8B7S10M9T8B7S")
    rs.templateSegments         should contain theSameElementsInOrderAs Seq(T(9), T(9))
    rs.molecularBarcodeSegments should contain theSameElementsInOrderAs Seq(M(10), M(10))
    rs.sampleBarcodeSegments    should contain theSameElementsInOrderAs Seq(B(8), B(8))
    rs.skipSegments             should contain theSameElementsInOrderAs Seq(S(7), S(7))
  }

  "ReadStructure.withVariableLastSegment" should "convert the last segment to a variant length segment" in {
    ReadStructure("75T").withVariableLastSegment.toString   shouldBe "+T"
    ReadStructure("5M70T").withVariableLastSegment.toString shouldBe "5M+T"
    ReadStructure("+B").withVariableLastSegment.toString    shouldBe "+B"
    ReadStructure("5B+T").withVariableLastSegment.toString  shouldBe "5B+T"
  }

  it should "return itself unchanged when a non-terminal + segment is present" in {
    // A non-terminal + already makes the structure variable-length; the result must not have two + segments.
    ReadStructure("8B+M10T").withVariableLastSegment.toString  shouldBe "8B+M10T"
    ReadStructure("+B10T").withVariableLastSegment.toString    shouldBe "+B10T"
    ReadStructure("5T+M5B").withVariableLastSegment.toString   shouldBe "5T+M5B"
  }

  it should "correctly extract bases after withVariableLastSegment on a non-terminal + structure" in {
    // Simulates the DemuxFastqs.variableReadStructures path: calling withVariableLastSegment on a
    // structure that already has a non-terminal + must leave the structure unchanged so that extract
    // still resolves spans correctly.
    val rs        = ReadStructure("8B+M10T")
    val converted = rs.withVariableLastSegment
    converted.toString shouldBe "8B+M10T"  // must be identical — no extra + added
    val bases = "BBBBBBBBUUUUUUUUUUUUTTTTTTTTTT"
    val out   = converted.extract(bases)
    out should have size 3
    out(0).kind  shouldBe SampleBarcode;    out(0).bases shouldBe "BBBBBBBB"
    out(1).kind  shouldBe MolecularBarcode; out(1).bases shouldBe "UUUUUUUUUUUU"
    out(2).kind  shouldBe Template;         out(2).bases shouldBe "TTTTTTTTTT"
  }

  "ReadStructure.extract" should "extract the bases for each segment" in {
    val rs = ReadStructure("2T2B2M2C2S")
    rs.extract("AACCGGNNTT", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"
        case SampleBarcode    => r.bases shouldBe "CC"
        case MolecularBarcode => r.bases shouldBe "GG"
        case CellBarcode      => r.bases shouldBe "NN"
        case Skip             => r.bases shouldBe "TT"
      }
    }

    an[Exception] should be thrownBy rs.extract("AAAAAAA")

    // trailing variable segment absorbs the short read's remaining base
    rs.withVariableLastSegment.extract("AACCGGNNT", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"
        case SampleBarcode    => r.bases shouldBe "CC"
        case MolecularBarcode => r.bases shouldBe "GG"
        case CellBarcode      => r.bases shouldBe "NN"
        case Skip             => r.bases shouldBe "T"
      }
    }
    // trailing + can be zero-length: a read that ends exactly at the end of the fixed portion
    // produces an empty bases string for the + segment (this pins the "0 or more" semantics for +).
    rs.withVariableLastSegment.extract("AACCGGNN", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"
        case SampleBarcode    => r.bases shouldBe "CC"
        case MolecularBarcode => r.bases shouldBe "GG"
        case CellBarcode      => r.bases shouldBe "NN"
        case Skip             => r.bases shouldBe ""
      }
    }
  }

  it should "omit Skip segments from the result by default" in {
    val rs  = ReadStructure("2T2B2M2C2S")
    val out = rs.extract("AACCGGNNTT")
    out.map(_.segment.kind) should contain noElementsOf Seq(Skip)
    out.map(_.segment.kind) should contain theSameElementsInOrderAs Seq(Template, SampleBarcode, MolecularBarcode, CellBarcode)
    rs.extract("AACCGGNNTT", "1122334455").map(_.segment.kind) should contain noElementsOf Seq(Skip)
  }

  it should "include Skip segments when includeSkips is true" in {
    val rs = ReadStructure("2T2B2M2C2S")
    rs.extract("AACCGGNNTT", includeSkips = true).map(_.segment.kind).last shouldBe Skip
    rs.extract("AACCGGNNTT", "1122334455", includeSkips = true).map(_.segment.kind).last shouldBe Skip
  }

  "ReadStructure.extract(bases, quals)" should "extract the bases and qualities for each segment" in {
    val rs = ReadStructure("2T2B2M2C2S")
    rs.extract("AACCGGNNTT", "1122334455", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case CellBarcode      => r.bases shouldBe "NN"; r.quals shouldBe "44"
        case Skip             => r.bases shouldBe "TT"; r.quals shouldBe "55"
      }
    }
    an[Exception] should be thrownBy rs.extract("AAAAAAA", "AAAAAAA")

    // the last segment is truncated to match the read length
    rs.withVariableLastSegment.extract("AACCGGNNT", "112233445", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case CellBarcode      => r.bases shouldBe "NN"; r.quals shouldBe "44"
        case Skip             => r.bases shouldBe "T";  r.quals shouldBe "5"
      }
    }
    // the last segment is zero-length (pins 0-or-more semantics for +)
    rs.withVariableLastSegment.extract("AACCGGNN", "11223344", includeSkips = true).foreach { r =>
      r.kind match {
        case Template         => r.bases shouldBe "AA"; r.quals shouldBe "11"
        case SampleBarcode    => r.bases shouldBe "CC"; r.quals shouldBe "22"
        case MolecularBarcode => r.bases shouldBe "GG"; r.quals shouldBe "33"
        case CellBarcode      => r.bases shouldBe "NN"; r.quals shouldBe "44"
        case Skip             => r.bases shouldBe ""  ; r.quals shouldBe ""
      }
    }
  }

  it should "require bases and quals to be the same length" in {
    an[IllegalArgumentException] should be thrownBy ReadStructure("2T2B").extract("AACC", "111")
    an[IllegalArgumentException] should be thrownBy ReadStructure("2T2B").extract("AACC", "11111")
    an[IllegalArgumentException] should be thrownBy ReadStructure("+T").extract("AACC", "111")
  }

  it should "reject reads that are longer than a fixed-length structure" in {
    val rs = ReadStructure("2T2B")
    an[Exception] should be thrownBy rs.extract("AACCXX")
    an[Exception] should be thrownBy rs.extract("AACCXX", "111111")
  }

  "ReadStructure.extract with non-terminal +" should "extract bases around a middle + segment" in {
    val rs    = ReadStructure("8B+M10T")
    val bases = "BBBBBBBBUUUUUUUUUUUUTTTTTTTTTT"
    val quals = "!!!!!!!!@@@@@@@@@@@@##########"
    bases.length shouldBe 30
    val out = rs.extract(bases, quals)
    out should have size 3
    out(0).kind shouldBe SampleBarcode;    out(0).bases shouldBe "BBBBBBBB";     out(0).quals shouldBe "!!!!!!!!"
    out(1).kind shouldBe MolecularBarcode; out(1).bases shouldBe "UUUUUUUUUUUU"; out(1).quals shouldBe "@@@@@@@@@@@@"
    out(2).kind shouldBe Template;         out(2).bases shouldBe "TTTTTTTTTT";   out(2).quals shouldBe "##########"
  }

  it should "extract bases around a leading + segment" in {
    val rs    = ReadStructure("+B10T")
    val bases = "BBBBBTTTTTTTTTT"
    val out   = rs.extract(bases)
    out should have size 2
    out(0).kind  shouldBe SampleBarcode; out(0).bases shouldBe "BBBBB"
    out(1).kind  shouldBe Template;      out(1).bases shouldBe "TTTTTTTTTT"
  }

  it should "extract bases with fixed segments before and after the +" in {
    val rs    = ReadStructure("10T8B+M10T")
    val bases = "TTTTTTTTTTBBBBBBBBUUUUUUUUUUUUTTTTTTTTTT"
    bases.length shouldBe 40
    val out   = rs.extract(bases)
    out.map(_.bases) shouldBe Seq("TTTTTTTTTT", "BBBBBBBB", "UUUUUUUUUUUU", "TTTTTTTTTT")
  }

  it should "handle a + segment with multiple fixed segments following it" in {
    val rs    = ReadStructure("8B+M5T5S")
    val bases = "BBBBBBBBUUUUUUUUUUUUTTTTTSSSSS"
    bases.length shouldBe 30
    val out   = rs.extract(bases, includeSkips = true)
    out.map(_.bases) shouldBe Seq("BBBBBBBB", "UUUUUUUUUUUU", "TTTTT", "SSSSS")
  }

  it should "accept reads where the + segment has zero width" in {
    val rs    = ReadStructure("8B+M10T")
    val bases = "BBBBBBBBTTTTTTTTTT"
    bases.length shouldBe 18
    val out   = rs.extract(bases)
    out.map(_.bases) shouldBe Seq("BBBBBBBB", "", "TTTTTTTTTT")
  }

  "ReadStructure.hasFixedLength / fixedLength" should "reflect whether a + segment is present" in {
    val fixed    = ReadStructure("10T8B")
    val variable = ReadStructure("10T+B")
    val middle   = ReadStructure("8B+M10T")
    fixed.hasFixedLength shouldBe true
    fixed.fixedLength shouldBe 18
    variable.hasFixedLength shouldBe false
    middle.hasFixedLength shouldBe false
    an[Exception] should be thrownBy variable.fixedLength
    an[Exception] should be thrownBy middle.fixedLength
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
    ReadStructure("10T10B10B10S10M10C").length shouldBe 6
    ReadStructure("8B+M10T").length shouldBe 3
  }

  "ReadStructure.apply(idx: Int)" should "return the segment at the 0-based index" in {
    ReadStructure("1T")(0) shouldBe T(1)
    ReadStructure("1B")(0) shouldBe B(1)
    ReadStructure("1M")(0) shouldBe M(1)
    ReadStructure("1S")(0) shouldBe S(1)
    ReadStructure("101T")(0) shouldBe T(101)

    ReadStructure("5B101T")(0) shouldBe B(5)
    ReadStructure("5B101T")(1) shouldBe T(101)

    ReadStructure("123456789T")(0) shouldBe T(123456789)

    ReadStructure("10T10B10B10S10M")(0) shouldBe T(10)
    ReadStructure("10T10B10B10S10M")(1) shouldBe B(10)
    ReadStructure("10T10B10B10S10M")(2) shouldBe B(10)
    ReadStructure("10T10B10B10S10M")(3) shouldBe S(10)
    ReadStructure("10T10B10B10S10M")(4) shouldBe M(10)

    an[Exception] should be thrownBy ReadStructure("101T")(1)
  }

  "ReadSegment.apply(length, char)" should "create a new ReadSegment" in {
    ReadSegment(2, 'T') shouldBe T(2)
    ReadSegment(2, 'M') shouldBe M(2)
    ReadSegment(2, 'S') shouldBe S(2)
    ReadSegment(2, 'B') shouldBe B(2)
    an[Exception] should be thrownBy ReadSegment(2, 'G')
  }

  "SegmentType.toString" should "return the String value of the segment code letter" in {
    SegmentType.values.foreach { kind => kind.toString shouldBe String.valueOf(kind.code) }
  }
}
