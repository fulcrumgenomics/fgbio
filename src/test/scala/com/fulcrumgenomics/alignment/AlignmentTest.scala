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

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.{TextCigarCodec, CigarOperator => Op}

class AlignmentTest extends UnitSpec {
  /////////////////////////////////////////////////////////////////////////////
  // Tests for the Cigar class
  /////////////////////////////////////////////////////////////////////////////
  "Cigar.apply(String)" should "parse valid cigars" in {
    Cigar("75M")   shouldBe Cigar(IndexedSeq(CigarElem(Op.M, 75)))
    Cigar("5S70M") shouldBe Cigar(IndexedSeq(CigarElem(Op.S, 5),CigarElem(Op.M, 70)))
    Cigar("5=2X30=3I35=") shouldBe Cigar(IndexedSeq(CigarElem(Op.EQ, 5),CigarElem(Op.X, 2),CigarElem(Op.EQ, 30),CigarElem(Op.I, 3),CigarElem(Op.EQ, 35)))
  }

  it should "throw an exception if given an empty or null string" in {
    an[Exception] shouldBe thrownBy { Cigar(null.asInstanceOf[String]) } // don't do this!
    an[Exception] shouldBe thrownBy { Cigar("") } // don't do this!
  }

  Seq("  ", "M75", "75", "*", "NotACigar", "75Y", "M", "-5M").foreach { cigar =>
    it should s"throw an exception if given an invalid Cigar string '$cigar'" in {
      an[Exception] shouldBe thrownBy { Cigar(cigar) }
    }
  }

  "Cigar.coalesce" should "return this if called on a cigar that doesn't need coalescing" in {
    val cigar = Cigar("75M")
    cigar.coalesce shouldBe cigar
    cigar.coalesce eq cigar shouldBe true
  }

  Seq(("10M10M", "20M"), ("10M10I10M", "10M10I10M"), ("10S10S10S10S10M", "40S10M")).foreach { case (raw, expected) =>
      it should s"coalesce operators in $raw to $expected" in { Cigar(raw).coalesce shouldBe Cigar(expected)}
  }

  "Cigar.truncateToQueryLength" should "not truncate a cigar that is already shorter than the given length" in {
    Cigar("75M").truncateToQueryLength(100).toString() shouldBe "75M"
  }

  Seq(("60M", "50M"), ("10H50M", "10H50M"), ("25M10I25M", "25M10I15M"), ("25M10D25M", "25M10D25M"), ("50M10S","50M")).foreach { case (raw, expected) =>
      it should s"truncate $raw to $expected" in { Cigar(raw).truncateToQueryLength(50).toString() shouldBe expected}
  }

  "Cigar.truncateToTargetLength" should "not truncate a cigar that is already shorter than the given length" in {
    Cigar("75M").truncateToTargetLength(100).toString() shouldBe "75M"
  }

  Seq(("60M", "50M"), ("10H50M", "10H50M"), ("25M10I25M", "25M10I25M"), ("25M10D25M", "25M10D15M"), ("50M10S","50M")).foreach { case (raw, expected) =>
      it should s"truncate $raw to $expected" in { Cigar(raw).truncateToTargetLength(50).toString() shouldBe expected}
  }

  "Cigar.alignedBases" should "return 100 for cigar 100M" in {
    Cigar("100M").alignedBases shouldBe 100
  }

  Seq("100M50S", "50S100M", "25S100M25S", "50H100M", "100M50H", "25S25H100M", "100M25S25H", "10H15S100M15S10H").foreach { cigar =>
    it should s"return 100 bases for cigars with clipping: $cigar" in { Cigar(cigar).alignedBases shouldBe 100 }
  }

  Seq("25M10I25M", "25M10D25M", "10I50M", "50M10I", "10I20M10I30M10I", "10I20M10D30M10I").foreach { cigar =>
    it should s"return 50 bases for cigars with indels: $cigar" in { Cigar(cigar).alignedBases shouldBe 50 }
  }

  "Cigar.clippedBases" should "return 0 if not bases are clipped" in {
    Cigar("100M").clippedBases shouldBe 0
  }

  Seq("100M50S", "50S100M", "25S100M25S", "50H100M", "100M50H", "25S25H100M", "100M25S25H", "10H15S100M15S10H").foreach { cigar =>
      it should s"return 50 bases for cigar $cigar" in { Cigar(cigar).clippedBases shouldBe 50 }
  }

  "Cigar.apply(HtsJdkCigar)" should "translate cigars" in {
    val text = "10H10S30M5I30M5D40M10S10S10H"
    val htsjdkCigar = TextCigarCodec.decode("10H10S30M5I30M5D40M10S10S10H")
    val cigar = Cigar(htsjdkCigar)
    cigar.toString() shouldBe text
  }

  "Cigar.isPrefixOf" should "return true when passed itself" in {
    Seq("10S65M", "75M", "40M2I35M", "35M2D35M").map(Cigar(_)).foreach { cigar => cigar.isPrefixOf(cigar) shouldBe true }
  }

  Seq(("75M","100M"), ("5M1I5M", "5M1I"), ("5M1I5M", "5M1I50M")).foreach { case (shorter, longer) =>
    it should s"return true for '$shorter'.isPrefixOf('$longer')" in { Cigar(shorter).isPrefixOf(Cigar(longer)) }
    it should s"return false for '$longer'.isPrefixOf('$shorter')" in { Cigar(longer).isPrefixOf(Cigar(shorter)) }
  }

  Seq(("3M1I3M","3M1D3M"), ("3M1I3M", "4M1I2M"), ("3S1I3M","3M1I3M")).foreach { case (cigar1, cigar2) =>
    it should s"return false for '$cigar1'.isPrefixOf('$cigar2')" in { Cigar(cigar1).isPrefixOf(Cigar(cigar2)) }
    it should s"return false for '$cigar2'.isPrefixOf('$cigar1')" in { Cigar(cigar2).isPrefixOf(Cigar(cigar1)) }
  }

  "Cigar.reverse" should "reverse the elements in the cigar" in {
    Cigar("75M").reverse.toString shouldBe "75M"
    Cigar("10M2D20M").reverse.toString shouldBe "20M2D10M"
  }

  /////////////////////////////////////////////////////////////////////////////
  // Tests for the Alignment class
  /////////////////////////////////////////////////////////////////////////////
  "Alignment.paddedString" should "produce a very simple string given a very simple alignment" in {
    val expected =
      """
      +AACCGGTT
      +||||||||
      +AACCGGTT
       """.stripMargin('+').trim.linesIterator.toSeq

    val alignment = Alignment(expected.head.replace("-", ""), expected.last.replace("-", ""), 1, 1, Cigar("8M"), 1)
    alignment.paddedString() shouldBe expected
  }

  it should "handle an alignment that has a single mismatch in it" in {
    val expected =
      """
      +AACCGGTT
      +||||||.|
      +AACCGGGT
       """.stripMargin('+').trim.linesIterator.toSeq

    Seq("8M", "6=1X2=").foreach { cigar =>
      val alignment = Alignment(expected.head.replace("-", ""), expected.last.replace("-", ""), 1, 1, Cigar("8M"), 1)
      alignment.paddedString() shouldBe expected
    }
  }

  it should "handle an alignment with some insertions and deletions in it" in {
    val expected =
      """
      +AA--GGGTAAACC-GGGTTT
      +||  ||||||||| || |||
      +AACCGGGTAAACCCGG-TTT
       """.stripMargin('+').trim.linesIterator.toSeq

    val alignment = Alignment(expected.head.replace("-", ""), expected.last.replace("-", ""), 1, 1, Cigar("2=2D9=1D2=1I3="), 1)
    alignment.paddedString() shouldBe expected
  }

  it should "handle a messy alignment with indels and mismatches" in {
    val expected =
      """
      +AA--GGGGGAACC-GGGTTT
      +||  |||..|||| || |||
      +AACCGGGTAAACCCGG-TTT
       """.stripMargin('+').trim.linesIterator.toSeq

    val alignment = Alignment(expected.head.replace("-", ""), expected.last.replace("-", ""), 1, 1, Cigar("2M2D9M1D2M1I3M"), 1)
    alignment.paddedString() shouldBe expected
  }

  Seq("5S10M", "10M5H", "5M50N5M", "50P10M").foreach { cigar =>
    it should s"throw an exception with unsupported operator contained in cigar $cigar" in {
      val cig       = Cigar(cigar)
      val alignment = Alignment("AAAAAAAAAA", "AAAAAAAAAA", 1, 1, cig, 1)
      an [Exception] shouldBe thrownBy { alignment.paddedString() }
    }
  }

  it should "use alternative characters if asked to" in {
    val expected =
      """
      |AA..GGGGGAACC.GGGTTT
      |++--+++##++++-++-+++
      |AACCGGGTAAACCCGG.TTT
       """.stripMargin.trim.linesIterator.toSeq

    val alignment = Alignment(expected.head.replace(".", ""), expected.last.replace(".", ""), 1, 1, Cigar("2M2D9M1D2M1I3M"), 1)
    alignment.paddedString(matchChar='+', mismatchChar='#', gapChar='-', padChar='.') shouldBe expected
  }

  "Alignment.subByTarget" should "yield appropriate sub-alignment when all bases are aligned" in {
    val sub = Alignment("AAAAAAAAAA", "AAAAAAAAAA", 1, 1, Cigar("10M"), 10).subByTarget(5, 6)
    sub.targetStart shouldBe 5
    sub.queryStart  shouldBe 5
    sub.cigar.toString() shouldBe "2M"
  }

  it should "remove cigar elements that are outside of the requested sub-region" in {
    val sub = Alignment("AAAAAAAAAA", "ATAAAAAATA", 1, 1, Cigar("1=1X6=1X1="), 10).subByTarget(5, 6)
    sub.targetStart shouldBe 5
    sub.queryStart  shouldBe 5
    sub.cigar.toString() shouldBe "2="
  }

  it should "work correctly when the alignment start and end are not 1" in {
    val sub = Alignment("AAAAAAAAAA", "ATAAAAAATA", 3, 3, Cigar("6="), 6).subByTarget(5, 6)
    sub.targetStart shouldBe 5
    sub.queryStart  shouldBe 5
    sub.cigar.toString() shouldBe "2="
  }

  it should "work with insertions (which don't count as consuming target)" in {
    val sub = Alignment("AAAAATTAAAAA", "AAAAAAAAAA", 1, 1, Cigar("5=2I5="), 10).subByTarget(3, 8)
    sub.targetStart shouldBe 3
    sub.queryStart  shouldBe 3
    sub.cigar.toString() shouldBe "3=2I3="
  }

  it should "work with deletions (which count as consuming target)" in {
    val sub = Alignment("AAAAAAAAAA", "AAAAATTAAAAA", 1, 1, Cigar("5=2D5="), 10).subByTarget(3, 8)
    sub.targetStart shouldBe 3
    sub.queryStart  shouldBe 3
    sub.cigar.toString() shouldBe "3=2D1="
  }

  it should "drop insertions and deletions that are adjacent to the desired region" in {
    val sub = Alignment("AACCAAAAAAAA", "AAAAAAAATTAA", 1, 1, Cigar("2=2I6=2D2="), 10).subByTarget(3, 8)
    sub.targetStart shouldBe 3
    sub.queryStart  shouldBe 5
    sub.cigar.toString() shouldBe "6="
  }

  it should "fail if the start or end is outside of the alignment" in {
    an[Exception] shouldBe thrownBy { Alignment("AAAAA", "TAAAT", 2, 2, Cigar("3="), 0).subByTarget(1, 3) }
    an[Exception] shouldBe thrownBy { Alignment("AAAAA", "TAAAT", 2, 2, Cigar("3="), 0).subByTarget(2, 5) }
  }

  it should "handle a real world test case" in {
    val query  = "GGCCAGAGTCCACAGATTAACCAGGGGATATGCTAGAAA"
    val target =    "CAGAGGCCACAGATTAACCAGGGGATATGCTAGAAA"
    val ali = Alignment(query, target, 1, 1, Cigar("3I5=1X30="), 0)
    val sub = ali.subByTarget(1, 20)
    sub.queryStart shouldBe 4
    sub.targetStart shouldBe 1
    sub.cigar.toString shouldBe "5=1X14="
  }

  "Alignment.subByQuery" should "yield appropriate sub-alignment when all bases are aligned" in {
    val sub = Alignment("AAAAAAAAAA", "AAAAAAAAAA", 1, 1, Cigar("10M"), 10).subByQuery(5, 6)
    sub.targetStart shouldBe 5
    sub.queryStart  shouldBe 5
    sub.cigar.toString() shouldBe "2M"
  }

  it should "fail if the start or end is outside of the alignment" in {
    an[Exception] shouldBe thrownBy { Alignment("TCGAAAAGGA", "AAAA", 4, 1, Cigar("4="), 0).subByQuery(1, 6) }
    an[Exception] shouldBe thrownBy { Alignment("TCGAAAAGGA", "AAAA", 4, 1, Cigar("4="), 0).subByQuery(5, 9) }
  }

  it should "not remove deletions that are 1bp away from the desired end of the sub-region" in {
    val alignment = Alignment("AAAAAAAAAA", "AAAAATTTTTAAAAA", 1, 1, Cigar("5=5D5="), 0)
    alignment.subByQuery(1,  5).cigar.toString() shouldBe "5="
    alignment.subByQuery(1,  6).cigar.toString() shouldBe "5=5D1="
    alignment.subByQuery(5, 10).cigar.toString() shouldBe "1=5D5="
    alignment.subByQuery(6, 10).cigar.toString() shouldBe "5="

    // Sub-by target sees deltions as consuming, so they should always be included
    alignment.subByTarget(1,  5).cigar.toString() shouldBe "5="
    alignment.subByTarget(1,  6).cigar.toString() shouldBe "5=1D"
    alignment.subByTarget(5, 10).cigar.toString() shouldBe "1=5D"
    alignment.subByTarget(6, 11).cigar.toString() shouldBe "5D1="
    alignment.subByTarget(5, 10).cigar.toString() shouldBe "1=5D"
    alignment.subByTarget(5, 11).cigar.toString() shouldBe "1=5D1="
    alignment.subByTarget(6, 10).cigar.toString() shouldBe "5D"
    alignment.subByTarget(6, 11).cigar.toString() shouldBe "5D1="
  }
}

class CigarStatsTest extends UnitSpec {

  "CigarStats" should "compute various stats about a cigar" in {
    CigarStats("100M")             shouldBe CigarStats(100, 100, 100,  0,  0,  0)
    CigarStats("50M10I50M")        shouldBe CigarStats(110, 100, 100,  0,  0,  0)
    CigarStats("50M10D50M")        shouldBe CigarStats(100, 110, 100,  0,  0,  0)
    CigarStats("10H100M")          shouldBe CigarStats(100, 100, 100, 10,  0,  0)
    CigarStats("100M10H")          shouldBe CigarStats(100, 100, 100, 10,  0,  0)
    CigarStats("10S100M")          shouldBe CigarStats(110, 100, 100, 10, 10,  0)
    CigarStats("100M10S")          shouldBe CigarStats(110, 100, 100, 10,  0, 10)
    CigarStats("10H10S100M")       shouldBe CigarStats(110, 100, 100, 20, 10,  0)
    CigarStats("100M10S10H")       shouldBe CigarStats(110, 100, 100, 20,  0, 10)
    CigarStats("1H2S3M4D5I6M7S8H") shouldBe CigarStats(23,   13,   9, 18,  2,  7)
  }

  it should "fail on an invalid cigar" in {
    an[Exception] should be thrownBy CigarStats("100M10S100M")
  }
}
