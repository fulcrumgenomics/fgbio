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

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{DovetailPairs, SamBuilder, UnitSpec}
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools.{CigarOperator, TextCigarCodec}
import org.scalatest.OptionValues

import scala.annotation.nowarn

class SamRecordTest extends UnitSpec with OptionValues {
  "SamRecord.isFrPair" should "return false for reads that are not part of a mapped pair" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addFrag(start=100).exists(_.isFrPair) shouldBe false
    builder.addPair(start1=100, start2=100, unmapped2=true).exists(_.isFrPair) shouldBe false
  }

  it should "return false for pairs that are mapped to different chromosomes" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200).foreach { rec =>
      rec.mateRefIndex = rec.refIndex + 1
      rec.isFrPair shouldBe false
    }
  }

  it should "return false for mapped FF, RR and RF pairs" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200, strand1=Plus,  strand2=Plus).exists(_.isFrPair)  shouldBe false
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Minus).exists(_.isFrPair) shouldBe false
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Plus).exists(_.isFrPair)  shouldBe false
  }

  it should "return true on an actual FR pair" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200, strand1=Plus, strand2=Minus).forall(_.isFrPair)  shouldBe true
  }

  "SamRecord.pairOrientation(r1, r2)" should "return None for unmapped or cross-chromosomal pairs" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    val Seq(a, b) = builder.addPair(start1=100, start2=200, unmapped2=true).toSeq
    SamRecord.pairOrientation(a, b) shouldBe None
    SamRecord.pairOrientation(b, a) shouldBe None
    val Seq(c, d) = builder.addPair(start1=100, start2=200).toSeq
    d.refIndex = c.refIndex + 1
    SamRecord.pairOrientation(c, d) shouldBe None
  }

  it should "return FR, RF and TANDEM regardless of argument order" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    Seq(
      (Plus,  Minus, PairOrientation.FR),
      (Minus, Plus,  PairOrientation.RF),
      (Plus,  Plus,  PairOrientation.TANDEM),
      (Minus, Minus, PairOrientation.TANDEM),
    ).foreach { case (s1, s2, expected) =>
      val Seq(r1, r2) = builder.addPair(start1=100, start2=200, strand1=s1, strand2=s2).toSeq
      SamRecord.pairOrientation(r1, r2).value shouldBe expected
      SamRecord.pairOrientation(r2, r1).value shouldBe expected
    }
  }

  it should "return FR when the 5' ends coincide and RF when they are one base apart" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    val Seq(a, b) = builder.addPair(start1=100, start2=91).toSeq // both 5' ends at 100
    SamRecord.pairOrientation(a, b).value shouldBe PairOrientation.FR
    val Seq(c, d) = builder.addPair(start1=101, start2=91).toSeq // forward 5' end at 101, reverse 5' end at 100
    SamRecord.pairOrientation(c, d).value shouldBe PairOrientation.RF
  }

  it should "return FR for a dovetail pair without MC tags, where the per-record orientation disagrees" in {
    val (r1, r2) = DovetailPairs.pair(withMateCigar=false)
    r1.isFrPair shouldBe false // the forward-strand record falls back to TLEN and gets this wrong
    r2.isFrPair shouldBe true
    SamRecord.pairOrientation(r1, r2).value shouldBe PairOrientation.FR
    SamRecord.pairOrientation(r2, r1).value shouldBe PairOrientation.FR
  }

  it should "use the records' own positions rather than stale mate fields" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=200).toSeq
    r2.mateStart = 300 // stale: r1 now appears to start after r2 ends
    r2.pairOrientation shouldBe PairOrientation.RF
    SamRecord.pairOrientation(r1, r2).value shouldBe PairOrientation.FR
  }

  "SamRecord.isFrPair(r1, r2)" should "return true only for paired FR records" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=200).toSeq
    SamRecord.isFrPair(r1, r2) shouldBe true
    SamRecord.isFrPair(r2, r1) shouldBe true
    val Seq(r3, r4) = builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Plus).toSeq
    SamRecord.isFrPair(r3, r4) shouldBe false
    r1.paired = false
    SamRecord.isFrPair(r1, r2) shouldBe false
  }

  it should "return true for a dovetail pair without MC tags regardless of argument order" in {
    val (r1, r2) = DovetailPairs.pair(withMateCigar=false)
    SamRecord.isFrPair(r1, r2) shouldBe true
    SamRecord.isFrPair(r2, r1) shouldBe true
  }

  "SamRecord.cigar" should "always return the latest/correct cigar, when when SAMRecord methods are used" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10, cigar="50M").get

    rec.cigar.toString shouldBe          "50M"
    rec.asSam.getCigarString shouldBe    "50M"
    rec.asSam.getCigar.toString shouldBe "50M"

    // Test with setting via SamRecord.cigar =
    rec.cigar = Cigar("10S40M")
    rec.cigar.toString shouldBe          "10S40M"
    rec.asSam.getCigarString shouldBe    "10S40M"
    rec.asSam.getCigar.toString shouldBe "10S40M"

    // Test with setting via SAMRecord.setCigarString
    rec.asSam.setCigarString("25M2D25M")
    rec.cigar.toString shouldBe          "25M2D25M"
    rec.asSam.getCigarString shouldBe    "25M2D25M"
    rec.asSam.getCigar.toString shouldBe "25M2D25M"

    // Test with setting via SAMRecord.setCigar()
    rec.asSam.setCigar(TextCigarCodec.decode("10H50M"))
    rec.cigar.toString shouldBe          "10H50M"
    rec.asSam.getCigarString shouldBe    "10H50M"
    rec.asSam.getCigar.toString shouldBe "10H50M"
  }

  it should "work with empty cigars" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10, unmapped=true, cigar="50M").get

    rec.cigar = Cigar("50M")
    rec.cigar = Cigar.empty
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null

    rec.cigar = Cigar("50M")
    rec.asSam.setCigar(null)
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null

    rec.cigar = Cigar("50M")
    rec.asSam.setCigarString(null)
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null
  }

  "SamRecord" should "provide easy and safe access to extended attributes" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10).get
    rec("ax") = 10

    // apply
    rec[Int]("ax") shouldBe 10
    Option(rec[Int]("bx")).isEmpty shouldBe true

    // get
    rec.get[Int]("ax").value shouldBe 10
    rec.get[Int]("bx").isEmpty shouldBe true

    // attributes
    rec.attributes.toList should contain theSameElementsAs Seq(("RG", "A"), ("ax", 10))

    // getOrElse
    rec.getOrElse[Int]("ax", 5) shouldBe 10
    rec.getOrElse[Int]("bx", 5) shouldBe 5

    // contains
    rec.contains("ax") shouldBe true
    rec.contains("bx") shouldBe false

    // update
    rec("ax") = 5
    rec[Int]("ax") shouldBe 5
    rec("bx") = 10
    rec[Int]("bx") shouldBe 10
  }

  it should "provide access to the record's primary, secondary, and supplementary status" in {
    val rec = new SamBuilder(readLength=50).addFrag(start=10).get

    rec.primary       = true
    rec.supplementary = false
    rec.primary       shouldBe true
    rec.secondary     shouldBe false
    rec.supplementary shouldBe false

    rec.secondary     = true
    rec.supplementary = false
    rec.primary       shouldBe false
    rec.secondary     shouldBe true
    rec.supplementary shouldBe false

    rec.primary       = true
    rec.supplementary = true
    rec.primary       shouldBe true
    rec.secondary     shouldBe false
    rec.supplementary shouldBe true

    rec.primary       = false
    rec.supplementary = true
    rec.primary       shouldBe false
    rec.secondary     shouldBe true
    rec.supplementary shouldBe true
  }

  "SamRecord.clone()" should "clone the record" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10, cigar="50M", attrs=Map("AB" -> 123)).get
    val clone = rec.clone()

    clone.paired shouldBe rec.paired
    clone.refIndex shouldBe rec.refIndex
    clone.attributes shouldBe rec.attributes

    clone.mapped = false
    clone("AB") = null
    clone("BA") = "Hello"
    clone.bases = ("A" * 50).getBytes

    rec.mapped shouldBe true
    rec[Int]("AB") shouldBe 123
    rec.get[String]("BA") shouldBe None
    rec.basesString should not be clone.basesString
  }

  "SamRecord.mateCigar" should "raise an exception for a fragment without a mate" in {
    val Some(rec) = new SamBuilder(readLength = 100).addFrag(start = 1, cigar = "100M")
    an[IllegalArgumentException] shouldBe thrownBy { rec.mateCigar }
  }

  it should "raise an exception for a paired read if the mate is unmapped" in {
    val Seq(rec1, _): Seq[SamRecord] = new SamBuilder(readLength = 100).addPair(start1 = 1, unmapped2 = true)
    an[IllegalArgumentException] shouldBe thrownBy { rec1.mateCigar }
  }

  it should "return None for a paired read without the 'MC' tag set" in {
    val Seq(rec1, rec2): Seq[SamRecord] = new SamBuilder(readLength = 100).addPair(start1 = 1, start2 = 200)
    rec1.remove("MC")
    rec2.remove("MC")
    rec1.mateCigar shouldBe None
    rec2.mateCigar shouldBe None
  }

  it should "return the mates cigar for a paired read with the 'MC' tag set" in {
    val Seq(rec1, rec2): Seq[SamRecord] = new SamBuilder(readLength = 100).addPair(
      start1 = 1,
      start2 = 200,
      cigar1 = "50M1I49M",
      cigar2 = "49M1I50M"
    )
    rec1.mateCigar.value shouldBe Cigar("49M1I50M")
    rec2.mateCigar.value shouldBe Cigar("50M1I49M")
  }

  "SamRecord.mateUnclippedStart/End" should "return None if no mate cigar set" in {
    val builder = new SamBuilder(readLength=50)
    val recs = builder.addPair(start1=10, start2=50)

    recs.foreach { r =>
      r("MC") = null // Clear mate cigar
      r.mateUnclippedStart.isEmpty shouldBe true
      r.mateUnclippedEnd.isEmpty shouldBe true
    }
  }

  it should "return the mate's unclipped start/end" in {
    val builder = new SamBuilder(readLength=50)
    val Seq(r1, r2) = builder.addPair(start1=10, start2=50, cigar1="5S45M10H", cigar2="10S40M5H")

    r1.start shouldBe 10
    r2.mateStart shouldBe r1.start
    r2.start shouldBe 50
    r1.mateStart shouldBe r2.start

    r1.end shouldBe 54
    r2.mateEnd.value shouldBe r1.end
    r2.end shouldBe 89
    r1.mateEnd.value shouldBe r2.end

    r1.unclippedStart shouldBe 5
    r2.mateUnclippedStart.value shouldBe r1.unclippedStart
    r2.unclippedStart shouldBe 40
    r1.mateUnclippedStart.value shouldBe r2.unclippedStart

    r1.unclippedEnd shouldBe 64
    r2.mateUnclippedEnd.value shouldBe r1.unclippedEnd
    r2.unclippedEnd shouldBe 94
    r1.mateUnclippedEnd.value shouldBe r2.unclippedEnd
  }

  "SamRecord.mateUnSoftClippedStart/End" should "return None if no mate cigar set" in {
    val builder = new SamBuilder(readLength=50)
    val recs = builder.addPair(start1=10, start2=50)

    recs.foreach { r =>
      r("MC") = null // Clear mate cigar
      r.mateUnSoftClippedStart.isEmpty shouldBe true
      r.mateUnSoftClippedEnd.isEmpty shouldBe true
    }
  }

  it should "return the mate's un-soft-clipped start/end" in {
    val builder = new SamBuilder(readLength=51)
    val Seq(r1, r2) = builder.addPair(start1=10, start2=50, cigar1="10H5S45M1S10H", cigar2="10H10S40M1S5H")

    r1.start shouldBe 10
    r2.mateStart shouldBe r1.start
    r2.start shouldBe 50
    r1.mateStart shouldBe r2.start

    r1.end shouldBe 54
    r2.mateEnd.value shouldBe r1.end
    r2.end shouldBe 89
    r1.mateEnd.value shouldBe r2.end

    // SAM R1 aligned bases start is at 10 and then there are 5 soft-clip bases to the left, so 10 - 5 = 5
    r1.unSoftClippedStart shouldBe 5
    r2.mateUnSoftClippedStart.value shouldBe r1.unSoftClippedStart
    // SAM R2 aligned bases start is at 50 and then there are 10 soft-clip bases to the left, so 50 - 10 = 40
    r2.unSoftClippedStart shouldBe 40
    r1.mateUnSoftClippedStart.value shouldBe r2.unSoftClippedStart

    // SAM R1 aligned bases end is at (10 + 45) - 1 = 54 and there is 1 soft-clip base to the right so "end" is 55
    r1.unSoftClippedEnd shouldBe 55
    r2.mateUnSoftClippedEnd.value shouldBe r1.unSoftClippedEnd
    // SAM R2 aligned bases end is at (50 + 40) - 1 = 89 and there is 1 soft-clip base to the right so "end" is 90
    r2.unSoftClippedEnd shouldBe 90
    r1.mateUnSoftClippedEnd.value shouldBe r2.unSoftClippedEnd
  }

  "SamRecord.matesOverlap" should "raise exceptions if any pre-requisite check isn't true" in {
    an[IllegalArgumentException] shouldBe thrownBy { new SamBuilder(readLength = 100).addFrag(start = 1).value.matesOverlap }
    an[IllegalArgumentException] shouldBe thrownBy { new SamBuilder(readLength = 100).addPair(start1 = 1, unmapped2 = true).foreach(_.matesOverlap) }
    an[IllegalArgumentException] shouldBe thrownBy { new SamBuilder(readLength = 100).addPair(start2 = 1, unmapped1 = true).foreach(_.matesOverlap) }
  }

  it should "return false when both reads in a pair are mapped but to different contigs" in {
    new SamBuilder(readLength = 100).addPair(start1 = 1, start2 = 1, contig = 0, contig2 = Some(1))
      .foreach(_.matesOverlap.value shouldBe false)
  }

  it should "return false when both reads in a pair are mapped to the same contig but do not overlap" in {
    new SamBuilder(readLength = 100).addPair(start1 = 1, start2 = 101, contig = 0, contig2 = Some(0))
      .foreach(_.matesOverlap.value shouldBe false)
  }

  it should "return true when reads overlap by one base pair in FR orientation" in {
    new SamBuilder(readLength = 100).addPair(start1 = 1, start2 = 100, contig = 0, contig2 = Some(0))
      .foreach(_.matesOverlap.value shouldBe true)
  }

  it should "return true when reads overlap by one base pair in RF orientation" in {
    new SamBuilder(readLength = 100).addPair(start1 = 100, start2 = 1, contig = 0, contig2 = Some(0))
      .foreach(_.matesOverlap.value shouldBe true)
  }

  it should "return true when reads overlap completely" in {
    new SamBuilder(readLength = 100).addPair(start1 = 1, start2 = 1, contig = 0, contig2 = Some(0))
      .foreach(_.matesOverlap.value shouldBe true)
  }

  it should "return true when mates overlap and we know the start of the mate, and None for when we don't know the end of the mate" in {
    val builder = new SamBuilder(readLength = 100)
    val Seq(rec1, rec2) = builder.addPair(start1 = 10, start2 = 1, contig = 0, contig2 = Some(0))
    rec1.remove("MC")
    rec2.remove("MC")
    rec1.matesOverlap shouldBe None // Mate's start is not enclosed by rec, and mate's end cannot be determined
    rec2.matesOverlap.value shouldBe true // Mate's start is enclosed by rec, regardless of where mate end is
  }

  /** A single test case for [[SamRecord.readPosAtReferencePos()]] and/or [[SamRecord.referencePosAtReadPos()]], for a
    * fragment with the given cigar starting at reference position one.
    *
    * If `returnPriorBaseIfInserted` is `true`, then [[SamRecord.readPosAtReferencePos()]] is not tested.  Similarly,
    * if either `returnPriorBaseIfDeleted` or `returnPriorBaseIfSkipped` is `true`, then
    * [[SamRecord.referencePosAtReadPos()]] is not tested.
    *
    * @param cigar the cigar
    * @param readPos the read position; empty if the reference position does not map to a read base
    * @param refPos the reference position; empty if the read position does not map to a reference base
    * @param returnPriorBaseIfInserted passed to [[SamRecord.referencePosAtReadPos()]]
    * @param returnPriorBaseIfDeleted passed to [[SamRecord.readPosAtReferencePos()]]
    * @param returnPriorBaseIfSkipped passed to [[SamRecord.readPosAtReferencePos()]]
    */
  case class TestCase(cigar: String, readPos: Option[Int], refPos: Option[Int],
                      returnPriorBaseIfInserted: Boolean = false,
                      returnPriorBaseIfDeleted: Boolean = false,
                      returnPriorBaseIfSkipped: Boolean = false) {
    require(readPos.isDefined || refPos.isDefined)
    val doRefPosAtReadPos: Boolean = !returnPriorBaseIfDeleted && !returnPriorBaseIfSkipped
    val doReadPosAtRefPos: Boolean = !returnPriorBaseIfInserted
    require(doRefPosAtReadPos || doReadPosAtRefPos)
  }
  object TestCase {
    def apply(cigar: String, readPos: Int, refPos: Int): TestCase = {
      TestCase(cigar=cigar, readPos=Some(readPos), refPos=Some(refPos))
    }
  }

  /** Test cases for [[SamRecord.readPosAtReferencePos()]] and [[SamRecord.referencePosAtReadPos()]]. */
  private val testCases = Seq(
    // all mapped bases
    TestCase("10M",          10,       10),
    TestCase("10M",           1,        1),
    TestCase("3=1X6=",        4,        4), // sequence match/mismatch operators are aligned bases
    // leading soft-clipping
    TestCase("3S9M",    Some(1),     None), // no reference base, start of soft-clipping
    TestCase("3S9M",    Some(3),     None), // no reference base, end of soft-clipping
    TestCase("3S9M",    Some(1),     None, returnPriorBaseIfInserted=true), // soft-clipping should not return a reference base even with returnPriorBaseIfInserted=true
    TestCase("3S9M",          4,        1), // start of aligned bases
    TestCase("3S9M",          5,        2), // second aligned base
    TestCase("3S9M",         10,        7), // end of aligned bases
    // leading hard-clipping
    TestCase("3H9M",          1,        1), // start of aligned bases
    TestCase("3H3S9M",        4,        1), // start of aligned bases
    TestCase("3H3S9M",  Some(3),     None), // no reference base, end of soft-clipping
    // leading padding
    TestCase("3P9M",          1,        1), // start of aligned bases
    TestCase("3P9M",          9,        9), // end of aligned bases
    // leading skip
    TestCase("3N9M",       None,  Some(1)), // no reference base, start of skip
    TestCase("3N9M",       None,  Some(3)), // no reference base, end of skip
    TestCase("3N9M",       None,  Some(1), returnPriorBaseIfSkipped=true), // no read base precedes the skip
    TestCase("3N9M",          1,        4), // start of aligned bases
    TestCase("3N9M",          9,       12), // end of aligned bases
    // leading deletion
    TestCase("2D9M",       None,  Some(2), returnPriorBaseIfDeleted=true), // no read base precedes the deletion
    // deletions
    TestCase("4M1D6M",        4,        4), // before the deletion
    TestCase("4M1D6M",     None,  Some(5)), // in a deletion
    TestCase("4M1D6M",        5,        6), // after the deletion
    TestCase("4M1D6M",  Some(4),   Some(5), returnPriorBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M1D6M",     None,   Some(5), returnPriorBaseIfSkipped=true), // a deletion is not a skip
    TestCase("4M2D6M",        4,        4), // before the deletion
    TestCase("4M2D6M",     None,  Some(5)), // in a deletion
    TestCase("4M2D6M",     None,  Some(6)), // in a deletion
    TestCase("4M2D6M",  Some(4),   Some(5), returnPriorBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M2D6M",  Some(4),   Some(6), returnPriorBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M2D6M",        5,        7), // after the deletion
    TestCase("20M10D20M",    21,       31), // first base of the deletion
    TestCase("4M2I2D6M", Some(4),  Some(5), returnPriorBaseIfDeleted=true), // returns the last aligned base, not the inserted base
    TestCase("3S2D6M",     None,   Some(1), returnPriorBaseIfDeleted=true), // only soft-clipped bases precede the deletion
    // insertions
    TestCase("4M1I6M",        4,        4), // before the insertion
    TestCase("4M1I6M",  Some(5),     None), // in an insertion
    TestCase("4M1I6M",  Some(5),   Some(4), returnPriorBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M1I6M",        6,        5), // after the insertion
    TestCase("4M2I6M",        4,        4), // before the insertion
    TestCase("4M2I6M",  Some(5),     None), // in an insertion
    TestCase("4M2I6M",  Some(6),     None), // in an insertion
    TestCase("4M2I6M",  Some(5),   Some(4), returnPriorBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M2I6M",  Some(6),   Some(4), returnPriorBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M2I6M",       7,         5), // after the insertion
    TestCase("20M10I20M",    31,       21), // first base of the insertion
    TestCase("4M2D2I6M", Some(5),  Some(4), returnPriorBaseIfInserted=true), // returns the last aligned base, not the deleted base
    TestCase("2I9M",     Some(1),     None, returnPriorBaseIfInserted=true), // no aligned base precedes the insertion
    TestCase("9M2I",     Some(11), Some(9), returnPriorBaseIfInserted=true), // trailing insertion
    // trailing soft-clipping
    TestCase("9M3S",          1,        1), // first aligned base
    TestCase("9M3S",          8,        8), // second to last aligned base
    TestCase("9M3S",          9,        9), // last aligned base
    TestCase("9M3S",   Some(10),     None), // no reference base, first base of soft-clipping
    TestCase("9M3S",   Some(12),     None), // no reference base, last base of soft-clipping
    TestCase("9M3S",   Some(10),     None, returnPriorBaseIfInserted=true), // soft-clipping should not return a reference base even with returnPriorBaseIfInserted=true
    // trailing hard-clipping
    TestCase("9M3H",          1,        1), // start of aligned bases
    TestCase("9M3S3H",        9,        9), // end of aligned bases
    TestCase("9M3S3H", Some(10),     None), // no reference base, first base of soft-clipping
    TestCase("9M3S3H", Some(12),     None), // no reference base, last base of soft-clipping
    // trailing padding
    TestCase("9M3P",          1,        1), // start of aligned bases
    TestCase("9M3P",          9,        9), // end of aligned bases
    // trailing skip
    TestCase("9M3N",       None, Some(10)), // no reference base, start of skip
    TestCase("9M3N",       None, Some(12)), // no reference base, end of skip
    TestCase("9M3N",       None, Some(10), returnPriorBaseIfSkipped=true), // no read base follows the skip
    TestCase("9M3N",          9,        9), // end of aligned bases
    // trailing deletion
    TestCase("9M2D",       None, Some(10), returnPriorBaseIfDeleted=true), // no read base follows the deletion
    // reference skip
    TestCase("5M6N5M", Some(5),  Some(10), returnPriorBaseIfSkipped=true), // returns the previous read base
    TestCase("5M6N5M", None,     Some(10), returnPriorBaseIfDeleted=true), // a skip is not a deletion
    TestCase("4M6N4M", Some(4),  Some(4)),
    TestCase("4M6N4M", None,     Some(5)), // in a skip region
    TestCase("4M6N4M", Some(5),  Some(11)),
    TestCase("5M6N5M", None,     Some(10)), // does not return the read base
  )

  "SamRecord.readPosAtReferencePos" should "throw an exception when the reference position is outside the alignment" in {
    val frag1 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M").value
    (the[IllegalArgumentException] thrownBy frag1.readPosAtReferencePos(pos=0)).getMessage should include ("before the start of the alignment")
    (the[IllegalArgumentException] thrownBy frag1.readPosAtReferencePos(pos=11)).getMessage should include ("after the end of the alignment")

    // request a base in the trailing hard-clipping should throw an exception
    val frag2 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M10H").value
    an[IllegalArgumentException] should be thrownBy frag2.readPosAtReferencePos(pos=11)

    // request a base in the trailing soft-clipping should throw an exception
    val frag3 = new SamBuilder(readLength=20).addFrag(start=1, cigar="10M10S").value
    an[IllegalArgumentException] should be thrownBy frag3.readPosAtReferencePos(pos=11)

    // request a base in the leading soft-clipping should throw an exception
    val frag4 = new SamBuilder(readLength=13).addFrag(start=4, cigar="3S10M").value
    an[IllegalArgumentException] should be thrownBy frag4.readPosAtReferencePos(pos=3)
  }

  it should "throw an exception when the record is unmapped" in {
    val frag = new SamBuilder(readLength=10).addFrag(unmapped=true).value
    an[IllegalArgumentException] should be thrownBy frag.readPosAtReferencePos(pos=1)

    // a record with an alignment that is flagged as unmapped should also throw
    val flagged = new SamBuilder(readLength=10).addFrag(start=100, cigar="10M").value
    flagged.unmapped = true
    the[IllegalArgumentException] thrownBy flagged.readPosAtReferencePos(pos=105) should have message (
      s"requirement failed: Read is not mapped: ${flagged}"
    )
  }

  // now run the test cases
  testCases.filter(_.doReadPosAtRefPos).foreach { testCase =>
    testCase.refPos.foreach { pos =>
      it should s"return read position ${testCase.readPos} for reference position ${testCase.refPos} with $testCase" in {
        val frag    = new SamBuilder(readLength=Cigar(testCase.cigar).lengthOnQuery).addFrag(start=1, cigar=testCase.cigar).value
        val readPos = frag.readPosAtReferencePos(
          pos                      = pos,
          returnPriorBaseIfDeleted = testCase.returnPriorBaseIfDeleted,
          returnPriorBaseIfSkipped = testCase.returnPriorBaseIfSkipped
        )
        withClue(testCase) {
          readPos shouldBe testCase.readPos
          // htsjdk returns the prior base for both deletions and skips, or for neither
          if (testCase.returnPriorBaseIfDeleted == testCase.returnPriorBaseIfSkipped) {
            readPos.getOrElse(0) shouldBe frag.asSam.getReadPositionAtReferencePosition(pos, testCase.returnPriorBaseIfDeleted)
          }
        }
      }
    }
  }

  "SamRecord.referencePosAtReadPos" should "throw an exception when the read position is outside the read" in {
    val frag1 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M").value
    (the[IllegalArgumentException] thrownBy frag1.referencePosAtReadPos(pos=0)).getMessage should include ("less than one")
    (the[IllegalArgumentException] thrownBy frag1.referencePosAtReadPos(pos=11)).getMessage should include ("greater than the number of read bases")

    // request a base in the trailing hard-clipping should throw an exception
    val frag2 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M10H").value
    an[IllegalArgumentException] should be thrownBy frag2.referencePosAtReadPos(pos=11)
  }

  it should "throw an exception when the record is unmapped" in {
    val frag = new SamBuilder(readLength=10).addFrag(unmapped=true).value
    an[IllegalArgumentException] should be thrownBy frag.referencePosAtReadPos(pos=1)

    // a record with an alignment that is flagged as unmapped should also throw
    val flagged = new SamBuilder(readLength=10).addFrag(start=100, cigar="10M").value
    flagged.unmapped = true
    the[IllegalArgumentException] thrownBy flagged.referencePosAtReadPos(pos=5) should have message (
      s"requirement failed: Read is not mapped: ${flagged}"
    )
  }

  it should "return None for every read position when the read is entirely clipped" in {
    Seq("10S", "76S", "10H66S", "76S10H", "10H66S10H").foreach { cigar =>
      withClue(s"cigar=$cigar") {
        val frag = new SamBuilder(readLength=Cigar(cigar).lengthOnQuery).addFrag(start=5, cigar=cigar).value
        Range.inclusive(1, frag.length).foreach { pos =>
          frag.referencePosAtReadPos(pos=pos) shouldBe None
          frag.referencePosAtReadPos(pos=pos, returnPriorBaseIfInserted=true) shouldBe None
        }
      }
    }
  }

  "SamRecord.readPosAtReferencePos and SamRecord.referencePosAtReadPos" should "handle a read that is entirely clipped (no aligned bases)" in {
    @nowarn("cat=deprecation") def readPosAtRefPos(rec: SamRecord, pos: Int, returnLastBaseIfDeleted: Boolean): Int = {
      rec.readPosAtRefPos(pos=pos, returnLastBaseIfDeleted=returnLastBaseIfDeleted)
    }
    @nowarn("cat=deprecation") def refPosAtReadPos(rec: SamRecord, pos: Int): Int = rec.refPosAtReadPos(pos=pos)

    Seq("76S", "10H66S", "76S10H", "10H66S10H").foreach { cigar =>
      withClue(s"cigar=$cigar") {
        val frag = new SamBuilder(readLength=Cigar(cigar).lengthOnQuery).addFrag(start=5, cigar=cigar).value
        // no reference bases are aligned, so every reference position is outside the alignment
        Range.inclusive(frag.start - 1, frag.start + 1).foreach { pos =>
          an[IllegalArgumentException] should be thrownBy frag.readPosAtReferencePos(pos=pos)
          Seq(true, false).foreach { returnLastBaseIfDeleted =>
            readPosAtRefPos(frag, pos, returnLastBaseIfDeleted) shouldBe frag.asSam.getReadPositionAtReferencePosition(pos, returnLastBaseIfDeleted)
          }
        }
        Range.inclusive(0, frag.length + 1).foreach { pos =>
          refPosAtReadPos(frag, pos) shouldBe frag.asSam.getReferencePositionAtReadPosition(pos)
        }
      }
    }
  }

  // now run the test cases
  testCases.filter(_.doRefPosAtReadPos).foreach { testCase =>
    testCase.readPos.foreach { pos: Int =>
      it should s"return ref position ${testCase.refPos} for read position ${testCase.readPos} with $testCase" in {
        val frag   = new SamBuilder(readLength=Cigar(testCase.cigar).lengthOnQuery).addFrag(start=1, cigar=testCase.cigar).value
        val refPos = frag.referencePosAtReadPos(pos=pos, returnPriorBaseIfInserted=testCase.returnPriorBaseIfInserted)
        withClue(testCase) {
          refPos shouldBe testCase.refPos
          if (!testCase.returnPriorBaseIfInserted) {
            refPos.getOrElse(0) shouldBe frag.asSam.getReferencePositionAtReadPosition(pos)
          }
        }
      }
    }
  }

  /** Generates random, valid cigars along with fragments aligned with them, for comparing against htsjdk.  Each cigar
    * has optional leading and trailing clipping, and one or more alignment blocks with an optional gap before each one
    * and after the last one.  A gap is a single `I`, `D`, `N` or `P` element, or two adjacent elements of different
    * types drawn from `I`, `D` and `N`. */
  private def randomFragments(seed: Int, count: Int): Seq[SamRecord] = {
    val random    = new scala.util.Random(seed)
    val alignment = IndexedSeq("M", "=", "X")
    val indels    = IndexedSeq("I", "D", "N")
    def length(max: Int): Int = random.nextInt(max) + 1
    def gap(allowPadding: Boolean): String = random.nextInt(4) match {
      case 0 => ""
      case 1 => s"${length(3)}${(if (allowPadding) indels :+ "P" else indels)(random.nextInt(if (allowPadding) 4 else 3))}"
      case _ =>
        val Seq(first, second) = random.shuffle(indels).take(2)
        s"${length(3)}$first${length(3)}$second"
    }

    Range(0, count).map { _ =>
      val leading  = Seq(
        if (random.nextBoolean()) s"${length(3)}H" else "",
        if (random.nextBoolean()) s"${length(3)}S" else "",
      ).mkString
      val trailing = Seq(
        if (random.nextBoolean()) s"${length(3)}S" else "",
        if (random.nextBoolean()) s"${length(3)}H" else "",
      ).mkString
      val middle = Range(0, length(4)).map { _ => gap(allowPadding=true) + s"${length(4)}${alignment(random.nextInt(alignment.length))}" }.mkString
      val cigar  = Cigar(leading + middle + gap(allowPadding=false) + trailing)
      new SamBuilder(readLength=cigar.lengthOnQuery).addFrag(start=10, cigar=cigar.toString()).value
    }
  }

  /** Returns the operator of the cigar element that contains each read position, indexed by 0-based read offset. */
  private def operatorAtReadOffset(cigar: Cigar): IndexedSeq[CigarOperator] = {
    cigar.elems.filter(_.operator.consumesReadBases()).flatMap(e => IndexedSeq.fill(e.length)(e.operator))
  }

  /** Returns the operator of the cigar element that contains each reference offset from the alignment start. */
  private def operatorAtRefOffset(cigar: Cigar): IndexedSeq[CigarOperator] = {
    cigar.elems.filter(_.operator.consumesReferenceBases()).flatMap(e => IndexedSeq.fill(e.length)(e.operator))
  }

  "SamRecord.readPosAtReferencePos" should "agree with htsjdk on random cigars for every combination of flags" in {
    randomFragments(seed=42, count=2000).foreach { frag =>
      val operators = operatorAtRefOffset(frag.cigar)
      Range.inclusive(frag.start, frag.end).foreach { pos =>
        val aligned = frag.asSam.getReadPositionAtReferencePosition(pos, false)
        val prior   = frag.asSam.getReadPositionAtReferencePosition(pos, true) // htsjdk treats skips as deletions
        for (returnPriorBaseIfDeleted <- Seq(true, false); returnPriorBaseIfSkipped <- Seq(true, false)) {
          val returnPrior = operators(pos - frag.start) match {
            case CigarOperator.D => returnPriorBaseIfDeleted
            case CigarOperator.N => returnPriorBaseIfSkipped
            case _               => false
          }
          val expected = if (aligned != 0) aligned else if (returnPrior) prior else 0
          val actual   = frag.readPosAtReferencePos(pos, returnPriorBaseIfDeleted, returnPriorBaseIfSkipped)
          withClue(s"cigar=${frag.cigar} pos=$pos deleted=$returnPriorBaseIfDeleted skipped=$returnPriorBaseIfSkipped") {
            actual.getOrElse(0) shouldBe expected
          }
        }
      }
    }
  }

  "SamRecord.referencePosAtReadPos" should "agree with htsjdk on random cigars with and without returnPriorBaseIfInserted" in {
    randomFragments(seed=43, count=2000).foreach { frag =>
      val operators = operatorAtReadOffset(frag.cigar)
      val refPos    = Range.inclusive(0, frag.cigar.lengthOnQuery).map(p => if (p == 0) 0 else frag.asSam.getReferencePositionAtReadPosition(p))
      Range.inclusive(1, frag.cigar.lengthOnQuery).foreach { pos =>
        // For an inserted base, the reference position of the last aligned base before it, if any
        val priorIfInserted = {
          if (operators(pos - 1) != CigarOperator.I) 0
          else Range(pos - 1, 0, -1).map(refPos).find(_ != 0).getOrElse(0)
        }
        withClue(s"cigar=${frag.cigar} readPos=$pos") {
          frag.referencePosAtReadPos(pos).getOrElse(0) shouldBe refPos(pos)
          frag.referencePosAtReadPos(pos, returnPriorBaseIfInserted=true).getOrElse(0) shouldBe (if (refPos(pos) != 0) refPos(pos) else priorIfInserted)
        }
      }
    }
  }

  "SamRecord.readPosAtReferencePos and SamRecord.referencePosAtReadPos" should "reflect changes to the cigar and start" in {
    val frag = new SamBuilder(readLength=10).addFrag(start=100, cigar="10M").value
    frag.readPosAtReferencePos(pos=105).value shouldBe 6
    frag.referencePosAtReadPos(pos=6).value shouldBe 105

    frag.start = 102
    frag.readPosAtReferencePos(pos=105).value shouldBe 4
    frag.referencePosAtReadPos(pos=6).value shouldBe 107

    frag.cigar = "3S7M"
    frag.readPosAtReferencePos(pos=105).value shouldBe 7
    frag.referencePosAtReadPos(pos=6).value shouldBe 104

    frag.asSam.setCigarString("2M1I7M") // via the htsjdk API
    frag.readPosAtReferencePos(pos=105).value shouldBe 5
    frag.referencePosAtReadPos(pos=3) shouldBe None
  }

  "SamRecord.readPosAtRefPos and SamRecord.refPosAtReadPos" should "agree with htsjdk on any position, including outside the alignment" in {
    @nowarn("cat=deprecation") def readPosAtRefPos(rec: SamRecord, pos: Int, returnLastBaseIfDeleted: Boolean): Int = {
      rec.readPosAtRefPos(pos=pos, returnLastBaseIfDeleted=returnLastBaseIfDeleted)
    }
    @nowarn("cat=deprecation") def refPosAtReadPos(rec: SamRecord, pos: Int): Int = rec.refPosAtReadPos(pos=pos)

    val unmapped = new SamBuilder(readLength=10).addFrag(unmapped=true).value
    (randomFragments(seed=7, count=500) :+ unmapped).foreach { frag =>
      withClue(s"cigar=${frag.cigar} start=${frag.start}") {
        Range.inclusive(-1, frag.start + frag.cigar.lengthOnTarget + 2).foreach { pos =>
          Seq(true, false).foreach { returnLastBaseIfDeleted =>
            withClue(s"pos=$pos returnLastBaseIfDeleted=$returnLastBaseIfDeleted") {
              readPosAtRefPos(frag, pos, returnLastBaseIfDeleted) shouldBe frag.asSam.getReadPositionAtReferencePosition(pos, returnLastBaseIfDeleted)
            }
          }
        }
        Range.inclusive(-1, frag.cigar.lengthOnQuery + 2).foreach { pos =>
          withClue(s"readPos=$pos") { refPosAtReadPos(frag, pos) shouldBe frag.asSam.getReferencePositionAtReadPosition(pos) }
        }
      }
    }
  }
}
