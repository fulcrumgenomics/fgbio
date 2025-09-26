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
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.TextCigarCodec
import org.scalatest.OptionValues

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

  /** A single test case for [[SamRecord.readPosAtReferencePos()]] and/or [[SamRecord.referencePosAtReadPos()]].
    *
    * If `returnLastBaseIfInserted` is `true`, then [[SamRecord.referencePosAtReadPos()]] should not be tested.  Similarly,
    * if either `returnLastBaseIfDeleted` or `returnLastBaseIfSkipped` is `true`, then [[SamRecord.readPosAtReferencePos()]]
    * should not be tested.
    *
    * @param cigar the cigar
    * @param readPos the read position; empty if the reference position does not map to a read base
    * @param refPos the reference position; empty if the read position does not map to a reference base
    * @param returnLastBaseIfInserted passed to [[SamRecord.readPosAtReferencePos()]]
    * @param returnLastBaseIfDeleted passed to [[SamRecord.referencePosAtReadPos()]]
    * @param returnLastBaseIfSkipped passed to [[SamRecord.referencePosAtReadPos()]]
    */
  case class TestCase(cigar: String, readPos: Option[Int], refPos: Option[Int],
                      returnLastBaseIfInserted: Boolean = false,
                      returnLastBaseIfDeleted: Boolean = false,
                      returnLastBaseIfSkipped: Boolean = false) {
    require(readPos.isDefined || refPos.isDefined)
    val doRefPosAtReadPos: Boolean = !returnLastBaseIfDeleted && !returnLastBaseIfSkipped
    val doReadPosAtRefPos: Boolean = !returnLastBaseIfInserted
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
    // leading soft-clipping
    TestCase("3S9M",    Some(1),     None), // no reference base, start of soft-clipping
    TestCase("3S9M",    Some(3),     None), // no reference base, end of soft-clipping
    TestCase("3S9M",          4,        1), // start of aligned bases
    TestCase("3S9M",         10,        7), // end of aligned bases
    // leading hard-clipping
    TestCase("3H9M",          1,        1), // start of aligned bases
    TestCase("3H3S9M",        4,        1), // start of aligned bases
    TestCase("3H3S9M",  Some(3),     None), // no reference base, end of soft-clipping
    // leading padding
    TestCase("3P9M",          1,        1), // start of aligned bases
    TestCase("3P9M",          9,        9), // end of aligned bases
    // leading skip
    TestCase("3N9M",       None,  Some(1)), // no reference base, start of padding
    TestCase("3N9M",       None,  Some(3)), // no reference base, end of padding
    TestCase("3N9M",          1,        4), // start of aligned bases
    TestCase("3N9M",          9,       12), // end of aligned bases
    // deletions
    TestCase("4M1D6M",        4,        4), // before the deletion
    TestCase("4M1D6M",     None,  Some(5)), // in a deletion
    TestCase("4M1D6M",        5,        6), // after the deletion
    TestCase("4M1D6M",  Some(4),   Some(5), returnLastBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M2D6M",        4,        4), // before the deletion
    TestCase("4M2D6M",     None,  Some(5)), // in a deletion
    TestCase("4M2D6M",     None,  Some(6)), // in a deletion
    TestCase("4M2D6M",  Some(4),   Some(5), returnLastBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M2D6M",  Some(4),   Some(6), returnLastBaseIfDeleted=true), // returns the previous mapped base
    TestCase("4M2D6M",        5,        7), // after the deletion
    TestCase("20M10D20M",    21,       31), // first base of the deletion
    // insertions
    TestCase("4M1I6M",        4,        4), // before the insertion
    TestCase("4M1I6M",  Some(5),     None), // in an insertion
    TestCase("4M1I6M",  Some(5),   Some(4), returnLastBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M1I6M",        6,        5), // after the insertion
    TestCase("4M2I6M",        4,        4), // before the insertion
    TestCase("4M2I6M",  Some(5),     None), // in an insertion
    TestCase("4M2I6M",  Some(6),     None), // in an insertion
    TestCase("4M2I6M",  Some(5),   Some(4), returnLastBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M2I6M",  Some(6),   Some(4), returnLastBaseIfInserted=true), // returns the previous mapped base
    TestCase("4M2I6M",       7,         5), // after the insertion
    TestCase("20M10I20M",    31,       21), // first base of the insertion
    // trailing soft-clipping
    TestCase("9M3S",          1,        1), // first aligned base
    TestCase("9M3S",          9,        9), // last aligned base
    TestCase("9M3S",   Some(10),     None), // no reference base, first base of soft-clipping
    TestCase("9M3S",   Some(12),     None), // no reference base, last base of soft-clipping
    // trailing hard-clipping
    TestCase("9M3H",          1,        1), // start of aligned bases
    TestCase("9M3S3H",        9,        9), // end of aligned bases
    TestCase("9M3S3H", Some(10),     None), /// no reference base, first base of soft-clipping
    TestCase("9M3S3H", Some(12),     None), // no reference base, last base of soft-clipping
    // trailing padding
    TestCase("9M3P",          1,        1), // start of aligned bases
    TestCase("9M3P",          9,        9), // end of aligned bases
    // trailing skip
    TestCase("9M3N",       None, Some(10)), // no reference base, start of padding
    TestCase("9M3N",       None, Some(12)), // no reference base, end of padding
    TestCase("9M3N",          9,        9), // end of aligned bases
    // reference skip
    TestCase("5M6N5M", Some(5),  Some(10), returnLastBaseIfSkipped=true), // returns the previous read base
    TestCase("4M6N4M", Some(4),  Some(4)),
    TestCase("4M6N4M", None,     Some(5)), // returns the previous read base
    TestCase("4M6N4M", Some(5),  Some(11)),
    TestCase("5M6N5M", None,     Some(10)), // does not return the read base
  )

  "SamRecord.readPosAtReferencePos" should "throw an exception when the read position is OOB" in {
    val frag1 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M").value
    an[Exception] should be thrownBy frag1.readPosAtReferencePos(pos=0)
    an[Exception] should be thrownBy frag1.readPosAtReferencePos(pos=11)

    // request a base in the trailing hard-clipping should throw an exception
    val frag2 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M10H").value
    an[Exception] should be thrownBy frag2.readPosAtReferencePos(pos=11)

    // request a base in the trailing soft-clipping should throw an exception
    val frag3 = new SamBuilder(readLength=20).addFrag(start=1, cigar="10M10S").value
    an[Exception] should be thrownBy frag3.readPosAtReferencePos(pos=11)
  }

  // now run the test cases
  testCases.filter(_.doReadPosAtRefPos).foreach { testCase =>
    testCase.refPos.foreach { pos =>
      it should s"return read position ${testCase.readPos} for reference position ${testCase.refPos} with cigar ${testCase.cigar}" in {
        val frag    = new SamBuilder(readLength=Cigar(testCase.cigar).lengthOnQuery).addFrag(start=1, cigar=testCase.cigar).value
        val readPos = frag.readPosAtReferencePos(pos=pos, returnLastBaseIfDeleted=testCase.returnLastBaseIfDeleted, testCase.returnLastBaseIfSkipped)
        withClue(testCase) {
          readPos shouldBe testCase.readPos
          if (!testCase.returnLastBaseIfSkipped) {
            readPos.getOrElse(0) shouldBe frag.asSam.getReadPositionAtReferencePosition(pos, testCase.returnLastBaseIfDeleted)
          }
        }
      }
    }
  }

  "SamRecord.referencePosAtReadPos" should "throw an exception when the read position is OOB" in {
    val frag1 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M").value
    an[Exception] should be thrownBy frag1.referencePosAtReadPos(pos=0)
    an[Exception] should be thrownBy frag1.referencePosAtReadPos(pos=11)

    // request a base in the trailing hard-clipping should throw an exception
    val frag2 = new SamBuilder(readLength=10).addFrag(start=1, cigar="10M10H").value
    an[Exception] should be thrownBy frag2.referencePosAtReadPos(pos=11)
  }

  // now run the test cases
  testCases.filter(_.doRefPosAtReadPos).foreach { testCase =>
    testCase.readPos.foreach { pos: Int =>
      it should s"return ref position ${testCase.refPos} for read position ${testCase.readPos} with cigar ${testCase.cigar}" in {
        val frag   = new SamBuilder(readLength=Cigar(testCase.cigar).lengthOnQuery).addFrag(start=1, cigar=testCase.cigar).value
        val refPos = frag.referencePosAtReadPos(pos=pos, returnLastBaseIfInserted=testCase.returnLastBaseIfInserted)
        withClue(testCase) {
          refPos shouldBe testCase.refPos
          if (!testCase.returnLastBaseIfInserted) {
            refPos.getOrElse(0) shouldBe frag.asSam.getReferencePositionAtReadPosition(pos)
          }
        }
      }
    }
  }

}
