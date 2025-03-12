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

  "SamRecord.refBasesFromMd" should "preserve ambiguous (N) reference bases" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAA", cigar="3M", start=1, attrs=Map("MD" -> "1N1")).value
    new String(rec.refBasesFromMd) shouldBe "ANA"
  }

  it should "omit soft-clipped and inserted bases, and include deleted bases" in {
    val builder = new SamBuilder(readLength=9)
    val rec     = builder.addFrag(bases="TTACGGGAC", cigar="2S2M1I1M1D3M", start=1, attrs=Map("MD" -> "2C0^T3")).value
    new String(rec.refBasesFromMd) shouldBe "ACCTGAC"
  }

  it should "return Ns for reference skips" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="ACGT", cigar="2M3N2M", start=1, attrs=Map("MD" -> "4")).value
    new String(rec.refBasesFromMd) shouldBe "ACNNNGT"
    new String(rec.refBasesFromMd(includeSkippedBases=false)) shouldBe "ACGT"
  }

  it should "accept any letter for mismatched and deleted reference bases" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="ACTA", cigar="2M1I1D1M", start=1, attrs=Map("MD" -> "1R0^Y1")).value
    new String(rec.refBasesFromMd) shouldBe "ARYA"
  }

  Seq("0\u01411", "0\u03b11", "0\u04111").foreach { md =>
    it should s"reject non-ASCII mismatch letters in MD $md" in {
      val rec = new SamBuilder(readLength=3).addFrag(bases="AAT", cigar="1M1I1M", start=1, attrs=Map("MD" -> md)).value
      an[IllegalArgumentException] should be thrownBy rec.refBasesFromMd
    }
  }

  Seq("1^\u01411", "1^A\u01411", "1^A0^\u01411").foreach { md =>
    it should s"reject non-ASCII deletion letters in MD $md" in {
      val cigar = if (md == "1^\u01411") "1M1D1M" else "1M1D1D1M"
      val rec = new SamBuilder(readLength=2).addFrag(bases="AT", cigar=cigar, start=1, attrs=Map("MD" -> md)).value
      an[IllegalArgumentException] should be thrownBy rec.refBasesFromMd
    }
  }

  it should "retain every upper- and lower-case ASCII MD letter and zero-padded match counts" in {
    val rec = new SamBuilder(readLength=2).addFrag(bases="AT", cigar="1M1D1M", start=1).value
    (('A' to 'Z') ++ ('a' to 'z')).foreach { base =>
      rec("MD") = s"00${base}0^${base}01"
      new String(rec.refBasesFromMd) shouldBe s"${base}${base}T"
    }
  }

  it should "reject non-ASCII MD digits" in {
    val bases = "A" * 1586
    val rec = new SamBuilder(readLength=bases.length).addFrag(bases=bases, cigar="1586M", start=1, attrs=Map("MD" -> "\u0662")).value
    an[IllegalArgumentException] should be thrownBy rec.refBasesFromMd
    rec("MD") = "1586"
    new String(rec.refBasesFromMd) shouldBe bases
  }

  it should "accept a trailing zero after a mismatch" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAC", cigar="3M", start=1, attrs=Map("MD" -> "2A0")).value
    new String(rec.refBasesFromMd) shouldBe "AAA"
  }

  it should "accept adjacent deletions written as one run or as separate runs (as htsjdk and samtools write them)" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="ACGT", cigar="2M1D1D2M", start=1, attrs=Map("MD" -> "2^T0^G2")).value
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    // the MD tag that htsjdk writes for this alignment
    val expected = rec.asSam.deepCopy()
    htsjdk.samtools.util.SequenceUtil.calculateMdAndNmTags(expected, "ACTGGT".getBytes, true, false)
    expected.getStringAttribute("MD") shouldBe "2^T0^G2"
    rec("MD") = "2^TG2"
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    // a single deletion of two bases, given as one run or as separate runs
    rec.cigar = "2M2D2M"
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    rec("MD") = "2^T0^G2"
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    // only a zero may separate the runs within a deletion
    rec("MD") = "2^T1^G1"
    intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage should include ("expected a deleted base")
  }

  it should "accept zero-padded counts between MD deletion runs for equivalent CIGARs" in {
    val cases = Seq(
      ("1M2D2M",   "1^A", "^C2", "AACAT"),
      ("1M1D1D2M", "1^A", "^C2", "AACAT"),
      ("2D3M",     "0^A", "^C3", "ACAAT"),
      ("3M2D",     "3^A", "^C0", "AATAC")
    )
    cases.foreach { case (cigar, prefix, suffix, expected) =>
      Seq("0", "00", "000").foreach { count =>
        val md = prefix + count + suffix
        withClue(s"cigar=$cigar MD=$md: ") {
          val rec = new SamBuilder(readLength=3).addFrag(
            bases="AAT", cigar=cigar, start=1, attrs=Map("MD" -> md)
          ).value
          new String(rec.refBasesFromMd) shouldBe expected
        }
      }
    }
  }

  it should "reject nonzero or malformed MD counts within a deletion" in {
    Seq("1", "01", "10", "4294967296", "00C", "00^^").foreach { separator =>
      val md = s"1^A$separator^C2"
      withClue(s"MD=$md: ") {
        val rec = new SamBuilder(readLength=3).addFrag(
          bases="AAT", cigar="1M2D2M", start=1, attrs=Map("MD" -> md)
        ).value
        an[IllegalArgumentException] should be thrownBy rec.refBasesFromMd
      }
    }
  }

  it should "accept deletions separated by an insertion" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="ACG", cigar="1M1D1I1D1M", start=1, attrs=Map("MD" -> "1^T0^A1")).value
    new String(rec.refBasesFromMd) shouldBe "ATAG"
  }

  for (operator <- Seq("P", "I"); md <- Seq("1^AC1", "1^A0^C1", "1^A00^C1")) {
    it should s"keep an MD deletion run across $operator with MD=$md" in {
      val bases = if (operator == "I") "TAG" else "TG"
      val rec = new SamBuilder(readLength=bases.length).addFrag(bases=bases, cigar=s"1M1D1${operator}1D1M",
        start=1, attrs=Map("MD" -> md)).value
      new String(rec.refBasesFromMd) shouldBe "TACG"
    }
  }

  it should "reject an '=' read base that the MD tag reports as matching" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="A=A", cigar="2M1D1M", start=1, attrs=Map("MD" -> "2^C1")).value
    intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage should include ("for an '=' read base")
    // an '=' read base that is a mismatch in the MD tag is not needed
    rec("MD") = "1G0^C1"
    new String(rec.refBasesFromMd) shouldBe "AGCA"
  }

  it should "reject an MD tag that does not match the CIGAR" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAA", cigar="3M", start=1, attrs=Map("MD" -> "4")).value
    def failure(md: String): String = { rec("MD") = md; intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage }
    failure("4") should include ("the MD tag is longer than the alignment")
    failure("2") should include ("expected a match or mismatch")
    failure("1^A2") should include ("expected a match or mismatch")
    failure("1AC1") should include ("expected a number between mismatched or deleted bases")
    rec.cigar = "1M1D2M"
    failure("1A2") should include ("expected a deletion")
    failure("1^2") should include ("expected a deleted base")
    failure("0A^C2") should include ("expected a number between mismatched or deleted bases")
    failure("1^CA1") should include ("expected a number between mismatched or deleted bases")
  }

  it should "require the MD tag" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAA", cigar="3M", start=1).value
    intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage should include ("without the MD tag")
  }

  it should "reject an MD tag that is not a string" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAA", cigar="3M", start=1, attrs=Map("MD" -> 3)).value
    intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage should include ("MD tag is not a string")
  }

  it should "reject a match count in the MD tag that overflows" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAA", cigar="3M", start=1, attrs=Map("MD" -> "4294967299")).value
    intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage should include ("too large")
    rec("MD") = "0003"
    new String(rec.refBasesFromMd) shouldBe "AAA"
  }

  it should "reject an MD tag that contradicts the = and X operators of the CIGAR" in {
    val builder = new SamBuilder(readLength=3)
    val rec     = builder.addFrag(bases="AAT", cigar="1X1I1M", start=1, attrs=Map("MD" -> "2")).value
    def failure(md: String): String = { rec("MD") = md; intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage }
    failure("2") should include ("expected a mismatch for an X operator")
    failure("0A1") should include ("expected a mismatch for an X operator")
    failure("0a1") should include ("expected a mismatch for an X operator")
    rec("MD") = "0C1"
    new String(rec.refBasesFromMd) shouldBe "CT"
    rec.cigar = "1=1I1M"
    failure("0C1") should include ("expected a match for an = operator")
    failure("0N1") should include ("expected a match for an = operator")
    rec("MD") = "2"
    new String(rec.refBasesFromMd) shouldBe "AT"
    // an M operator accepts either
    rec.cigar = "1M1I1M"
    new String(rec.refBasesFromMd) shouldBe "AT"
    rec("MD") = "0C1"
    new String(rec.refBasesFromMd) shouldBe "CT"
    // an N may be counted as either a match or a mismatch
    rec.bases = "NAT"
    rec.cigar = "1X1I1M"
    rec("MD") = "2"
    new String(rec.refBasesFromMd) shouldBe "NT"
    rec.cigar = "1=1I1M"
    rec("MD") = "0N1"
    new String(rec.refBasesFromMd) shouldBe "NT"
    failure("0A1") should include ("expected a match for an = operator")
    rec.cigar = "1X1I1M"
    rec("MD") = "0N1"
    new String(rec.refBasesFromMd) shouldBe "NT"
  }

  it should "require a number before and after each mismatched base and run of deleted bases" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="ACGT", cigar="4M", start=1, attrs=Map("MD" -> "4")).value
    def failure(md: String): String = { rec("MD") = md; intercept[IllegalArgumentException] { rec.refBasesFromMd }.getMessage }
    failure("C3") should include ("expected a number at the start of the MD tag")
    failure("3C") should include ("expected a number after the last mismatched or deleted base")
    rec.cigar = "2M1D1D2M"
    failure("2^T^G2") should include ("expected a number between mismatched or deleted bases")
    rec("MD") = "2^T0^G2"
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    rec("MD") = "2^TG2"
    new String(rec.refBasesFromMd) shouldBe "ACTGGT"
    rec.cigar = "4M1D"
    failure("4^T") should include ("expected a number after the last mismatched or deleted base")
    rec("MD") = "4^T0"
    new String(rec.refBasesFromMd) shouldBe "ACGTT"
    rec.cigar = "1D4M"
    failure("^A4") should include ("expected a number at the start of the MD tag")
    rec("MD") = "0^A4"
    new String(rec.refBasesFromMd) shouldBe "AACGT"
  }
}
