/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Sequences.OffsetAndLength

/**
  * Tests for the Sequences utility object.
  */
class SequencesTest extends UnitSpec {
  "Sequences.countMismatches" should "return 0 when comparing a string with itself" in {
    val seq = "ACGTAGCTGACTGCA"
    Sequences.countMismatches(seq, seq) shouldBe 0
  }

  it should "throw an exception when comparing sequences of different lengths" in {
    an[IllegalArgumentException] should be thrownBy Sequences.countMismatches("ACGT", "ACGTGG")
  }

  it should "correct count mismatches" in {
    Sequences.countMismatches("ACGACTATATGCAT", "acgactatatgcat") shouldBe 0
    Sequences.countMismatches("ACGACTATATGCAT", "acgactaNatgcat") shouldBe 1
    Sequences.countMismatches("ACGACTATATGCAT", "acgactataGTcat") shouldBe 2
    Sequences.countMismatches("ACGACTATATGCAT", "NNNNNNNNNNNNNN") shouldBe 14
  }

  "Sequences.longestHomopolymer" should "fail when passed a null String" in {
    an[Throwable] shouldBe thrownBy { Sequences.longestHomopolymer(null) }
  }

  it should "return zero when passed a zero length string" in {
    Sequences.longestHomopolymer("") shouldBe OffsetAndLength(0, 0)
  }

  it should "correctly find the single longest homopolymer" in {
    Sequences.longestHomopolymer("A") shouldBe OffsetAndLength(0,1)
    Sequences.longestHomopolymer("AAAA") shouldBe OffsetAndLength(0,4)
    Sequences.longestHomopolymer("ATTTGCGAT") shouldBe OffsetAndLength(1,3)
    Sequences.longestHomopolymer("GGGGCGTGC") shouldBe OffsetAndLength(0,4)
    Sequences.longestHomopolymer("ACGTACGTCCCCC") shouldBe OffsetAndLength(8,5)
    Sequences.longestHomopolymer("AAAACCCGGGGAAAAA") shouldBe OffsetAndLength(11,5)
    Sequences.longestHomopolymer("ATTCCCGG") shouldBe OffsetAndLength(3,3)
  }

  it should "report the earliest homopolymer of the longest length when there are multiple" in {
    Sequences.longestHomopolymer("AAAAACCCCCGGGGGTTTTT") shouldBe OffsetAndLength(0,5)
    Sequences.longestHomopolymer("ACGTAACCGGTTAAACCCGGGTTT") shouldBe OffsetAndLength(12,3)
  }

  "Sequences.maxDinuc" should "correctly count the longest dinucleotide in sequences" in {
    Sequences.longestDinuc("A")        shouldBe OffsetAndLength(0, 0)
    Sequences.longestDinuc("ACGT")     shouldBe OffsetAndLength(0, 2)
    Sequences.longestDinuc("ACACGTT")  shouldBe OffsetAndLength(0, 4)
    Sequences.longestDinuc("AGTGTGT")  shouldBe OffsetAndLength(1, 6)
    Sequences.longestDinuc("GGGGGGGG") shouldBe OffsetAndLength(0, 8)
    Sequences.longestDinuc("AGCGtagCGCGCgcGCTCTCTatCGCGCA") shouldBe OffsetAndLength(6, 10)
  }

  "Sequences.complement" should "return the complement of sequences" in {
    Sequences.complement("AAA") shouldBe "TTT"
    Sequences.complement("ACAC") shouldBe "TGTG"
    Sequences.complement("") shouldBe ""
    Sequences.complement("AACCGGTGTG") shouldBe "TTGGCCACAC"
  }

  "Sequences.revcomp" should "return the reverse complement of sequences" in {
    Sequences.revcomp("AAA")        shouldBe "TTT"
    Sequences.revcomp("ACAC")       shouldBe "GTGT"
    Sequences.revcomp("")           shouldBe ""
    Sequences.revcomp("AACCGGTGTG") shouldBe "CACACCGGTT"
    Sequences.revcomp("acacNNNN")   shouldBe "NNNNgtgt"
    Sequences.revcomp("NRG")        shouldBe "CYN"
  }

  "Sequences.compatible" should "do the right thing for pairs of bases" in {
    val bases: Seq[Byte] = Seq('A', 'C', 'G', 'T', 'a', 'c', 'g', 't')
    for (b1 <- bases; b2 <- bases) {
      Sequences.compatible(b1, b2) shouldBe b1.toChar.toUpper == b2.toChar.toUpper
      Sequences.compatible(b1, b2) shouldBe b1.toChar.toUpper == b2.toChar.toUpper
      Sequences.compatible(b1, b2) shouldBe b1.toChar.toUpper == b2.toChar.toUpper
      Sequences.compatible(b1, b2) shouldBe b1.toChar.toUpper == b2.toChar.toUpper
    }

    // Check T/U compatability
    Sequences.compatible('T', 'U') shouldBe true
    Sequences.compatible('u', 'T') shouldBe true

    // Spot check some abiguity codes
    Sequences.compatible('H', 'C') shouldBe true
    Sequences.compatible('R', 'A') shouldBe true
    Sequences.compatible('S', 'G') shouldBe true
    bases.foreach(b => Sequences.compatible(b, 'N'.toByte) shouldBe true)

    // Check that incompatible combinations yield false
    Sequences.compatible('A', 'C') shouldBe false
    Sequences.compatible('M', 'K') shouldBe false
  }

  "Sequences.iupacCode()" should "return the same base given when given only one base" in {
    Seq('A', 'C', 'G', 'T', 'R', 'H', 'N').map(_.toByte).foreach { b =>
      Sequences.iupacCode(Seq(b)) shouldBe b
    }
  }

  it should "give the correct combination of regular bases" in {
    Sequences.iupacCode(Seq('A', 'C').map(_.toByte)) shouldBe 'M'.toByte
    Sequences.iupacCode(Seq('T', 'G').map(_.toByte)) shouldBe 'K'.toByte
    Sequences.iupacCode(Seq('A', 'C', 'T').map(_.toByte)) shouldBe 'H'.toByte
  }

  it should "handle duplicate entries" in {
    Sequences.iupacCode(Seq('A', 'C', 'A', 'C', 'A').map(_.toByte)) shouldBe 'M'.toByte
  }

  it should "handle IUPAC codes in the input" in {
    Sequences.iupacCode(Seq('A', 'N').map(_.toByte)) shouldBe 'N'.toByte
    Sequences.iupacCode(Seq('R', 'S').map(_.toByte)) shouldBe 'V'.toByte
  }

  "Sequences.homopolymers" should "return zeros for empty string" in {
    val zeroBase = ""
    Sequences.homopolymers(zeroBase, 1) shouldBe 0
  }

  it should "work with lower case sequences" in {
    val lowercase = "ggggccccttttatatatgatag"
    Sequences.homopolymers(lowercase, 4) shouldBe 3
  }

  it should "work with mixed case sequences" in {
    val mixedCase = "gGGgcCCcttTTAtatatGAtag" // 10/23
    Sequences.homopolymers(mixedCase, 4) shouldBe 3
  }

  it should "calculate number of homopolymers correctly" in {
    val smallKmers = "AAATTTAAATTTAAAAATTTTTAATT"
    Sequences.homopolymers(smallKmers, 1) shouldBe 8
    Sequences.homopolymers(smallKmers, 3) shouldBe 6
    Sequences.homopolymers(smallKmers, 5) shouldBe 2
  }

  "Sequences.gcContent" should "produce zeros for empty string" in {
    val zeroBase = "" // 0
    Sequences.gcContent(zeroBase) shouldBe 0
  }

  it should "work with lower case sequences" in {
    val lowercase = "ggggccccttttatatatgatag" // 10/23
    Sequences.gcContent(lowercase) shouldBe (10.0 / 23)
  }

  it should "work with mixed case sequences" in {
    val mixedCase = "gGGgcCCcttTTAtatatGAtag" // 10/23
    Sequences.gcContent(mixedCase) shouldBe (10.0 / 23)
  }

  it should "calculate Gc content correctly" in {
    val noGc = "tatatatatatatatatata" // 0
    val tenBaseAllGc = "CGGCGGGGGG"
    val tenBaseHalfGc = "CGCGCATATT"
    Sequences.gcContent(noGc) shouldBe 0
    Sequences.gcContent(tenBaseAllGc) shouldBe 1.0
    Sequences.gcContent(tenBaseHalfGc) shouldBe 0.5
  }
}
