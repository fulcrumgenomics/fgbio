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

import com.fulcrumgenomics.alignment.Mode.{Global, Glocal, Local}
import com.fulcrumgenomics.commons.util.NumericCounter
import com.fulcrumgenomics.testing.UnitSpec

class AlignerTest extends UnitSpec {
  /** Upper-cases and remove display-related characters from a string. */
  def s(str: String): String = str.filterNot(ch => ch == '-').toUpperCase

  /** Asserts a couple of things about an alignment that should be constant for global alignments. */
  def assertValidGlobalAlignment(alignment: Alignment): Unit = {
    alignment.queryStart  shouldBe 1
    alignment.queryEnd    shouldBe alignment.query.length
    alignment.targetStart shouldBe 1
    alignment.targetEnd   shouldBe alignment.target.length
  }

  /** Asserts a couple of things about an alignment that should be constant for glocal alignments. */
  def assertValidGlocalAlignment(alignment: Alignment): Unit = {
    alignment.queryStart  shouldBe 1
    alignment.queryEnd    shouldBe alignment.query.length
    alignment.targetStart should be >= 1
    alignment.targetEnd   should be <= alignment.target.length
  }

  /** Asserts a couple of things about an alignment that should be constant for local alignments. */
  def assertValidLocalAlignment(alignment: Alignment): Unit = {
    alignment.queryStart  should be >= 1
    alignment.queryEnd    should be <= alignment.query.length
    alignment.targetStart should be >= 1
    alignment.targetEnd   should be <= alignment.target.length
    alignment.score       should be >= 0
  }

  "Aligner.align(Global)" should "align two identical sequences with all matches" in {
    val result = Aligner(1, -1, -3, -1).align(s("ACGTAACC"), s("ACGTAACC"))
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "8="
    result.score shouldBe 8
  }

  it should "align two sequences with a single mismatch in them" in {
    val q = s("AACCGGTT")
    val t = s("AACCGtTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "5=1X2="
    result.score shouldBe 6
  }

  it should "align two sequences with a single small deletion in the query sequence" in {
    val q = s("AACC-GTT")
    val t = s("AACCGGTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "4=1D3="
    result.score shouldBe 7 - 4
  }

  it should "align two sequences with a single small insertion in the query sequence" in {
    val q = s("AACCGGGTT")
    val t = s("AACC-GGTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "4=1I4="
    result.score shouldBe 8 - 4
  }

  it should "align two sequences with compensating insertions and deletions" in {
    val q = s("AAACGCGCGCGCG-TT")
    val t = s("-AACGCGCGCGCGTTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I12=1D2="
    result.score shouldBe 14 - 4 - 4
  }

  it should "align two sequences with a leading insertion" in {
    val q = s("ATTTTTTTTTTT")
    val t = s( "TTTTTTTTTTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11="
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a trailing insertion" in {
    val q = s("TTTTTTTTTTTA")
    val t = s("TTTTTTTTTTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "11=1I"
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a leading deletion" in {
    val q = s( "TTTTTTTTTTT")
    val t = s("ATTTTTTTTTTT")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1D11="
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a trailing deletion" in {
    val q = s("TTTTTTTTTTT")
    val t = s("TTTTTTTTTTTA")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "11=1D"
    result.score shouldBe 11 - 4
  }

  it should "align two sequences preferring a 2bp insertion and mismatch vs two small insertions" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: a leading and trailing 1b insertion ("1I11=1I") has score: 11 - 4 - 4 = 3, where as a leading 2bp
    // insertion with the last base a mismatch ("2I10=1X") has score: 10 - 4 - 1 - 1 = 4.
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 4 - 1 - 1
  }

  it should "align two sequences preferring two small insertions vs a 2bp insertion and mismatch" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: a leading 2bp insertion with the last base a mismatch ("2I10=1X") has score: 10 - 4 - 1 - 2 = 3, whereas a
    // leading and trailing 1b insertion ("1I11=1I") has score: 11 - 4 - 4 = 3
    val result = Aligner(1, -3, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11=1I"
    result.score shouldBe 11 - 4 - 4
  }

  it should "align two sequences preferring a mismatch and a 2bp insertion vs two small deletions when they have the same score" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: both have the same score, but we must consistently choose one over the other.
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 4 - 1 - 1
  }

  it should "align a homopolymer insertion and left-justify the insertion" in {
    val q = s("GTTTTTTTTTTA")
    val t = s("GTTTTTTTTTA")
    // NB: both have the same score, but we must consistently choose one over the other.
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1=1I10="
    result.score shouldBe 11 - 4
  }

  it should "align a simple triplet insertion and left-justify the insertion" in {
    val q = s("GACGACGACGACGA")
    val t = s("GACGACGACGA")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3I11="
    result.score shouldBe 11 - 4 - 1 - 1
  }

  it should "align a simple triplet insertion with leading matches and left-justify the insertion" in {
    val q = s("TTTGACGACGACGACGA")
    val t = s("TTTGACGACGACGA")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3=3I11="
    result.score shouldBe 14 - 4 - 1 - 1
  }

  it should "align a simple triplet deletion with leading matches and left-justify the deletion" in {
    val q = s("TTTGACGACGACGA")
    val t = s("TTTGACGACGACGACGA")
    val result = Aligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3=3D11="
    result.score shouldBe 14 - 4 - 1 - 1
  }

  it should "prefer a mismatch over an insertion and deletion when the mismatch score is less than than two times the gap open + gap extend scores" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = Aligner(1, -3, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2=1X3="
    result.score shouldBe 5 - 3
  }

  it should "prefer a mismatch over an insertion and deletion when the alignments have the same score" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = Aligner(1, -4, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2=1X3="
    result.score shouldBe 5 - 4
  }

  it should "prefer an insertion and deletion over a mismatch when the mismatch score is greater than two times the gap open + gap extend scores" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = Aligner(1, -5, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I2=1D3="
    result.score shouldBe 5 - 2 - 2
  }

  it should "prefer one insertion when the gap open penalty is much larger than the gap extend penalty" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    val result = Aligner(1, -5, -100, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 100 - 1 - 1 - 5
  }

  it should "prefer two small insertions when the gap open penalty is much less than the gap extend penalty" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    val result = Aligner(1, -5, -1, -100).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11=1I"
    result.score shouldBe 11 - 100 - 1 - 100 - 1
  }

  // Glocal alignment tests
  "Aligner.align(Glocal)" should "align two identical sequences with all matches" in {
    val result = Aligner(1, -1, -3, -1, mode=Glocal).align(s("ACGTAACC"), s("ACGTAACC"))
    assertValidGlocalAlignment(result)
    result.cigar.toString() shouldBe "8="
    result.score shouldBe 8
  }

  it should "correctly align an identical subsequence" in {
    val result = Aligner(1, -1, -3, -1, mode=Glocal).align(s("CCGG"), s("AACCGGTT"))
    assertValidGlocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 3
    result.cigar.toString() shouldBe "4="
    result.score shouldBe 4
  }

  it should "correctly align a sub-sequence with a mismatch" in {
    val result = Aligner(1, -1, -3, -1, mode=Glocal).align(s("-------CGCGTCGTATACGTCGTT"), s("AAGATATCGCGTCGTATACGTCGTA"))
    assertValidGlocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 8
    result.cigar.toString() shouldBe "17=1X"
    result.score shouldBe 16
  }

  it should "correctly align a sub-sequence with a gap" in {
    val result = Aligner(1, -1, -3, -1, mode=Glocal).align(s("CGCGCGCG"), s("AACGCGACGCGTT"))
    assertValidGlocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 3
    result.cigar.toString() shouldBe "4=1D4="
    result.score shouldBe 4
  }

  it should "align a query that is longer than the target, creating an insertion" in {
    val result = Aligner(1, -1, -3, -1, mode=Glocal).align("AAAAGGGGTTTT", "AAAATTTT")
    assertValidGlocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 1
    result.cigar.toString shouldBe "4=4I4="
  }

  it should "align a query with leading and trailing gaps" in {
    val q = s("-------------------GGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG---------------------------")
    val t = s("AGGGCTATAGACTGCTAGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAATGAGCTATTAGTCATGACGCTTTT")
    val result = Aligner(1, -3, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "19D54=27D"
    result.score shouldBe 54 - (1*19 + 3) - (1*27 + 3)
  }

  "Aligner.align(Local)" should "align two identical sequences with all matches" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("ACGTAACC"), s("ACGTAACC"))
    assertValidLocalAlignment(result)
    result.cigar.toString() shouldBe "8="
    result.score shouldBe 8
  }

  it should "correctly align an identical subsequence (query in target)" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("CCGG"), s("AACCGGTT"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 3
    result.queryEnd shouldBe 4
    result.targetEnd shouldBe 6
    result.cigar.toString() shouldBe "4="
    result.score shouldBe 4
  }

  it should "correctly align an identical subsequence (target in query)" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("AACCGGTT"), s("CCGG"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 3
    result.targetStart shouldBe 1
    result.cigar.toString() shouldBe "4="
    result.score shouldBe 4
  }

  it should "correctly align a sub-sequence with a leading mismatch" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("AGCGTCGTATACGTCGTA-------"), s("CGCGTCGTATACGTCGTAAAGATAT"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 2
    result.targetStart shouldBe 2
    result.queryEnd shouldBe 18
    result.targetEnd shouldBe 18
    result.cigar.toString() shouldBe "17="
    result.score shouldBe 17
  }


  it should "correctly align a sub-sequence with a trailing mismatch" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("-------CGCGTCGTATACGTCGTT"), s("AAGATATCGCGTCGTATACGTCGTA"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 8
    result.queryEnd shouldBe 17
    result.targetEnd shouldBe 24
    result.cigar.toString() shouldBe "17="
    result.score shouldBe 17
  }

  it should "correctly align a sub-sequence with a gap in the query" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("CCGCGCGCGC"), s("AACCGCGACGCGCTT"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 3
    result.queryEnd shouldBe 10
    result.targetEnd shouldBe 13
    result.cigar.toString() shouldBe "5=1D5="
    result.score shouldBe 6
  }

  it should "correctly align a sub-sequence with a gap in the target" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("AACCGCGACGCGCTT"), s("CCGCGCGCGC"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 3
    result.targetStart shouldBe 1
    result.queryEnd shouldBe 13
    result.targetEnd shouldBe 10
    result.cigar.toString() shouldBe "5=1I5="
    result.score shouldBe 6
  }

  it should "prefer a match over an indel" in {
    val result = Aligner(1, -1, -3, -1, mode=Local).align(s("CGCGCGCG"), s("AACGCGACGCGTT"))
    assertValidLocalAlignment(result)
    result.queryStart shouldBe 1
    result.targetStart shouldBe 3
    result.queryEnd shouldBe 4
    result.targetEnd shouldBe 6
    result.cigar.toString() shouldBe "4="
    result.score shouldBe 4
  }

  it should "generate a cigar string that gives matches for Us and IUPAC codes" in {
    val query  = "AGCGCUGACGUCGUUGACnrg"
    val target = "AGCaCTGACGTCGTTGACGGG"
    val result = Aligner(5, -3, -2, -4, mode=Global).align(query, target)
    result.cigar.toString() shouldBe "3=1X17="

    val Seq(q, aln, t) = result.paddedString()
    aln shouldBe "|||.|||||||||||||||||"
  }

  /** Timing test - change "ignore" to "it" to enable. */
  ignore should "perform lots of glocal alignments" in {
    val count = 25000
    val aligner = Aligner(5, -4, -2, -5, Mode.Glocal)
    val query  =                                       "AcATCTTTCGCATcGACTGAC-TTG".filterNot(_ == '-').toUpperCase.getBytes
    val target = "GCCAGGGACCGTTTCAGACAGATATTTGCCTGGTGGATAGATCTT-CGCATGGACTGACTTTGACGATACTGCTAATTTTTTTTATAGCCTTTGCCTTGTT".filterNot(_ == '-').getBytes
    val counter = new NumericCounter[Double]()

    Range(0, 10).foreach { _ =>
      val startTime = System.currentTimeMillis()
      Range(0, count).foreach { _ => val aln = aligner.align(query, target) }
      val endTime = System.currentTimeMillis()
      val total = (endTime - startTime) / 1000.0
      System.out.println(s"Total time: $total.  Average time: ${total / count}")
      counter.count(total)
    }

    System.out.println(s"Run median=${counter.median()}, mean=${counter.mean()}")
  }
}
