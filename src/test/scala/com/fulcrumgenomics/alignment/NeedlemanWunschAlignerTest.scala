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

class NeedlemanWunschAlignerTest extends UnitSpec {
  /** Upper-cases and remove display-related characters from a string. */
  def s(str: String): String = str.filterNot(ch => ch == '-').toUpperCase

  /** Asserts a couple of things about an alignment that should be constant for global alignments. */
  def assertValidGlobalAlignment(alignment: Alignment): Unit = {
    alignment.queryStart shouldBe 1
    alignment.queryEnd   shouldBe  alignment.query.length

    alignment.targetStart shouldBe 1
    alignment.targetEnd   shouldBe alignment.target.length
  }

  "NeedlemanWunschAligner.align" should "align two identical sequences with all matches" in {
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(s("ACGTAACC"), s("ACGTAACC"))
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "8="
    result.score shouldBe 8
  }

  it should "align two sequences with a single mismatch in them" in {
    val q = s("AACCGGTT")
    val t = s("AACCGtTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "5=1X2="
    result.score shouldBe 6
  }

  it should "align two sequences with a single small deletion in the query sequence" in {
    val q = s("AACC-GTT")
    val t = s("AACCGGTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "4=1D3="
    result.score shouldBe 7 - 4
  }

  it should "align two sequences with a single small insertion in the query sequence" in {
    val q = s("AACCGGGTT")
    val t = s("AACC-GGTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "4=1I4="
    result.score shouldBe 8 - 4
  }

  it should "align two sequences with compensating insertions and deletions" in {
    val q = s("AAACGCGCGCGCG-TT")
    val t = s("-AACGCGCGCGCGTTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I12=1D2="
    result.score shouldBe 14 - 4 - 4
  }

  it should "align two sequences with a leading insertion" in {
    val q = s("ATTTTTTTTTTT")
    val t = s( "TTTTTTTTTTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11="
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a trailing insertion" in {
    val q = s("TTTTTTTTTTTA")
    val t = s("TTTTTTTTTTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "11=1I"
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a leading deletion" in {
    val q = s( "TTTTTTTTTTT")
    val t = s("ATTTTTTTTTTT")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1D11="
    result.score shouldBe 11 - 4
  }

  it should "align two sequences with a trailing deletion" in {
    val q = s("TTTTTTTTTTT")
    val t = s("TTTTTTTTTTTA")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "11=1D"
    result.score shouldBe 11 - 4
  }

  it should "align two sequences preferring a 2bp insertion and mismatch vs two small insertions" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: a leading and trailing 1b insertion ("1I11=1I") has score: 11 - 4 - 4 = 3, where as a leading 2bp
    // insertion with the last base a mismatch ("2I10=1X") has score: 10 - 4 - 1 - 1 = 4.
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 4 - 1 - 1
  }

  it should "align two sequences preferring two small insertions vs a 2bp insertion and mismatch" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: a leading 2bp insertion with the last base a mismatch ("2I10=1X") has score: 10 - 4 - 1 - 2 = 3, whereas a
    // leading and trailing 1b insertion ("1I11=1I") has score: 11 - 4 - 4 = 3
    val result = NeedlemanWunschAligner(1, -3, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11=1I"
    result.score shouldBe 11 - 4 - 4
  }

  it should "align two sequences preferring a mismatch and a 2bp insertion vs two small deletions when they have the same score" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    // NB: both have the same score, but we must consistently choose one over the other.
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 4 - 1 - 1
  }

  it should "align a homopolymer insertion and left-justify the insertion" in {
    val q = s("GTTTTTTTTTTA")
    val t = s("GTTTTTTTTTA")
    // NB: both have the same score, but we must consistently choose one over the other.
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1=1I10="
    result.score shouldBe 11 - 4
  }

  it should "align a simple triplet insertion and left-justify the insertion" in {
    val q = s("GACGACGACGACGA")
    val t = s("GACGACGACGA")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3I11="
    result.score shouldBe 11 - 4 - 1 - 1
  }

  it should "align a simple triplet insertion with leading matches and left-justify the insertion" in {
    val q = s("TTTGACGACGACGACGA")
    val t = s("TTTGACGACGACGA")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3=3I11="
    result.score shouldBe 14 - 4 - 1 - 1
  }

  it should "align a simple triplet deletion with leading matches and left-justify the deletion" in {
    val q = s("TTTGACGACGACGA")
    val t = s("TTTGACGACGACGACGA")
    val result = NeedlemanWunschAligner(1, -1, -3, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "3=3D11="
    result.score shouldBe 14 - 4 - 1 - 1
  }

  it should "prefer a mismatch over an insertion and deletion when the mismatch score is less than than two times the gap open + gap extend scores" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = NeedlemanWunschAligner(1, -3, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2=1X3="
    result.score shouldBe 5 - 3
  }

  it should "prefer a mismatch over an insertion and deletion when the alignments have the same score" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = NeedlemanWunschAligner(1, -4, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2=1X3="
    result.score shouldBe 5 - 4
  }

  it should "prefer an insertion and deletion over a mismatch when the mismatch score is greater than two times the gap open + gap extend scores" in {
    val q = s("AAACCC")
    val t = s("AACCCC")
    // NB: could be either "1I2=1D3=" or "2=1X3="
    val result = NeedlemanWunschAligner(1, -5, -1, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I2=1D3="
    result.score shouldBe 5 - 2 - 2
  }

  it should "prefer one insertion when the gap open penalty is much larger than the gap extend penalty" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    val result = NeedlemanWunschAligner(1, -5, -100, -1).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "2I10=1X"
    result.score shouldBe 10 - 100 - 1 - 1 - 5
  }

  it should "prefer two small insertions when the gap open penalty is much less than the gap extend penalty" in {
    val q = s("ATTTTTTTTTTTA")
    val t = s( "TTTTTTTTTTT")
    val result = NeedlemanWunschAligner(1, -5, -1, -100).align(q, t)
    assertValidGlobalAlignment(result)
    result.cigar.toString() shouldBe "1I11=1I"
    result.score shouldBe 11 - 100 - 1 - 100 - 1
  }
}
