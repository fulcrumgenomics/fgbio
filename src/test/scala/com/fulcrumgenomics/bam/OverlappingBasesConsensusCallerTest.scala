/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore

class OverlappingBasesConsensusCallerTest extends UnitSpec {

  private def caller(maskDisagreements: Boolean = false,
                     maxQualOnAgreement: Boolean = false): OverlappingBasesConsensusCaller = {
    new OverlappingBasesConsensusCaller(maskDisagreements=maskDisagreements, maxQualOnAgreement=maxQualOnAgreement)
  }

  private def quals(q: Int, rl: Int): String = (33 + q).toChar.toString * rl

  private val minQual: String = quals(q=PhredScore.MinValue, rl=1)
  private val q10: String = quals(q=10, rl=1)
  private val q20: String = quals(q=20, rl=1)
  private val q30: String = quals(q=30, rl=1)

  "OverlappingBasesConsensusCaller.call" should "not change a matching single base overlap" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="A"*10, start1=1, start2=10, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(1, 0, 0)
    r1.basesString shouldBe "A"*10
    r1.qualsString shouldBe q10*9  + q30
    r2.basesString shouldBe "A"*10
    r2.qualsString shouldBe q30 + q20*9
  }

  it should "pick the base with the highest quality (r1) with a single base overlap (disagreement)" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10, start1=1, start2=10, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(1, 1, 0)
    r1.basesString shouldBe "A"*9  + "C"
    r1.qualsString shouldBe q10*10
    r2.basesString shouldBe "C"*10
    r2.qualsString shouldBe q10 + q20*9
  }

  it should "pick the base with the highest quality (r2) with a single base overlap (disagreement)" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10, start1=1, start2=10, quals1=quals(q=20, rl=10), quals2=quals(q=10, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(1, 0, 1)
    r1.basesString shouldBe "A"*10
    r1.qualsString shouldBe q20*9 + q10
    r2.basesString shouldBe "A" + "C"*9
    r2.qualsString shouldBe q10*10
  }

  it should "consensus call fully overlapping reads" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10, start1=1, start2=1, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(10, 10, 0)
    r1.basesString shouldBe "C"*10
    r1.qualsString shouldBe q10*10
    r2.basesString shouldBe "C"*10
    r2.qualsString shouldBe q10*10
  }

  it should "consensus call reads who extend past their mate" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10, start1=5, start2=1, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(6, 6, 0)
    r1.basesString shouldBe "C"*6 + "A"*4
    r1.qualsString shouldBe q10*10
    r2.basesString shouldBe "C"*10
    r2.qualsString shouldBe q20*4 + q10*6
  }

  it should "consensus call read pairs with indels" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10,
      start1=1, start2=1,
      // REF: B-BB-BBBBBBB
      // R1:  A-AAAAAA-AAA
      // R2:  CCCC-C-CCCCC
      // R1': C-CCACAC-CCC (8 bases updated)
      // R2': CCCC-C-CCCCC (0 bases updated)
      cigar1="3M1I3M1D3M",
      cigar2="1M1I3M1D5M",
      quals1=quals(q=10, rl=10),
      quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller().call(r1, r2) shouldBe CorrectionStats(8, 8, 0)
    r1.basesString shouldBe "CCCACACCCC"
    r1.qualsString shouldBe quals(q=10, rl=10)
    r2.basesString shouldBe "CCCCCCCCCC"
    r2.qualsString shouldBe q10 + q20 + q10*4 + q20 + q10*3
  }

  it should "mask the base and quality with a single base overlap (disagreement) when onlyMaskDisagreements=true" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="C"*10, start1=1, start2=10, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller(maskDisagreements=true).call(r1, r2) shouldBe CorrectionStats(1, 0, 0)
    r1.basesString shouldBe "A"*9  + "N"
    r1.qualsString shouldBe q10*9 + minQual
    r2.basesString shouldBe "N" + "C"*9
    r2.qualsString shouldBe minQual + q20*9
  }

  it should "pick the maximum base quality with a single base overlap (agreement) when maxQualOnAgreement=true" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(
      bases1="A"*10, bases2="A"*10, start1=1, start2=10, quals1=quals(q=10, rl=10), quals2=quals(q=20, rl=10))
    r1.matesOverlap.contains(true) shouldBe true
    r2.matesOverlap.contains(true) shouldBe true

    caller(maxQualOnAgreement=true).call(r1, r2) shouldBe CorrectionStats(1, 0, 0)
    r1.basesString shouldBe "A"*10
    r1.qualsString shouldBe q10*9 + q20
    r2.basesString shouldBe "A"* 10
    r2.qualsString shouldBe q20*10
  }
}
