/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.identifyprimers

import com.fulcrumgenomics.testing.UnitSpec

class PrimerPairMatchTypeTest extends UnitSpec {
  import PrimerPairMatchType._

  private val pairOneForward1 = new Primer("pair1", "primer1", "GATTACA", "chr1", 1, 7, true)
  private val pairOneReverse1 = new Primer("pair1", "primer2", "GATTACA", "chr1", 1, 57, false)
  private val pairOneForward2 = new Primer("pair1", "primer3", "GATTACA", "chr1", 1, 7, true)
  private val pairTwoForward1 = new Primer("pair2", "primer4", "GATTACA", "chr2", 1, 7, true)
  private val pairTwoReverse1 = new Primer("pair2", "primer5", "GATTACA", "chr2", 1, 57, false)
  private val pairThreeReverse1 = new Primer("pair3", "primer6", "GATTACA", "chr2", 1, 67, false)

  private def toPrimerMatch(p: Primer): Option[PrimerMatch] = Some(LocationBasedPrimerMatch(p, 1))

  "PrimerPairMatchType.apply" should "return NoMatch if no matches for either primer" in {
    PrimerPairMatchType(None, None) shouldBe NoMatch
  }

  it should "return Single if only one primer of the pair had a match" in {
    PrimerPairMatchType(toPrimerMatch(pairOneForward1), None) shouldBe Single
    PrimerPairMatchType(None, toPrimerMatch(pairOneReverse1)) shouldBe Single
  }

  it should "return Canonical if from from the same primer pair id but the primers were different strands" in {
    PrimerPairMatchType(toPrimerMatch(pairOneForward1), toPrimerMatch(pairOneReverse1)) shouldBe Canonical
  }

  it should "return NonCanonical if the same primer pair id and primers on the same strand" in {
    PrimerPairMatchType(toPrimerMatch(pairOneForward1), toPrimerMatch(pairOneForward2)) shouldBe NonCanonical
  }

  it should "return SelfDimer if the same primer pair id and the same primer" in {
    PrimerPairMatchType(toPrimerMatch(pairOneForward1), toPrimerMatch(pairOneForward1)) shouldBe SelfDimer
  }

  it should "return CrossDimer if the different primer pair id and primers on different chromosomes" in {
    PrimerPairMatchType(toPrimerMatch(pairOneForward1), toPrimerMatch(pairTwoReverse1)) shouldBe CrossDimer
  }

  it should "return CrossDimer if the different primer pair id and primers on the same chromosome and outside the insert length" in {
    PrimerPairMatchType(toPrimerMatch(pairTwoForward1), toPrimerMatch(pairThreeReverse1), 68, 100) shouldBe CrossDimer
    PrimerPairMatchType(toPrimerMatch(pairTwoForward1), toPrimerMatch(pairThreeReverse1), 0, 66) shouldBe CrossDimer
  }

  it should "return NonCanonical if the different primer pair id and primers on the same chromosome and within the insert length" in {
    PrimerPairMatchType(toPrimerMatch(pairTwoForward1), toPrimerMatch(pairThreeReverse1)) shouldBe NonCanonical
    PrimerPairMatchType(toPrimerMatch(pairTwoForward1), toPrimerMatch(pairThreeReverse1), 0, 100) shouldBe NonCanonical
    PrimerPairMatchType(toPrimerMatch(pairTwoForward1), toPrimerMatch(pairThreeReverse1), 67, 67) shouldBe NonCanonical
  }

  "PrimerPairMatchType.insertLength" should "give the insert length between two primers" in {
    PrimerPairMatchType.insertLength(LocationBasedPrimerMatch(pairTwoForward1, 1), LocationBasedPrimerMatch(pairThreeReverse1, 1)) shouldBe 67
  }
}
