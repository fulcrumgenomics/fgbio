/*
 * The MIT License
 *
 * Copyright (c) $year Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.LogDouble
import com.fulcrumgenomics.util.LogDouble._
import com.fulcrumgenomics.util.LogProbability._
import com.fulcrumgenomics.util.PhredScore._


/**
  * Tests for ConsensusCaller.
  */
class ConsensusCallerTest extends UnitSpec {

  private def qualtiesEqual(actual: Seq[LogDouble], expected: Seq[Int]): Boolean = {
    if (actual.length != expected.length) false
    else actual.zip(expected.map(_.fromPhredScore)).exists { case (act, exp) => Math.abs(act.logValue - exp.logValue) > 0.00001 }
  }

  "ConsensusCaller.adjustBaseQualities" should "cap base qualities" in {
    val quals = ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5).map(_.fromPhredScore),
      maxBaseQuality   = 10.fromPhredScore,
      baseQualityShift = 0,
      errorRatePostUmi = ZeroProbability
    )
    qualtiesEqual(quals, Seq(10, 10, 10, 5))
  }

  it should "shift base qualities" in {
    val quals = ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5).map(_.fromPhredScore),
      maxBaseQuality   = Int.MaxValue.fromPhredScore,
      baseQualityShift = 10,
      errorRatePostUmi = ZeroProbability
    )
    qualtiesEqual(quals, Seq(10, 5, 0, 0))
  }

  it should "scale base qualities using the post-umi error rate" in {
    val quals = ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5).map(_.fromPhredScore),
      maxBaseQuality   = Int.MaxValue.fromPhredScore,
      baseQualityShift = 0,
      errorRatePostUmi = 10.fromPhredScore
    )
    qualtiesEqual(quals, Seq(9, 8, 7, 4))
  }

  it should "cap, shift, and scale base qualities" in {
    val quals = ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5).map(_.fromPhredScore),
      maxBaseQuality   = 10.fromPhredScore,
      baseQualityShift = 5,
      errorRatePostUmi = 10.fromPhredScore
    )
    qualtiesEqual(quals, Seq(7, 7, 4, 0))
  }

  "ConsensusCaller.consensusCalls" should "produce a consensus from one read" in {
    val bases = "GATTACA"
    val quals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      baseStrings             = Seq(bases),
      qualSeqs                = Seq(quals),
      errorRatePreUmi         = ZeroProbability,
      minReads                = 1,
      minConsensusBaseQuality = 0.fromPhredScore
    )
    cBases shouldBe bases
    cQuals.map(_.toPhredScoreInt) should contain theSameElementsInOrderAs quals.map(_.toPhredScoreInt)
  }

  it should "produce a consensus from two reads" in {
    val bases = "GATTACA"
    val quals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      baseStrings = Seq(bases, bases),
      qualSeqs = Seq(quals, quals),
      minReads = 1,
      minConsensusBaseQuality = 0.fromPhredScore
    )
    val err = 10.fromPhredScore / 3.0.toLogDouble
    val ok  = 10.fromPhredScore.oneMinus()
    val numerator = ok * ok
    val denominator = numerator + (3.0.toLogDouble * err * err)
    val expectedQual = (numerator / denominator).oneMinus()
    val expectedQuals = quals.map(q => expectedQual.toPhredScoreInt)
    cBases shouldBe bases
    cQuals.map(_.toPhredScoreInt) should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a consensus from three reads, with one disagreement" in {
    val bases = "GATTACA"
    val otherBases = "GATTTCA"
    val quals = Array(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      baseStrings = Seq(bases, bases, otherBases),
      qualSeqs = Seq(quals, quals, quals),
      errorRatePreUmi = ZeroProbability,
      minReads = 1,
      minConsensusBaseQuality = 0.fromPhredScore
    )
    val err = 10.fromPhredScore / 3.0.toLogDouble
    val ok  = 10.fromPhredScore.oneMinus()
    val numeratorAgreement = ok * ok * ok
    val denominatorAgreement = numeratorAgreement + (3.0.toLogDouble * err * err * err)
    val agreementQual = (numeratorAgreement / denominatorAgreement).oneMinus()

    val numeratorDisagreement = ok * ok * err
    val denominatorDisagreement = numeratorDisagreement + (err * err * ok) + (2.0.toLogDouble * err * err * err)
    val disagreementQual = (numeratorDisagreement / denominatorDisagreement).oneMinus()

    val expectedQuals = bases.zip(otherBases).map {
      case (left, right) =>
        if (left == right) agreementQual
        else disagreementQual
  }.map(_.toPhredScoreInt)

    cBases shouldBe bases
    cQuals.map(_.toPhredScoreInt) should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a consensus from two reads of differing lengths" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      baseStrings = Seq("GATTACA", "GATTAC"),
      qualSeqs = Seq(quals, quals.slice(0, quals.length-1)),
      minReads = 2,
      minConsensusBaseQuality = 0.fromPhredScore
    )

    val err = 10.fromPhredScore / 3.0.toLogDouble
    val ok = 10.fromPhredScore.oneMinus()
    val numerator = ok * ok
    val denominator = numerator + (3.0.toLogDouble * err * err)
    val newQual = (numerator / denominator).oneMinus().toPhredScoreChar.toString

    cBases shouldBe "GATTACN"
    cQuals.map(_.toPhredScoreChar).mkString shouldBe ((newQual * 6) + PhredZeroChar.toString)
  }

  it should "mask bases with too low of a consensus quality" in {
    val bases = "GATTACA"
    val quals = Array(10, 10, 10, 10, 10, 10, 0)
    val expectedQuals = Array(10, 10, 10, 10, 10, 10, 1)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      baseStrings = Seq(bases),
      qualSeqs = Seq(quals.map(_.fromPhredScore)),
      errorRatePreUmi = ZeroProbability,
      minReads = 1,
      minConsensusBaseQuality = 10.fromPhredScore
    )
    cBases shouldBe "GATTACN"
    cQuals.map(_.toPhredScoreChar).mkString shouldBe expectedQuals.map(_.fromPhredScore.toPhredScoreChar).mkString
  }

  "ConsensusCaller.consensusFromStringBasesAndQualities" should "return None if there are not enough reads" in {
    ConsensusCaller.consensusFromStringBasesAndQualities(Seq.empty, ConsensusCallerOptions(minReads=1)) shouldBe None
    ConsensusCaller.consensusFromStringBasesAndQualities(Seq(("GATTACA", "IIIIIII")), ConsensusCallerOptions(minReads=2)) shouldBe None
  }

  it should "throw an exception if the bases and qualities are of a different length" in {
    an[IllegalArgumentException] should be thrownBy ConsensusCaller.consensusFromStringBasesAndQualities(Seq(("GATTACA", "I")))
    an[IllegalArgumentException] should be thrownBy ConsensusCaller.consensusFromStringBasesAndQualities(Seq(("G", "IIIIIII")))
  }

  it should "not return a consensus read if the mean consensus quality is too low" in {
    val call1 = ConsensusCaller.consensusFromStringBasesAndQualities(
      basesAndQualities = Seq(("GATTACA", "AAAAAAA")),
      options           = ConsensusCallerOptions(
        errorRatePreUmi = ZeroProbability,
        minReads=1,
        minMeanConsensusBaseQuality=255.fromPhredScore
      )
    )
    call1 shouldBe None

    val call2 = ConsensusCaller.consensusFromStringBasesAndQualities(
      basesAndQualities = Seq(("GATTACA", "AAAAAAA")),
      options           = ConsensusCallerOptions(
        errorRatePreUmi = ZeroProbability,
        minReads=1,
        minMeanConsensusBaseQuality=0.fromPhredScore
      )
    )
    call2 shouldBe 'defined
  }

  it should "apply the pre-umi-error-rate when it has probability zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val inputQualsString = inputQuals.map(_.toPhredScoreChar).mkString
    ConsensusCaller.consensusFromStringBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = ZeroProbability,
        errorRatePostUmi            = ZeroProbability,
        maxBaseQuality              = 255.fromPhredScore,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0.fromPhredScore,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0.fromPhredScore
      )
    ) match {
      case None => fail
      case Some(ConsensusRead(cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe inputQualsString
    }
  }

  it should "apply the pre-umi-error-rate when it has probability greater than zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val inputQualsString = inputQuals.map(_.toPhredScoreChar).mkString
    val outputQuals = inputQuals.map(phred => ConsensusCaller.probabilityOfErrorTwoTrials(phred, 10.fromPhredScore))
    val outputQualsString = outputQuals.map(_.toPhredScoreChar).mkString
    ConsensusCaller.consensusFromStringBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = 10.fromPhredScore,
        errorRatePostUmi            = ZeroProbability,
        maxBaseQuality              = 255.fromPhredScore,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0.fromPhredScore,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0.fromPhredScore
      )
    ) match {
      case None => fail
      case Some(ConsensusRead(cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe outputQualsString
    }
  }

  it should "apply the post-umi-error-rate when it has probability greater than zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.fromPhredScore)
    val inputQualsString = inputQuals.map(_.toPhredScoreChar).mkString
    val outputQuals = inputQuals.map(phred => ConsensusCaller.probabilityOfErrorTwoTrials(phred, 10.fromPhredScore))
    val outputQualsString = outputQuals.map(_.toPhredScoreChar).mkString
    ConsensusCaller.consensusFromStringBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = ZeroProbability,
        errorRatePostUmi            = 10.fromPhredScore,
        maxBaseQuality              = 255.fromPhredScore,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0.fromPhredScore,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0.fromPhredScore
      )
    ) match {
      case None => fail
      case Some(ConsensusRead(cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe outputQualsString
    }
  }

  // TODO: test requireConsensusForBothPairs
}
