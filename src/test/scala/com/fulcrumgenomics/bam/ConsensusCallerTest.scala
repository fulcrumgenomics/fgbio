/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.PhredValue
import com.fulcrumgenomics.util.PhredValue.{IntPhredValue, ZeroPhred, ZeroProbability}
import htsjdk.samtools.SAMUtils

/**
  * Tests for ConsensusCaller.
  */
class ConsensusCallerTest extends UnitSpec {

  private def phredToFastq(calls: Seq[PhredValue]): String = {
    SAMUtils.phredToFastq(calls.map(pv => Math.min(SAMUtils.MAX_PHRED_SCORE, pv.toInt)).map(_.toByte).toArray[Byte])
  }

  "ConsensusCaller.adjustBaseQualities" should "cap base qualities" in {
    ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5),
      maxBaseQuality   = 10,
      baseQualityShift = 0,
      errorRatePostUmi = ZeroProbability
    ) should contain theSameElementsInOrderAs Seq(10, 10, 10, 5).map(_.toPhredValue)
  }

  it should "shift base qualities" in {
    ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5),
      maxBaseQuality   = Int.MaxValue,
      baseQualityShift = 10,
      errorRatePostUmi = ZeroProbability
    ) shouldBe Seq(10, 5, 0, 0).map(_.toPhredValue)
  }

  it should "scale base qualities using the post-umi error rate" in {
    ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5),
      maxBaseQuality   = Int.MaxValue,
      baseQualityShift = 0,
      errorRatePostUmi = 10
    ).map(_.toInt) shouldBe Seq(9, 8, 7, 4)
  }

  it should "cap, shift, and scale base qualities" in {
    ConsensusCaller.adjustBaseQualities(
      quals=Seq(20, 15, 10, 5),
      maxBaseQuality   = 10,
      baseQualityShift = 5,
      errorRatePostUmi = 10
    ).map(_.toInt) shouldBe Seq(7, 7, 4, 0)
  }

  "ConsensusCaller.consensusCalls" should "produce a consensus from one read" in {
    val bases = "GATTACA"
    val quals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      bases                   = Seq(bases),
      quals                   = Seq(quals),
      errorRatePreUmi         = ZeroProbability,
      minReads                = 1,
      minConsensusBaseQuality = 0
    )
    cBases shouldBe bases
    cQuals.map(_.toInt) should contain theSameElementsInOrderAs quals.map(_.toInt)
  }

  it should "produce a consensus from two reads" in {
    val bases = "GATTACA"
    val quals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      bases = Seq(bases, bases),
      quals = Seq(quals, quals),
      minReads = 1,
      minConsensusBaseQuality = 0
    )
    val err = PhredValue(10) / 3.0
    val ok = PhredValue(10).inv()
    val numerator = ok * ok
    val denominator = numerator + (3.0 * err * err)
    val expectedQual = (numerator / denominator).inv()
    val expectedQuals = quals.map(q => expectedQual.toInt)
    cBases shouldBe bases
    cQuals.map(_.toInt) should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a consensus from three reads, with one disagreement" in {
    val bases = "GATTACA"
    val otherBases = "GATTTCA"
    val quals = Array(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      bases = Seq(bases, bases, otherBases),
      quals = Seq(quals, quals, quals),
      errorRatePreUmi = ZeroProbability,
      minReads = 1,
      minConsensusBaseQuality = 0
    )
    val err = PhredValue(10) / 3.0
    val ok = PhredValue(10).inv()
    val numeratorAgreement = ok * ok * ok
    val denominatorAgreement = numeratorAgreement + (3.0 * err * err * err)
    val agreementQual = (numeratorAgreement / denominatorAgreement).inv()

    val numeratorDisagreement = ok * ok * err
    val denominatorDisagreement = numeratorDisagreement + (err * err * ok) + (2.0 * err * err * err)
    val disagreementQual = (numeratorDisagreement / denominatorDisagreement).inv()

    val expectedQuals = bases.zip(otherBases).map {
      case (left, right) =>
        if (left == right) agreementQual
        else disagreementQual
    }.map(_.toInt)

    cBases shouldBe bases
    cQuals.map(_.toInt) should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a consensus from two reads of differing lengths" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      bases = Seq("GATTACA", "GATTAC"),
      quals = Seq(quals, quals.slice(0, quals.length-1)),
      minReads = 2,
      minConsensusBaseQuality = 0
    )

    val err = PhredValue(10) / 3.0
    val ok = PhredValue(10).inv()
    val numerator = ok * ok
    val denominator = numerator + (3.0 * err * err)
    val newQual = SAMUtils.phredToFastq((numerator / denominator).inv().toInt).toString

    cBases shouldBe "GATTACN"
    phredToFastq(cQuals) shouldBe ((newQual * 6) + SAMUtils.phredToFastq(ZeroPhred.toInt).toString)
  }

  it should "mask bases with too low of a consensus quality" in {
    val bases = "GATTACA"
    val quals = Array(10, 10, 10, 10, 10, 10, 0)
    val expectedQuals = Array(10, 10, 10, 10, 10, 10, 1)
    val (cBases, cQuals) = ConsensusCaller.consensusCalls(
      bases = Seq(bases),
      quals = Seq(quals.map(_.toPhredValue)),
      errorRatePreUmi = ZeroProbability,
      minReads = 1,
      minConsensusBaseQuality = 10
    )
    cBases shouldBe "GATTACN"
    phredToFastq(cQuals) shouldBe phredToFastq(expectedQuals.map(_.toPhredValue).toSeq)
  }

  "ConsensusCaller.fromBasesAndQualities" should "return None if there are not enough reads" in {
    ConsensusCaller.fromBasesAndQualities(Seq.empty, ConsensusCallerOptions(minReads=1)) shouldBe None
    ConsensusCaller.fromBasesAndQualities(Seq(("GATTACA", "IIIIIII")), ConsensusCallerOptions(minReads=2)) shouldBe None
  }

  it should "throw an exception if the bases and qualities are of a different length" in {
    an[IllegalArgumentException] should be thrownBy ConsensusCaller.fromBasesAndQualities(Seq(("GATTACA", "I")))
    an[IllegalArgumentException] should be thrownBy ConsensusCaller.fromBasesAndQualities(Seq(("G", "IIIIIII")))
  }

  it should "not return a consensus read if the mean consensus quality is too low" in {
    val call1 = ConsensusCaller.fromBasesAndQualities(
      basesAndQualities = Seq(("GATTACA", "AAAAAAA")),
      options           = ConsensusCallerOptions(
        errorRatePreUmi = ZeroProbability,
        minReads=1,
        minMeanConsensusBaseQuality=255
      )
    )
    call1 shouldBe None

    val call2 = ConsensusCaller.fromBasesAndQualities(
      basesAndQualities = Seq(("GATTACA", "AAAAAAA")),
      options           = ConsensusCallerOptions(
        errorRatePreUmi = ZeroProbability,
        minReads=1,
        minMeanConsensusBaseQuality=0
      )
    )
    call2 shouldBe 'defined
  }

  it should "apply the pre-umi-error-rate when it has probability zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val inputQualsString = SAMUtils.phredToFastq(inputQuals.map(_.value.toByte).toArray)
    ConsensusCaller.fromBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = ZeroProbability,
        errorRatePostUmi            = ZeroProbability,
        maxBaseQuality              = 255,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0
      )
    ) match {
      case None => fail
      case Some((cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe inputQualsString
    }
  }

  it should "apply the pre-umi-error-rate when it has probability greater than zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val inputQualsString = SAMUtils.phredToFastq(inputQuals.map(_.value.toByte).toArray)
    val outputQuals = inputQuals.map(phred => ConsensusCaller.probabilityOfErrorTwoTrials(phred, 10))
    val outputQualsString = SAMUtils.phredToFastq(outputQuals.map(_.value.toByte).toArray)
    ConsensusCaller.fromBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = 10,
        errorRatePostUmi            = ZeroProbability,
        maxBaseQuality              = 255,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0
      )
    ) match {
      case None => fail
      case Some((cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe outputQualsString
    }
  }

  it should "apply the post-umi-error-rate when it has probability greater than zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10).map(_.toPhredValue)
    val inputQualsString = SAMUtils.phredToFastq(inputQuals.map(_.value.toByte).toArray)
    val outputQuals = inputQuals.map(phred => ConsensusCaller.probabilityOfErrorTwoTrials(phred, 10))
    val outputQualsString = SAMUtils.phredToFastq(outputQuals.map(_.value.toByte).toArray)
    ConsensusCaller.fromBasesAndQualities(
      basesAndQualities=Seq(("GATTACA", inputQualsString)),
      options = ConsensusCallerOptions(
        errorRatePreUmi             = ZeroProbability,
        errorRatePostUmi            = 10,
        maxBaseQuality              = 255,
        baseQualityShift            = 0,
        minConsensusBaseQuality     = 0,
        minReads                    = 1,
        minMeanConsensusBaseQuality = 0
      )
    ) match {
      case None => fail
      case Some((cBases, cQuals)) =>
        cBases shouldBe "GATTACA"
        cQuals shouldBe outputQualsString
    }
  }
}
