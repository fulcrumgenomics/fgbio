/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import htsjdk.samtools.util.CloserUtil
import htsjdk.samtools.{SAMFileHeader, SAMRecordSetBuilder, SAMUtils}
import net.jafama.FastMath._

import scala.collection.JavaConversions._
import scala.collection.JavaConverters._

/**
  * Tests for ConsensusCaller.
  */
class VanillaUmiConsensusCallerTest extends UnitSpec {
  /** Helper function to make a set of consensus caller options. */
  def cco = VanillaUmiConsensusCallerOptions

  /** Helper function to make a consensus caller. */
  def cc(options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions()) = {
    new VanillaUmiConsensusCaller(Iterator.empty, new SAMFileHeader, options=options)
  }

  /** Helper function to make a SourceRead. */
  def src(bases: String, quals: TraversableOnce[Int]) = SourceRead(bases.getBytes(), quals.toArray.map(_.toByte))

  /** Helper function to make a SourceRead from bases and Phred-33 ascii quals. */
  def src(bases: String, quals: String) = SourceRead(bases.getBytes(), SAMUtils.fastqToPhred(quals))

  /**
    * Function to calculated the expected quality of a consensus base in non-log math, that should work for
    * modest values of Q and N.
    *
    * @param q The quality score of the correct bases (assumed all the same)
    * @param n The number of observations at that quality
    * @return the phred-scaled number (byte) of the consensus base
    */
  def expectedConsensusQuality(q: Int, n: Int): Byte = {
    val p   = BigDecimal(pow(10.0, q  / -10.0))
    val ok  = BigDecimal(1) - p
    val err = p / 3 // error could be one of three bases

    val numerator   = ok.pow(n)
    val denomenator = numerator + (err.pow(2) * 3)
    val pError = 1 - (numerator / denomenator)
    val phred = -10 * log10(pError.toDouble)
    phred.toByte
  }

  "VanillaUmiConsensusCaller.consensusCalls" should "produce a consensus from one read" in {
    val source = src("GATTACA", Seq(10, 10, 10, 10, 10, 10, 10))
    val caller = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, rawBaseQualityShift=0.toByte, minConsensusBaseQuality=0.toByte, minMeanConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(source))
    consensus shouldBe 'defined
    consensus.get.bases shouldBe source.bases
    consensus.get.quals shouldBe source.quals
  }

  it should "produce a consensus from two reads" in {
    val source    = src("GATTACA", Seq(10, 10, 10, 10, 10, 10, 10))
    val sources   = Seq(source, source)

    val expectedQual = expectedConsensusQuality(10, 2)
    val expectedQuals = source.quals.map(q => expectedQual)

    val consensus = cc(cco(minReads=1, minConsensusBaseQuality=0.toByte, rawBaseQualityShift=0.toByte)).consensusCall(sources)
    consensus shouldBe 'defined
    consensus.get.bases shouldBe source.bases
    consensus.get.quals should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a consensus from three reads, with one disagreement" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10)
    val source1 = src("GATTACA", quals)
    val source2 = src("GATTTCA", quals)
    val err = LogProbability.normalizeByScalar(LogProbability.fromPhredScore(10), 3)
    val ok  = LogProbability.not(LogProbability.fromPhredScore(10))
    val numeratorAgreement = LogProbability.and(Array(ok, ok, ok))
    val denominatorAgreement = LogProbability.or(numeratorAgreement, LogProbability.and(Array(LnThree, err, err, err)))
    val agreementQual = LogProbability.not(LogProbability.normalizeByLogProbability(numeratorAgreement, denominatorAgreement))

    val numeratorDisagreement =  LogProbability.and(Array(ok, ok, err))
    val denominatorDisagreement = LogProbability.or(Array(
      numeratorDisagreement,
      LogProbability.and(Array(err, err, ok)),
      LogProbability.and(Array(LnThree, err, err, err)))
    )
    val disagreementQual = LogProbability.not(numeratorDisagreement - denominatorDisagreement)

    val expectedQuals = source1.bases.zip(source2.bases).map {
      case (left, right) =>
        if (left == right) agreementQual
        else disagreementQual
    }.map(PhredScore.fromLogProbability)

    val caller = cc(cco(errorRatePreUmi=PhredScore.MaxValue, rawBaseQualityShift=0.toByte, minReads=1, minConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(source1, source1, source2))
    consensus shouldBe 'defined
    consensus.get.bases shouldBe source1.bases
    consensus.get.quals should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a shortened consensus from two reads of differing lengths" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10)
    val caller = cc(cco(minReads=2, minConsensusBaseQuality=0.toByte, rawBaseQualityShift=0.toByte))
    val consensus = caller.consensusCall(Seq(src("GATTACA", quals), src("GATTAC", quals.slice(0, quals.length-1))))

    val newQual = expectedConsensusQuality(10, 2)
    val newQuals = Array(newQual, newQual, newQual, newQual, newQual, newQual)

    consensus shouldBe 'defined
    consensus.get.baseString shouldBe "GATTAC"
    consensus.get.quals shouldBe newQuals
  }

  /** This is to test that we don't generate a ton of Q0 or Q2 Ns at the end when we drop below minReads, deflate
    * the average consensus base quality and then discard the consensus read.
    */
  it should "produce a consensus even when most of the bases have < minReads" in {
    val src1 = src("A" * 10, (1 to 10).map(q => 30))
    val src2 = src("A" * 20, (1 to 20).map(q => 30))

    val caller = cc(cco(minReads=2, minConsensusBaseQuality=10.toByte, rawBaseQualityShift=0.toByte, minMeanConsensusBaseQuality=25.toByte))
    val consensus = caller.consensusCall(Seq(src1, src2))

    consensus shouldBe 'defined
    consensus.get.baseString shouldBe "AAAAAAAAAA"
  }


  it should "mask bases with too low of a consensus quality" in {
    val bases = "GATTACA"
    val quals         = Array(10, 10, 10, 10, 10, 10, 5)
    val expectedQuals = Array(10, 10, 10, 10, 10, 10, 2).map(_.toByte)
    val caller    = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minConsensusBaseQuality=10.toByte,
                           minMeanConsensusBaseQuality=PhredScore.MinValue, rawBaseQualityShift=0.toByte))
    val consensus = caller.consensusCall(Seq(src(bases, quals)))
    consensus shouldBe 'defined
    consensus.get.baseString shouldBe "GATTACN"
    consensus.get.quals shouldBe expectedQuals
  }

  it should "return None if there are not enough reads" in {
    cc(cco(minReads=1)).consensusCall(Seq.empty) shouldBe None
    cc(cco(minReads=2)).consensusCall(Seq(src("GATTACA", Array(20,20,20,20,20,20,20)))) shouldBe None
  }

  it should "throw an exception if the bases and qualities are of a different length" in {
    an[AssertionError] should be thrownBy cc().consensusCall(Seq(src("GATTACA", Array(20))))
    an[AssertionError] should be thrownBy cc().consensusCall(Seq(src("G", Array(20,20,20,20,20))))
  }

  it should "not return a consensus read if the mean consensus quality is too low" in {
    val call1 = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minMeanConsensusBaseQuality=PhredScore.MaxValue))
                  .consensusCall(Seq(src("GATTACA", "AAAAAAA")))
    call1 shouldBe None

    val call2 = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minMeanConsensusBaseQuality=PhredScore.MinValue))
                  .consensusCall(Seq(src("GATTACA", "AAAAAAA")))
    call2 shouldBe 'defined
  }

  it should "apply the pre-umi-error-rate when it has probability zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10)
    val source = src("GATTACA", inputQuals)
    val opts = cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = PhredScore.MaxValue,
      maxRawBaseQuality           = PhredScore.MaxValue,
      rawBaseQualityShift         = 0.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1,
      minMeanConsensusBaseQuality = PhredScore.MinValue
    )

    cc(opts).consensusCall(Seq(source)) match {
      case None => fail
      case Some(consensus) =>
        consensus.baseString shouldBe "GATTACA"
        consensus.quals shouldBe inputQuals.map(_.toByte)
    }
  }

  it should "apply the pre-umi-error-rate when it has probability greater than zero" in {
    val caller = cc(cco(
      errorRatePreUmi             = 10.toByte,
      errorRatePostUmi            = PhredScore.MaxValue,
      maxRawBaseQuality           = PhredScore.MaxValue,
      rawBaseQualityShift         = 0.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1,
      minMeanConsensusBaseQuality = PhredScore.MinValue
    ))

    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10)
    val lnProbError = LogProbability.fromPhredScore(10)
    val outputQual  = PhredScore.fromLogProbability(LogProbability.probabilityOfErrorTwoTrials(lnProbError, lnProbError))
    val outputQuals = inputQuals.map(q => outputQual)
    caller.consensusCall(Seq(src("GATTACA", inputQuals))) match {
      case None => fail
      case Some(consensus) =>
        consensus.baseString shouldBe "GATTACA"
        consensus.quals      shouldBe outputQuals
    }
  }

  it should "apply the post-umi-error-rate when it has probability greater than zero" in {
    val caller = cc(cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = 10.toByte,
      maxRawBaseQuality           = PhredScore.MaxValue,
      rawBaseQualityShift         = 0.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1,
      minMeanConsensusBaseQuality = PhredScore.MinValue
    ))

    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10)
    val lnProbError = LogProbability.fromPhredScore(10)
    val outputQual  = PhredScore.fromLogProbability(LogProbability.probabilityOfErrorTwoTrials(lnProbError, lnProbError))
    val outputQuals = inputQuals.map(q => outputQual)

    caller.consensusCall(Seq(src("GATTACA", inputQuals))) match {
      case None => fail
      case Some(consensus) =>
        consensus.baseString shouldBe "GATTACA"
        consensus.quals      shouldBe outputQuals
    }
  }

  it should "apply the minInputBaseQuality appropriately" in {
    val sources    = Seq(src("GATTACA", Array.tabulate(7)(_ => 20)), src("GATTACA", Array.tabulate(7)(_ => 30)))
    val opts = cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = PhredScore.MaxValue,
      minInputBaseQuality         = 30.toByte,
      maxRawBaseQuality           = PhredScore.MaxValue,
      rawBaseQualityShift         = 0.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1,
      minMeanConsensusBaseQuality = PhredScore.MinValue
    )

    cc(opts).consensusCall(sources) match {
      case None => fail
      case Some(consensus) =>
        consensus.baseString shouldBe "GATTACA"
        consensus.quals.foreach(q => q shouldBe 30.toByte)
    }
  }


  it should "should create two consensus for two UMI groups" in {
    val builder = new SAMRecordSetBuilder()
    builder.addFrag("READ1", 0, 1, false).setAttribute(DefaultTag, "GATTACA")
    builder.addFrag("READ2", 0, 1, false).setAttribute(DefaultTag, "GATTACA")
    builder.addFrag("READ3", 0, 1, false).setAttribute(DefaultTag, "ACATTAG")
    builder.addFrag("READ4", 0, 1, false).setAttribute(DefaultTag, "ACATTAG")
    builder.getRecords.foreach { rec =>
      rec.setReadString("A" * rec.getReadLength)
      rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new VanillaUmiConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new VanillaUmiConsensusCallerOptions(
        minReads         = 1,
        errorRatePreUmi  = PhredScore.MaxValue,
        errorRatePostUmi = PhredScore.MaxValue
      )
    )
    consensusCaller.hasNext shouldBe true
    val calls = consensusCaller.toList
    consensusCaller.hasNext shouldBe false
    calls should have size 2
    CloserUtil.close(reader)
  }

  it should "should create two consensus for a read pair" in {
    val builder = new SAMRecordSetBuilder()
    builder.addPair("READ1", 0, 1, 1000)
    builder.getRecords.foreach {
      rec =>
        rec.setAttribute(DefaultTag, "GATTACA")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new VanillaUmiConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new VanillaUmiConsensusCallerOptions(
        minReads         = 1,
        errorRatePreUmi  = PhredScore.MaxValue,
        errorRatePostUmi = PhredScore.MaxValue
      )
    )
    consensusCaller.hasNext shouldBe true
    val calls = consensusCaller.toList
    consensusCaller.hasNext shouldBe false
    calls should have size 2
    calls.foreach { rec =>
      rec.getReadPairedFlag shouldBe true
    }
    calls.head.getFirstOfPairFlag shouldBe true
    calls.last.getSecondOfPairFlag shouldBe true
    calls.head.getReadName shouldBe calls.last.getReadName
    CloserUtil.close(reader)
  }

  it should "should create four consensus for two read pairs with different group ids" in {
    val builder = new SAMRecordSetBuilder()
    builder.addPair("READ1", 0, 1, 1000)
    builder.addPair("READ2", 1, 1, 1000)

    builder.getRecords.slice(0, 2).foreach {
      rec =>
        rec.setAttribute(DefaultTag, "GATTACA")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    builder.getRecords.slice(2, 4).foreach {
      rec =>
        rec.setAttribute(DefaultTag, "ACATTAG")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new VanillaUmiConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new VanillaUmiConsensusCallerOptions(
        minReads         = 1,
        errorRatePreUmi  = PhredScore.MaxValue,
        errorRatePostUmi = PhredScore.MaxValue
      )
    )
    consensusCaller.hasNext shouldBe true
    val calls = consensusCaller.toList
    consensusCaller.hasNext shouldBe false
    calls should have size 4
    calls.foreach { rec =>
      rec.getReadPairedFlag shouldBe true
    }
    calls.map(_.getFirstOfPairFlag) should contain theSameElementsInOrderAs Seq(true, false, true, false)
    calls.map(_.getSecondOfPairFlag) should contain theSameElementsInOrderAs Seq(false, true, false, true)
    calls(0).getReadName shouldBe calls(1).getReadName
    calls(2).getReadName shouldBe calls(3).getReadName
    calls(0).getReadName should not be calls(2).getReadName
    CloserUtil.close(reader)
  }
}
