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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.umi.UmiConsensusCaller.SourceRead
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import htsjdk.samtools.util.CloserUtil
import htsjdk.samtools.{SAMRecordSetBuilder, SAMUtils}
import net.jafama.FastMath._
import org.scalatest.OptionValues

/**
  * Tests for ConsensusCaller.
  */
class VanillaUmiConsensusCallerTest extends UnitSpec with OptionValues {
  /** Helper function to make a set of consensus caller options. */
  def cco = VanillaUmiConsensusCallerOptions

  /** Helper function to make a consensus caller. */
  def cc(options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions()) = {
    new VanillaUmiConsensusCaller(options=options, readNamePrefix="testconsensus")
  }

  /** Helper function to make a SourceRead. */
  def src(bases: String, quals: TraversableOnce[Int]) = SourceRead(id="x", bases.getBytes(), quals.toArray.map(_.toByte))

  /** Helper function to make a SourceRead from bases and Phred-33 ascii quals. */
  def src(bases: String, quals: String) = SourceRead(id="x", bases.getBytes(), SAMUtils.fastqToPhred(quals))

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
    val caller = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minConsensusBaseQuality=0.toByte))
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

    val consensus = cc(cco(minReads=1, minConsensusBaseQuality=0.toByte)).consensusCall(sources)
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

    val caller = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(source1, source1, source2))
    consensus shouldBe 'defined
    consensus.get.bases shouldBe source1.bases
    consensus.get.quals should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a shortened consensus from two reads of differing lengths" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10)
    val caller = cc(cco(minReads=2, minConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(src("GATTACA", quals), src("GATTAC", quals.slice(0, quals.length-1))))

    val newQual = expectedConsensusQuality(10, 2)
    val newQuals = Array(newQual, newQual, newQual, newQual, newQual, newQual)

    consensus shouldBe 'defined
    consensus.get.baseString shouldBe "GATTAC"
    consensus.get.quals shouldBe newQuals
  }

  /** This is to test that we don't generate a ton of Q0 or Q2 Ns at the end when we drop below minReads, but
    * instead produce a shortened read.
    */
  it should "produce a consensus even when most of the bases have < minReads" in {
    val src1 = src("A" * 10, (1 to 10).map(q => 30))
    val src2 = src("A" * 20, (1 to 20).map(q => 30))

    val caller = cc(cco(minReads=2, minConsensusBaseQuality=10.toByte))
    val consensus = caller.consensusCall(Seq(src1, src2))

    consensus shouldBe 'defined
    consensus.get.baseString shouldBe "AAAAAAAAAA"
  }


  it should "mask bases with too low of a consensus quality" in {
    val bases = "GATTACA"
    val quals         = Array(10, 10, 10, 10, 10, 10, 5)
    val expectedQuals = Array(10, 10, 10, 10, 10, 10, 2).map(_.toByte)
    val caller    = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minInputBaseQuality=2.toByte, minConsensusBaseQuality=10.toByte))
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
    an[IllegalArgumentException] should be thrownBy cc().consensusCall(Seq(src("GATTACA", Array(20))))
    an[IllegalArgumentException] should be thrownBy cc().consensusCall(Seq(src("G", Array(20,20,20,20,20))))
  }

  it should "apply the pre-umi-error-rate when it has probability zero" in {
    val inputQuals = Seq(10, 10, 10, 10, 10, 10, 10)
    val source = src("GATTACA", inputQuals)
    val opts = cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = PhredScore.MaxValue,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1
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
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1
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
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1
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
    val builder = new SamRecordSetBuilder(readLength=7)
    builder.addFrag(start=100, baseQuality=20).foreach(_.setReadString("GATTACA"))
    builder.addFrag(start=100, baseQuality=30).foreach(_.setReadString("GATTACA"))
    val opts = cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = PhredScore.MaxValue,
      minInputBaseQuality         = 30.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1
    )

    cc(opts).consensusFromSamRecords(builder.toSeq) match {
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
    val consensusCaller = cc(cco(minReads=1, errorRatePreUmi = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.toIterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
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
    val consensusCaller = cc(cco(minReads = 1, errorRatePreUmi  = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.toIterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
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
    val consensusCaller = cc(cco(minReads=1, errorRatePreUmi  = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.toIterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
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

  it should "generate accurate per-read and per-base tags on the consensus reads" in {
    val builder = new SamRecordSetBuilder(readLength=10)
    builder.addFrag("READ1", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.addFrag("READ2", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.addFrag("READ3", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.addFrag("READ4", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.foreach { rec => rec.setReadString("A" * rec.getReadLength) }

    // Monkey with one of the four reads
    builder.find(_.getReadName == "READ2").foreach { r =>
      r.getReadBases()(5) = 'C'.toByte
      r.getBaseQualities()(7) = 5.toByte
    }

    val consensusCaller = cc(cco(minReads = 2, minInputBaseQuality = 20.toByte))
    val consensuses = consensusCaller.consensusReadsFromSamRecords(builder.toSeq)
    consensuses.size shouldBe 1

    val consensus = consensuses.head
    consensus.getReadBases.foreach(_ shouldBe 'A'.toByte)
    consensus.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadCount)  shouldBe Array[Short](4,4,4,4,4,4,4,3,4,4)
    consensus.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadErrors) shouldBe Array[Short](0,0,0,0,0,1,0,0,0,0)
    consensus.getIntegerAttribute(ConsensusTags.PerRead.RawReadCount)     shouldBe 4
    consensus.getIntegerAttribute(ConsensusTags.PerRead.MinRawReadCount)  shouldBe 3
    consensus.getFloatAttribute(ConsensusTags.PerRead.RawReadErrorRate) shouldBe (1 / 39.toFloat)
  }

  it should "not generate per-base tags when turned off" in {
    val builder = new SamRecordSetBuilder(readLength=10)
    builder.addFrag("READ1", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.addFrag("READ2", 0, 1, false).foreach(_.setAttribute(DefaultTag, "AAA"))
    builder.foreach { rec => rec.setReadString("A" * rec.getReadLength) }

    val caller = cc(options=cco(producePerBaseTags=false))
    val consensuses = caller.consensusReadsFromSamRecords(builder.toSeq)
    consensuses.size shouldBe 1

    val consensus = consensuses.head
    consensus.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadCount)  shouldBe null
    consensus.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadErrors) shouldBe null
    consensus.getIntegerAttribute(ConsensusTags.PerRead.RawReadCount)     shouldBe 2
    consensus.getIntegerAttribute(ConsensusTags.PerRead.MinRawReadCount)  shouldBe 2
    consensus.getFloatAttribute(ConsensusTags.PerRead.RawReadErrorRate) shouldBe 0.toFloat
  }

  "VanillaUmiConsensusCaller.filterToMostCommonAlignment" should "return all reads when all cigars are 50M" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="50M") }
    val recs = cc().filterToMostCommonAlignment(builder.toSeq)
    recs should have size 10
  }

  it should "return all reads when all cigars are complicated but the same" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="10M5D10M5I20M5S") }
    val recs = cc().filterToMostCommonAlignment(builder.toSeq)
    recs should have size 10
  }

  it should "return only the 50M reads (i.e. the most common alignment)" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to  3).foreach { i => builder.addFrag(start=100, cigar="25M1D25M") }
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="50M") }
    (1 to  3).foreach { i => builder.addFrag(start=100, cigar="25M2I23M") }

    val recs = cc().filterToMostCommonAlignment(builder.toSeq)
    recs should have size 10
    recs.map(_.getCigarString).distinct shouldBe Seq("50M")
  }

  it should "return only the reads with a single base deletion at base 25" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    // These should all be returned
    (1 to  5).foreach { i => builder.addFrag(start=100, cigar="25M1D25M") }
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="5S20M1D25M") }
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="5S20M1D20M5H") }
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="25M1D20M5S") }

    // These should not be!
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="25M2D25M") }
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="25M1I24M") }
    (1 to  2).foreach { i => builder.addFrag(start=100, cigar="20M1D5M1D25M") }

    val recs = cc().filterToMostCommonAlignment(builder.toSeq)
    recs should have size 11
    recs.map(_.getCigarString).distinct.sorted shouldBe Seq("25M1D25M", "5S20M1D25M", "5S20M1D20M5H", "25M1D20M5S").sorted
  }

  it should "calculate the # of errors relative to the most likely consensus call, even when the final call is an N" in {
    // NB: missing last base on the first read, which causes an N no-call, but errors should still be 1/3
    val call = cc(cco(minReads=4)).consensusCall(Seq(
      src("GATTACAN", Array(20,20,20,20,20,20,20,20)),
      src("GATTACAG", Array(20,20,20,20,20,20,20,20)),
      src("GATTACAG", Array(20,20,20,20,20,20,20,20)),
      src("GATTACAT", Array(20,20,20,20,20,20,20,20))
    )).value

    call.baseString shouldBe "GATTACAN"
    call.depths should contain theSameElementsInOrderAs Seq(4,4,4,4,4,4,4,3)
    call.errors should contain theSameElementsInOrderAs Seq(0,0,0,0,0,0,0,1)
  }
}
