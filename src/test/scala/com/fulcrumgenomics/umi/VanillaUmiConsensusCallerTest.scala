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

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
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
    val caller    = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minConsensusBaseQuality=10.toByte))
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
    val sources    = Seq(src("GATTACA", Array.tabulate(7)(_ => 20)), src("GATTACA", Array.tabulate(7)(_ => 30)))
    val opts = cco(
      errorRatePreUmi             = PhredScore.MaxValue,
      errorRatePostUmi            = PhredScore.MaxValue,
      minInputBaseQuality         = 30.toByte,
      minConsensusBaseQuality     = PhredScore.MinValue,
      minReads                    = 1
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

  it should "create a dummy read when not requiring pairs and one side has fewer reads for some reason" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=30)
    builder.addPair("q1", start1=100, start2=200)
    builder.addPair("q2", start1=100, start2=200)
    builder.addPair("q3", start1=100, start2=200)
    builder.foreach { rec =>
      rec.setReadString("A" * rec.getReadLength)
      rec.setAttribute("MI", "OneAndOnly")
    }

    // Create an iterator that will return all the R1s but only one of the R2s
    val iterator = builder.iterator.filterNot(r => r.getSecondOfPairFlag && r.getReadName > "q1")
    val caller = new VanillaUmiConsensusCaller(input=iterator, header=builder.header, options=cco(minReads=3, minInputBaseQuality=30.toByte, requireConsensusForBothPairs=false))
    val recs = caller.toSeq
    recs.size shouldBe 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("There should be an R1 in the consensus reads."))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("There should be an R2 in the consensus reads."))
    r1.getReadString shouldBe ("A" * 10)
    r2.getReadString shouldBe ("N" * 10)
  }

  "VanillaUmiConsensusCaller.filterToMostCommonAlignment" should "return all reads when all cigars are 50M" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="50M") }
    val recs = VanillaUmiConsensusCaller.filterToMostCommonAlignment(builder.toSeq)
    recs should have size 10
  }

  it should "return all reads when all cigars are complicated but the same" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="10M5D10M5I20M5S") }
    val recs = VanillaUmiConsensusCaller.filterToMostCommonAlignment(builder.toSeq)
    recs should have size 10
  }

  it should "return only the 50M reads (i.e. the most common alignment)" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    (1 to  3).foreach { i => builder.addFrag(start=100, cigar="25M1D25M") }
    (1 to 10).foreach { i => builder.addFrag(start=100, cigar="50M") }
    (1 to  3).foreach { i => builder.addFrag(start=100, cigar="25M2I23M") }

    val recs = VanillaUmiConsensusCaller.filterToMostCommonAlignment(builder.toSeq)
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

    val recs = VanillaUmiConsensusCaller.filterToMostCommonAlignment(builder.toSeq)
    recs should have size 11
    recs.map(_.getCigarString).distinct.sorted shouldBe Seq("25M1D25M", "5S20M1D25M", "5S20M1D20M5H", "25M1D20M5S").sorted
  }
}
