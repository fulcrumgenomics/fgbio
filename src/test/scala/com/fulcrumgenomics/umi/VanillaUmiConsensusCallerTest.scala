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

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.UmiConsensusCaller.SourceRead
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import htsjdk.samtools.SAMUtils
import htsjdk.samtools.util.CloserUtil
import org.scalatest.OptionValues
import org.apache.commons.math3.util.FastMath._

import scala.collection.mutable.ArrayBuffer

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
  def src(bases: String, quals: TraversableOnce[Int]) = SourceRead(id="x", bases.getBytes(), quals.toArray.map(_.toByte), Cigar(bases.length + "M"))

  /** Helper function to make a SourceRead from bases and Phred-33 ascii quals. */
  def src(bases: String, quals: String) = SourceRead(id="x", bases.getBytes(), SAMUtils.fastqToPhred(quals), Cigar(bases.length + "M"))

  /** Helper function to make a SourceRead from bases and Phred-33 ascii quals. */
  def src(cigar: String) = {
    val cig = Cigar(cigar)
    val len = cig.lengthOnQuery
    SourceRead(id="x", ("A"*len).getBytes, Array.tabulate(len)(_ => 30.toByte), cigar=cig)
  }

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
    val denominator = numerator + (err.pow(2) * 3)
    val pError = 1 - (numerator / denominator)
    val phred = -10 * log10(pError.toDouble)
    phred.toByte
  }

  "VanillaUmiConsensusCaller.consensusCalls" should "produce a consensus from one read" in {
    val source = src("GATTACA", Seq(10, 10, 10, 10, 10, 10, 10))
    val caller = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(source))
    consensus shouldBe 'defined
    consensus.value.bases shouldBe source.bases
    consensus.value.quals shouldBe source.quals
  }

  it should "produce a consensus from two reads" in {
    val source    = src("GATTACA", Seq(10, 10, 10, 10, 10, 10, 10))
    val sources   = Seq(source, source)

    val expectedQual = expectedConsensusQuality(10, 2)
    val expectedQuals = source.quals.map(q => expectedQual)

    val consensus = cc(cco(minReads=1, minConsensusBaseQuality=0.toByte)).consensusCall(sources)
    consensus shouldBe 'defined
    consensus.value.bases shouldBe source.bases
    consensus.value.quals should contain theSameElementsInOrderAs expectedQuals
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
    consensus.value.bases shouldBe source1.bases
    consensus.value.quals should contain theSameElementsInOrderAs expectedQuals
  }

  it should "produce a shortened consensus from two reads of differing lengths" in {
    val quals = Array(10, 10, 10, 10, 10, 10, 10)
    val caller = cc(cco(minReads=2, minConsensusBaseQuality=0.toByte))
    val consensus = caller.consensusCall(Seq(src("GATTACA", quals), src("GATTAC", quals.slice(0, quals.length-1))))

    val newQual = expectedConsensusQuality(10, 2)
    val newQuals = Array(newQual, newQual, newQual, newQual, newQual, newQual)

    consensus shouldBe 'defined
    consensus.value.baseString shouldBe "GATTAC"
    consensus.value.quals shouldBe newQuals
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
    consensus.value.baseString shouldBe "AAAAAAAAAA"
  }

  it should "downsample input reads so that each consensus is made from <= max reads" in {
    val r = src("AAAAAAAAAA", "##########")

    for (max <- Seq(3, 1000); n <- Range.inclusive(1, 10)) {
      val srcs      = Seq.tabulate(n)(_ => r)
      val caller    = cc(cco(minReads=1, maxReads=max))
      val consensus = caller.consensusCall(srcs)

      consensus match {
        case None                => fail("Consensus should have been generated")
        case Some(c) if n <= max => c.depths.forall(_ <= n)   shouldBe true
        case Some(c)             => c.depths.forall(_ <= max) shouldBe true
      }
    }
  }

  it should "mask bases with too low of a consensus quality" in {
    val bases = "GATTACA"
    val quals         = Array(10, 10, 10, 10, 10, 10, 5)
    val expectedQuals = Array(10, 10, 10, 10, 10, 10, 2).map(_.toByte)
    val caller    = cc(cco(errorRatePreUmi=PhredScore.MaxValue, minReads=1, minInputBaseQuality=2.toByte, minConsensusBaseQuality=10.toByte))
    val consensus = caller.consensusCall(Seq(src(bases, quals)))
    consensus shouldBe 'defined
    consensus.value.baseString shouldBe "GATTACN"
    consensus.value.quals shouldBe expectedQuals
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
    val builder = new SamBuilder(readLength=7)
    builder.addFrag(start=100, bases="GATTACA", quals="5555555", attrs=Map("MI" -> "1")) // Q20
    builder.addFrag(start=100, bases="GATTACA", quals="???????", attrs=Map("MI" -> "1")) // Q30
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
    val builder = new SamBuilder()
    builder.addFrag(name="READ1", start=1, attrs=Map(DefaultTag -> "GATTACA"), bases="A"*builder.readLength, quals="|"*builder.readLength)
    builder.addFrag(name="READ2", start=1, attrs=Map(DefaultTag -> "GATTACA"), bases="A"*builder.readLength, quals="|"*builder.readLength)
    builder.addFrag(name="READ3", start=1, attrs=Map(DefaultTag -> "ACATTAG"), bases="A"*builder.readLength, quals="|"*builder.readLength)
    builder.addFrag(name="READ4", start=1, attrs=Map(DefaultTag -> "ACATTAG"), bases="A"*builder.readLength, quals="|"*builder.readLength)

    val reader = builder.toSource
    val consensusCaller = cc(cco(minReads=1, errorRatePreUmi = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.iterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
    calls should have size 2
    CloserUtil.close(reader)
  }

  it should "should create two consensus for a read pair" in {
    val len = 100
    val builder = new SamBuilder(readLength = len)
    builder.addPair(name="READ1", start1=1, start2=1000, bases1="A"*len, bases2="A"*len, quals1="|"*len, quals2="|"*len, attrs=Map(DefaultTag -> "GATTACA"))
    val reader = builder.toSource
    val consensusCaller = cc(cco(minReads = 1, errorRatePreUmi  = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.iterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
    calls should have size 2
    calls.foreach { rec =>
      rec.paired shouldBe true
    }
    calls.head.firstOfPair shouldBe true
    calls.last.secondOfPair shouldBe true
    calls.head.name shouldBe calls.last.name
    CloserUtil.close(reader)
  }

  it should "should create four consensus for two read pairs with different group ids" in {
    val len     = 100
    val builder = new SamBuilder()
    builder.addPair("READ1", start1=1, start2=1000, bases1="A"*len, bases2="A"*len, quals1="|"*len, quals2="|"*len, attrs=Map(DefaultTag -> "GATTACA"))
    builder.addPair("READ2", start1=1, start2=1000, bases1="A"*len, bases2="A"*len, quals1="|"*len, quals2="|"*len, attrs=Map(DefaultTag -> "ACATTAG"))

    val reader = builder.toSource
    val consensusCaller = cc(cco(minReads=1, errorRatePreUmi  = PhredScore.MaxValue, errorRatePostUmi = PhredScore.MaxValue))
    val iterator = new ConsensusCallingIterator(reader.iterator, consensusCaller)

    iterator.hasNext shouldBe true
    val calls = iterator.toList
    iterator.hasNext shouldBe false
    calls should have size 4
    calls.foreach { rec =>
      rec.paired shouldBe true
    }
    calls.map(_.firstOfPair) should contain theSameElementsInOrderAs Seq(true, false, true, false)
    calls.map(_.secondOfPair) should contain theSameElementsInOrderAs Seq(false, true, false, true)
    calls(0).name shouldBe calls(1).name
    calls(2).name shouldBe calls(3).name
    calls(0).name should not be calls(2).name
    CloserUtil.close(reader)
  }

  it should "generate accurate per-read and per-base tags on the consensus reads" in {
    val builder = new SamBuilder(readLength=10)
    builder.addFrag("READ1", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "GAT-ACA"))
    builder.addFrag("READ2", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "GAT-ACA"))
    builder.addFrag("READ3", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "GAT-ACA"))
    builder.addFrag("READ4", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "CCC-TTT"))

    // Monkey with one of the four reads
    builder.find(_.name == "READ2").foreach { r =>
      r.bases(5) = 'C'.toByte
      r.quals(7) = 5.toByte
    }

    val consensusCaller = cc(cco(minReads = 2, minInputBaseQuality = 20.toByte))
    val consensuses = consensusCaller.consensusReadsFromSamRecords(builder.toSeq)
    consensuses.size shouldBe 1

    val consensus = consensuses.head
    consensus.bases.foreach(_ shouldBe 'A'.toByte)
    consensus[Array[Short]](ConsensusTags.PerBase.RawReadCount)  shouldBe Array[Short](4,4,4,4,4,4,4,3,4,4)
    consensus[Array[Short]](ConsensusTags.PerBase.RawReadErrors) shouldBe Array[Short](0,0,0,0,0,1,0,0,0,0)
    consensus[Int](ConsensusTags.PerRead.RawReadCount)     shouldBe 4
    consensus[Int](ConsensusTags.PerRead.MinRawReadCount)  shouldBe 3
    consensus[Float](ConsensusTags.PerRead.RawReadErrorRate) shouldBe (1 / 39.toFloat)
    consensus[String](DefaultTag) shouldBe "AAA"
    consensus[String](ConsensusTags.UmiBases) shouldBe "GAT-ACA"
  }

  it should "not generate per-base tags when turned off" in {
    val builder = new SamBuilder(readLength=10)
    builder.addFrag("READ1", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "GAT-ACA"))
    builder.addFrag("READ2", start=1, bases="A"*builder.readLength, attrs=Map(DefaultTag -> "AAA", ConsensusTags.UmiBases -> "GAT-ACA"))

    val caller = cc(options=cco(producePerBaseTags=false))
    val consensuses = caller.consensusReadsFromSamRecords(builder.toSeq)
    consensuses.size shouldBe 1

    val consensus = consensuses.head
    consensus.get[Array[Short]](ConsensusTags.PerBase.RawReadCount)  shouldBe None
    consensus.get[Array[Short]](ConsensusTags.PerBase.RawReadErrors) shouldBe None
    consensus[Int](ConsensusTags.PerRead.RawReadCount)     shouldBe 2
    consensus[Int](ConsensusTags.PerRead.MinRawReadCount)  shouldBe 2
    consensus[Float](ConsensusTags.PerRead.RawReadErrorRate) shouldBe 0.toFloat
    consensus[String](DefaultTag) shouldBe "AAA"
    consensus[String](ConsensusTags.UmiBases) shouldBe "GAT-ACA"
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

  "VanillaUmiConsensusCaller.filterToMostCommonAlignment" should "return all reads when all cigars are 50M" in {
    val srcs = (1 to 10).map { i => src(cigar="50M") }
    val recs = cc().filterToMostCommonAlignment(srcs)
    recs should have size 10
  }

  it should "return all reads when all cigars are complicated but the same" in {
    val srcs = (1 to 10).map { i => src(cigar="10M5D10M5I20M5S") }
    val recs = cc().filterToMostCommonAlignment(srcs)
    recs should have size 10
  }

  it should "return only the 50M reads (i.e. the most common alignment)" in {
    val buffer = new ArrayBuffer[SourceRead]
    (1 to  3).foreach { i => buffer += src(cigar="25M1D25M") }
    (1 to 10).foreach { i => buffer += src(cigar="50M") }
    (1 to  3).foreach { i => buffer += src(cigar="25M2I23M") }

    val recs = cc().filterToMostCommonAlignment(buffer)
    recs should have size 10
    recs.map(_.cigar.toString()).distinct shouldBe Seq("50M")
  }

  it should "return only the reads with a single base deletion at base 25" in {
    val buffer = new ArrayBuffer[SourceRead]
    // These should all be returned
    (1 to  5).foreach { i => buffer += src(cigar="25M1D25M") }
    (1 to  2).foreach { i => buffer += src(cigar="5S20M1D25M") }
    (1 to  2).foreach { i => buffer += src(cigar="5S20M1D20M5H") }
    (1 to  2).foreach { i => buffer += src(cigar="25M1D20M5S") }

    // These should not be!
    (1 to  2).foreach { i => buffer += src(cigar="25M2D25M") }
    (1 to  2).foreach { i => buffer += src(cigar="25M1I24M") }
    (1 to  2).foreach { i => buffer += src(cigar="20M1D5M1D25M") }

    val recs = cc().filterToMostCommonAlignment(buffer)
    recs should have size 11
    recs.map(_.cigar.toString).distinct.sorted shouldBe Seq("25M1D25M", "5S20M1D25M", "5S20M1D20M5H", "25M1D20M5S").sorted
  }

  it should "return all the reads that are compatible with a 2 base deletion at base 25" in {
    val expected = Seq("25M2D75M", "25M2D65M", "25M2D50M5S", "25M", "24M", "10M").sorted
    val others   = Seq("30M", "25M1D25M", "25M4D25M").sorted
    val input    = (expected ++ others).map(c => src(cigar=c))
    val output   = cc().filterToMostCommonAlignment(input)
    output.map(_.cigar.toString()).sorted shouldBe expected
  }

  it should "return a single read if a single read was given" in {
    val srcs = Seq(src(cigar="50M"))
    val recs = cc().filterToMostCommonAlignment(srcs)
    recs should have size 1
  }

  "VanillaConsensusCaller.toSourceRead" should "mask bases that are below the quality threshold" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, bases="AAAAAAAAAA").map { r => r.quals = Array[Byte](2,30,19,21,18,20,0,30,2,30); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=20.toByte, trim=false).value

    source.baseString shouldBe "NANANANANA"
    source.quals      shouldBe Array[Byte](2,30,2,21,2,20,2,30,2,30)
  }

  it should "trim the source read when the end is low-quality so that there are no trailing no-calls" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, bases="AAAAAAAAAA").map { r => r.quals = Array[Byte](30,30,30,30,30,30,2,2,2,2); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=20.toByte, trim=false).value

    source.baseString shouldBe "AAAAAA"
    source.quals      shouldBe Array[Byte](30,30,30,30,30,30)
  }

  it should "trim the source read when the end of the raw read is all Ns so that there are no trailing no-calls" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, bases="AAAAAANNNN").map { r => r.quals = Array[Byte](30,30,30,30,30,30,30,30,30,30); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=20.toByte, trim=false).value

    source.baseString shouldBe "AAAAAA"
    source.quals      shouldBe Array[Byte](30,30,30,30,30,30)
  }

  it should "trim the source read when the end of the raw read is all Ns and the read is mapped to the negative strand" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, strand=Minus, cigar="4S1M1D5M", bases="NNNNAAAAAA").map { r => r.quals = Array[Byte](30,30,30,30,30,30,30,30,30,30); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=20.toByte, trim=false).value

    source.baseString shouldBe "TTTTTT" // cos revcomp'd
    source.quals      shouldBe Array[Byte](30,30,30,30,30,30)
    source.cigar.toString() shouldBe "5M1D1M"
  }

  it should "trim the source read when the read length is shorter than the insert size" in {
    val builder     = new SamBuilder(readLength=50)
    val Seq(r1, r2) = builder.addPair(start1=11, start2=1, strand1=Plus, strand2=Minus).map { r => r.bases = "A"*10 + "C"*30 + "G"*10; r }
    val s1          = cc().toSourceRead(r1, minBaseQuality=2.toByte, trim=false).value
    val s2          = cc().toSourceRead(r2, minBaseQuality=2.toByte, trim=false).value

    s1.baseString shouldBe "A"*10 + "C"*30 // the first ten bases should be trimmed
    s1.cigar.toString shouldBe "40M"
    s2.baseString shouldBe "C"*10 + "G"*30 // the last ten bases should be trimmed
    s2.cigar.toString shouldBe "40M"
  }

  it should "trim based on insert size in the presence of soft-clipping" in {
    val builder     = new SamBuilder(readLength=50)
    val Seq(r1, r2) = builder.addPair(start1=20, start2=20, strand1=Plus, strand2=Minus, cigar1="10S35M5S", cigar2="12S30M8S")
      .map { r => r.bases = "A"*2 + "C"*46 + "G"*2; r }
    val s1          = cc().toSourceRead(r1, minBaseQuality=2.toByte, trim=false).value
    val s2          = cc().toSourceRead(r2, minBaseQuality=2.toByte, trim=false).value

    s1.baseString shouldBe "A"*2 + "C"*38 // trimmed by ten bases at the 3' end
    s1.cigar.toString shouldBe "10S30M"
    s2.baseString shouldBe "C"*2 + "G"*36 // trimmed by twelve bases at the 3' end
    s2.cigar.toString shouldBe "8S30M" // NB: cigar is reversed in toSourceRead
  }

  it should "not to trim based on insert size in the presence of soft-clipping" in {
    val builder     = new SamBuilder(readLength=142)
    val Seq(r1, r2) = builder.addPair(start1=545, start2=493, strand1=Minus, strand2=Plus, cigar1="47S72M23S", cigar2="46S96M")
      .map { r => r.bases = "A"*142; r }
    val s1          = cc().toSourceRead(r1, minBaseQuality=2.toByte, trim=false).value
    val s2          = cc().toSourceRead(r2, minBaseQuality=2.toByte, trim=false).value

    s1.baseString shouldBe "T"*142
    s1.cigar.toString shouldBe "23S72M47S"  // NB: cigar is reversed in toSourceRead
    s2.baseString shouldBe "A"*142
    s2.cigar.toString shouldBe "46S96M"
  }

  it should "not trim the source read based on insert size if the read is not an FR pair" in {
    val builder     = new SamBuilder(readLength=50)
    val Seq(r1, r2) = builder.addPair(start1=11, start2=1, strand1=Plus, strand2=Plus).map { r => r.bases = "A"*50; r }
    val s1          = cc().toSourceRead(r1, minBaseQuality=2.toByte, trim=false).value
    val s2          = cc().toSourceRead(r2, minBaseQuality=2.toByte, trim=false).value

    s1.baseString shouldBe "A"*50
    s1.cigar.toString shouldBe "50M"
    s2.baseString shouldBe "A"*50
    s2.cigar.toString shouldBe "50M"
  }

  it should "return None if the read is all low-quality or Ns" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, bases="NANANANANA").map { r => r.quals = Array[Byte](30,2,30,2,30,2,30,2,30,2); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=20.toByte, trim=false)
    source shouldBe None
  }

  it should "apply phred-style quality trimming to the read in addition to masking" in {
    val builder = new SamBuilder(readLength=10)
    val rec     = builder.addFrag(start=1, bases="AGCACGACGT").map { r => r.quals = Array[Byte](30,30,30,2,5,2,3,20,2,6); r}.value
    val source  = cc().toSourceRead(rec, minBaseQuality=15.toByte, trim=true).value
    source.baseString shouldBe "AGC"
    source.quals should have length 3
    source.cigar.toString() shouldBe "3M"
    source.sam shouldBe Some(rec)
  }
}
