/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.{SamOrder, SamSource}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.UmiConsensusCaller.{ConsensusKvMetric, RejectionReason}
import com.fulcrumgenomics.util.Metric
import org.scalatest.OptionValues

import java.nio.file.Paths

class CallCodecConsensusReadsTest extends UnitSpec with OptionValues {
  "CallCodecConsensusReads" should "throw an exception if the input file doesn't exist" in {
    an[Throwable] should be thrownBy {
      new CallCodecConsensusReads(input=Paths.get("/tmp/path/to/no/where/foo.bam"), output=Paths.get("/tmp")).execute()
    }
  }

  it should "throw an exception if the output file isn't writable" in {
    an[Throwable] should be thrownBy {
      val in = makeTempFile("in.", ".bam")
      val out = Paths.get("/tmp/path/to/no/where.bam")
      new CallCodecConsensusReads(input=in, output=out).execute()
    }
  }

  it should "throw an exception if numeric parameters are set incorrectly" in {
    val in  = makeTempFile("in.", ".bam")
    val out = makeTempFile("out.", ".bam")
    an[Exception] should be thrownBy { new CallCodecConsensusReads(input=in, output=out, errorRatePreUmi=0.toByte).execute() }
    an[Exception] should be thrownBy { new CallCodecConsensusReads(input=in, output=out, errorRatePostUmi=0.toByte).execute() }
    an[Exception] should be thrownBy { new CallCodecConsensusReads(input=in, output=out, minReadPairs=0).execute() }
    an[Exception] should be thrownBy { new CallCodecConsensusReads(input=in, output=out, minReadPairs=2, maxReadPairs=Some(1)).execute() }
    an[Exception] should be thrownBy { new CallCodecConsensusReads(input=in, output=out, minDuplexLength=0).execute() }
  }

  it should "have working CLP and arg annotations" in {
    checkClpAnnotations[CallCodecConsensusReads]
  }

  it should "run correctly on an empty input BAM" in {
    val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.TemplateCoordinate))
    val in  = builder.toTempFile()
    val out = makeTempFile("codec.", ".bam")
    new CallCodecConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq
    recs shouldBe Seq()
  }

  it should "run correctly on a BAM with some reads in it" in {
    val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(
      start1=100, start2=100, bases1="AC" * 15, bases2="AC" * 15, attrs=Map(("RX", "ACC-TGA"), ("MI", "hi"))
    )
    val in  = builder.toTempFile()

    val out = makeTempFile("codec.", ".bam")
    new CallCodecConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq
    recs should have size 1
  }

  it should "filter on supported duplex bases rather than reference overlap through the CLI path" in {
    for (minDuplexLength <- Seq(30, 40, 41)) {
      val builder = new SamBuilder(readLength=40, baseQuality=35, sort=Some(SamOrder.TemplateCoordinate))
      builder.addPair(start1=111, start2=101, cigar1="10S30M", cigar2="30M10S",
        bases1="AC" * 20, bases2="AC" * 20, attrs=Map(("RX", "ACC-TGA"), ("MI", "hi")))
      val in = builder.toTempFile()
      val out = makeTempFile("codec.", ".bam")
      val rej = makeTempFile("rejects.", ".bam")
      val stats = makeTempFile("stats.", ".txt")
      new CallCodecConsensusReads(input=in, output=out, rejects=Some(rej), stats=Some(stats),
        minDuplexLength=minDuplexLength).execute()

      val consensus = readBamRecs(out)
      val rejected = readBamRecs(rej)
      val metrics = Metric.read[ConsensusKvMetric](stats).map(m => m.key -> m.value).toMap
      if (minDuplexLength <= 40) {
        consensus should have size 1
        consensus.head.basesString shouldBe "AC" * 20
        rejected shouldBe empty
        metrics("consensus_duplex_bases_emitted").toString.toInt shouldBe 40
      }
      else {
        consensus shouldBe empty
        rejected should have size 2
        rejected.foreach { rec =>
          rec[String](UmiConsensusCaller.RejectReasonTag) shouldBe RejectionReason.R1R2OverlapTooShort.code
        }
        metrics("raw_reads_rejected_for_r1_r2_overlap_too_short").toString.toInt shouldBe 2
      }
    }
  }

  it should "emit statistics and a reject BAM" in {
    val specialCellTag = "XX"
    for (threads <- Seq(1, 4)) {
      val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.TemplateCoordinate))
      builder.addPair( // Should form a consensus
        start1=100, start2=100, bases1="AC" * 15, bases2="AC" * 15, attrs=Map(("RX", "ACC-TGA"), ("MI", "hi"), (specialCellTag, "AB"))
      )
      builder.addPair( // Too far apart to form a consensus
        name="x", start1=200, start2=500, bases1="AC" * 15, bases2="AC" * 15, attrs=Map(("RX", "ACC-TGA"), ("MI", "bye"), (specialCellTag, "AB"))
      )
      val in  = builder.toTempFile()

      val out = makeTempFile("codec.", ".bam")
      val rej = makeTempFile("rejects.", ".bam")
      val statsPath = makeTempFile("stats.", ".txt")
      val caller = new CallCodecConsensusReads(
        input       = in,
        output      = out,
        readGroupId ="ZZ",
        threads     = threads,
        cellTag     = Some(specialCellTag),
        rejects     = Some(rej),
        stats       = Some(statsPath)
      )
      caller.execute()

      val recs = readBamRecs(out)
      recs should have size 1
      recs.foreach { rec =>
        rec[String](specialCellTag) shouldBe "AB"
      }

      val rejectedRecs = readBamRecs(rej)
      rejectedRecs should have size 2
      rejectedRecs.foreach(_.name shouldBe "x")

      val metrics = Metric.read[ConsensusKvMetric](statsPath).toIndexedSeq
      val metric = metrics.find(_.key.contains(RejectionReason.R1R2OverlapTooShort.code)).value
      metric.value.toString.toInt shouldBe 2
    }
  }
}
