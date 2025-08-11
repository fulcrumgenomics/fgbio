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
import com.fulcrumgenomics.util.Io

import java.nio.file.Paths

class CallCodecConsensusReadsTest extends UnitSpec {
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

  it should "emit statistics and a reject BAM" in {
    for (threads <- Seq(1, 4)) {
      val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.TemplateCoordinate))
      builder.addPair( // Should form a consensus
        start1=100, start2=100, bases1="AC" * 15, bases2="AC" * 15, attrs=Map(("RX", "ACC-TGA"), ("MI", "hi"))
      )
      builder.addPair( // Too far apart to form a consensus
        name="x", start1=200, start2=500, bases1="AC" * 15, bases2="AC" * 15, attrs=Map(("RX", "ACC-TGA"), ("MI", "bye"))
      )
      val in  = builder.toTempFile()

      val out = makeTempFile("codec.", ".bam")
      val rej = makeTempFile("rejects.", ".bam")
      val statsPath = makeTempFile("stats.", ".txt")
      val caller = new CallCodecConsensusReads(input=in, output=out, readGroupId="ZZ", threads=threads, rejects=Some(rej), stats=Some(statsPath))
      caller.execute()

      val recs = readBamRecs(out)
      recs should have size 1

      val rejectedRecs = readBamRecs(rej)
      rejectedRecs should have size 2
      rejectedRecs.foreach(_.name shouldBe "x")

      val stats: Map[String, String] = Io.readLines(statsPath).drop(1)
        .map { line =>
          val tabIdx = line.indexOf('\t')
          (line.substring(0, tabIdx), line.substring(tabIdx + 1))
        }
        .toMap

      val key = stats.keys.find(_.contains("overlap too short")).getOrElse(fail("Couldn't find key in stats."))
      val count = stats(key).toInt
      count shouldBe 2
    }
  }
}
