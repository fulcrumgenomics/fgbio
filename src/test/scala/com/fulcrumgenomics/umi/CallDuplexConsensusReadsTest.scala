/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

import java.nio.file.Paths

import com.fulcrumgenomics.bam.api.{SamOrder, SamSource}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class CallDuplexConsensusReadsTest extends UnitSpec {
  private val MI = ConsensusTags.MolecularId
  private val RX = ConsensusTags.UmiBases

  "CallDuplexConsensusReads" should "throw an exception if the input file doesn't exist" in {
    an[Throwable] should be thrownBy {
      new CallDuplexConsensusReads(input=Paths.get("/tmp/path/to/no/where/foo.bam"), output=Paths.get("/tmp")).execute()
    }
  }

  it should "throw an exception if the output file isn't writable" in {
    an[Throwable] should be thrownBy {
      val in = makeTempFile("in.", ".bam")
      val out = Paths.get("/tmp/path/to/no/where.bam")
      new CallDuplexConsensusReads(input=in, output=out).execute()
    }
  }

  it should "throw an exception if either error rate is set too low" in {
      val in  = makeTempFile("in.", ".bam")
      val out = makeTempFile("out.", ".bam")
    an[Exception] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, errorRatePreUmi=0.toByte).execute() }
    an[Exception] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, errorRatePostUmi=0.toByte).execute() }
  }

  it should "have working CLP and arg annotations" in {
    checkClpAnnotations[CallDuplexConsensusReads]
  }

  it should "not generate a consensus if AB-R1s are not on the same strand ads BA-R2s" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="ab1", start1=100, start2=200, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Plus)
    builder.addPair(name="ba1", start1=200, start2=100, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Minus)

    val in  = builder.toTempFile()
    val out = makeTempFile("duplex.", ".bam")
    new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq

    reader.header.getReadGroups should have size 1
    reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
    recs should have size 0
  }

  it should "not generate a consensus if AB-R2s are not on the same strand ads BA-R1s" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="ab1", start1=100, start2=200, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Minus)
    builder.addPair(name="ba1", start1=200, start2=100, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Plus)

    val in  = builder.toTempFile()
    val out = makeTempFile("duplex.", ".bam")
    new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq

    reader.header.getReadGroups should have size 1
    reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
    recs should have size 0
  }

  it should "run successfully and create consensus reads" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="ab1", start1=100, start2=100, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    builder.addPair(name="ab2", start1=100, start2=100, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    builder.addPair(name="ab3", start1=100, start2=100, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    builder.addPair(name="ba1", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    builder.addPair(name="ba2", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    builder.addPair(name="ba3", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")

    // Add the original UMI bases to each read
    builder.foreach { rec =>
      val mi = rec[String](MI)
      // first of pair ABs and second of pair BAs
      if ((rec.firstOfPair && mi.endsWith("/A")) || (rec.secondOfPair && mi.endsWith("/B"))) {
        rec(RX) = "AAT-CCG"
      }
      else {
        rec(RX) = "CCG-AAT"
      }
    }

    val in  = builder.toTempFile()
    val out = makeTempFile("duplex.", ".bam")
    new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq

    reader.header.getReadGroups should have size 1
    reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
    recs should have size 2
    recs.foreach { rec =>
      rec[String](MI) shouldBe "1"
      rec[String](RX) shouldBe (if (rec.firstOfPair) "AAT-CCG" else "CCG-AAT")
    }
  }
}
