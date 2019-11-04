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

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._

/**
  * Tests for CallMolecularConsensusReads.
  *
  * This makes sure the tool runs end-to-end, and the majority of the tests that cover various options are covered
  * in [[VanillaUmiConsensusCallerTest]].
  */
class CallMolecularConsensusReadsTest extends UnitSpec {

  private def newBam = makeTempFile("call_molecular_consensus_reads_test.", ".bam")

  "CallMolecularConsensusReads" should "run end-to-end" in {
    val rlen    = 100
    val builder = new SamBuilder(baseQuality=30, readLength=rlen, readGroupId=Some("ABC"), sort=Some(SamOrder.TemplateCoordinate))
    val output  = newBam
    val rejects = newBam

    // Create 2000 paired end reads, where there are two pairs with the same coordinates and have the same group tag.
    Range(0, 1000).foreach { idx =>
      val attrs = Map(DefaultTag -> ("GATTACA:" + idx), ConsensusTags.UmiBases -> "ACGT-TGCA")
      builder.addPair(name=s"READ:" + 2*idx,   start1=1+idx, start2=1000000+idx, bases1="A"*rlen, bases2="T"*rlen, attrs=attrs)
      builder.addPair(name=s"READ:" + 2*idx+1, start1=1+idx, start2=1000000+idx, bases1="A"*rlen, bases2="T"*rlen, attrs=attrs)
    }

    // Run the tool
    new CallMolecularConsensusReads(input=builder.toTempFile(), output=output, minReads=1, rejects=Some(rejects), readGroupId="ABC").execute()

    // check we have no rejected records
    readBamRecs(rejects).isEmpty shouldBe true

    // we should have 1000 consensus paired end reads
    val records = readBamRecs(output)
    records.count { rec => rec.firstOfPair } shouldBe 1000
    records.count { rec => rec.secondOfPair } shouldBe 1000
    records.foreach { rec =>
      rec.readGroup.getId shouldBe "ABC"
      rec.basesString shouldBe "A" * 100
      rec.length shouldBe 100
      rec[String](DefaultTag).startsWith("GATTACA") shouldBe true
      rec[String](ConsensusTags.UmiBases) shouldBe "ACGT-TGCA"
    }
  }

  it should "run end-to-end on single-end data" in {
    val rlen    = 100
    val builder = new SamBuilder(baseQuality=30, readLength=rlen, readGroupId=Some("ABC"), sort=Some(SamOrder.TemplateCoordinate))
    val output  = newBam
    val rejects = newBam

    builder.addFrag(name="a1", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a"))
    builder.addFrag(name="a2", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a"))
    builder.addFrag(name="a3", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a"))

    builder.addFrag(name="b1", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACAC", "MI" -> "b"))
    builder.addFrag(name="b2", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACAC", "MI" -> "b"))

    // Run the tool
    new CallMolecularConsensusReads(input=builder.toTempFile(), output=output, minReads=1, rejects=Some(rejects), readGroupId="ABC").execute()

    // check we have no rejected records
    readBamRecs(rejects).isEmpty shouldBe true

    // we should have 1000 consensus paired end reads
    val records = readBamRecs(output)
    records.size shouldBe 2
    records.count { rec => !rec.paired } shouldBe 2

    records.foreach { rec =>
      rec.readGroup.getId shouldBe "ABC"
      rec.basesString shouldBe "A" * 100
      rec.length shouldBe 100
    }
  }
}
