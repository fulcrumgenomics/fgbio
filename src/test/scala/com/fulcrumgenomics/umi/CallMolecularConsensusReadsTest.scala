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
import com.fulcrumgenomics.umi.UmiConsensusCaller.{ConsensusKvMetric, RejectionReason}
import com.fulcrumgenomics.util.Metric
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
    val specialCellTag = "XX"
    val rlen    = 100
    val builder = new SamBuilder(baseQuality=30, readLength=rlen, readGroupId=Some("ABC"), sort=Some(SamOrder.TemplateCoordinate))
    val output  = newBam
    val rejects = newBam

    builder.addFrag(name="a1", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a", specialCellTag -> "AB"))
    builder.addFrag(name="a2", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a", specialCellTag -> "AB"))
    builder.addFrag(name="a3", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACGT", "MI" -> "a", specialCellTag -> "AB"))

    builder.addFrag(name="b1", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACAC", "MI" -> "b", specialCellTag -> "AB"))
    builder.addFrag(name="b2", start=100, bases="A"*rlen, attrs=Map("RX" -> "ACAC", "MI" -> "b", specialCellTag -> "AB"))

    // Run the tool
    new CallMolecularConsensusReads(
      input       = builder.toTempFile(),
      output      = output,
      minReads    = 1,
      cellTag     = Some(specialCellTag),
      rejects     = Some(rejects),
      readGroupId = "ABC"
    ).execute()

    // check we have no rejected records
    readBamRecs(rejects).isEmpty shouldBe true

    // we should have 2 consensus paired end reads
    val records = readBamRecs(output)
    records.size shouldBe 2
    records.count { rec => !rec.paired } shouldBe 2

    records.foreach { rec =>
      rec.readGroup.getId shouldBe "ABC"
      rec.basesString shouldBe "A" * 100
      rec.length shouldBe 100
      rec[String](specialCellTag) shouldBe "AB"
    }
  }

  /** Builds a tag family of three read pairs at the same coordinates, two with both ends mapped and one whose R2 is
    * unmapped.
    *
    * The two mapped R2s carry `G`s while the unmapped R2 carries `T`s, so the consensus R2 shows which reads
    * contributed to it.  R2s are placed on the negative strand by [[SamBuilder]], so `toSourceRead` reverse
    * complements them and a consensus built from the mapped R2s alone reads as `C`s. */
  private def halfMappedFamily(rlen: Int): SamBuilder = {
    val builder = new SamBuilder(baseQuality=30, readLength=rlen, readGroupId=Some("ABC"), sort=Some(SamOrder.TemplateCoordinate))
    val attrs   = Map(DefaultTag -> "GATTACA:1")
    builder.addPair(name="mapped:1",     start1=100, start2=300, bases1="A"*rlen, bases2="G"*rlen, attrs=attrs)
    builder.addPair(name="mapped:2",     start1=100, start2=300, bases1="A"*rlen, bases2="G"*rlen, attrs=attrs)
    builder.addPair(name="halfmapped:3", start1=100, start2=300, bases1="A"*rlen, bases2="T"*rlen, unmapped2=true, attrs=attrs)
    builder
  }

  Seq(2, 3).foreach { maxReads =>
    it should f"cap R1 at max-reads and build R2 from the mapped reads only with --max-reads=$maxReads" in {
      val rlen    = 100
      val output  = newBam
      val rejects = newBam
      val stats   = makeTempFile("call_molecular_consensus_reads_test.", ".txt")

      new CallMolecularConsensusReads(
        input       = halfMappedFamily(rlen).toTempFile(),
        output      = output,
        minReads    = 1,
        maxReads    = Some(maxReads),
        rejects     = Some(rejects),
        stats       = Some(stats),
        readGroupId = "ABC"
      ).execute()

      // Exactly one consensus read pair is produced.
      val records = readBamRecs(output)
      records.size shouldBe 2
      val r1 = records.find(_.firstOfPair).value
      val r2 = records.find(_.secondOfPair).value

      // All three R1s are mapped and eligible, so R1 is capped at max-reads.
      r1[Int](ConsensusTags.PerRead.RawReadCount) shouldBe math.min(3, maxReads)
      r1.basesString shouldBe "A" * rlen

      // Only the two mapped R2s may contribute, whether or not the cap binds.  The unmapped R2's `C`s must not reach
      // the consensus, and it must not be counted towards the consensus depth.
      r2[Int](ConsensusTags.PerRead.RawReadCount) shouldBe 2
      r2.basesString shouldBe "C" * rlen

      // The unmapped R2 must not have been used, and must be accounted for rather than vanishing:
      // it is written to the rejects BAM and counted under its own rejection reason.
      val rejected = readBamRecs(rejects)
      rejected.map(rec => (rec.name, rec.secondOfPair)) should contain theSameElementsAs Seq(("halfmapped:3", true))

      val metrics = Metric.read[ConsensusKvMetric](stats).map(m => m.key -> m.value.toString).toMap
      metrics(s"raw_reads_rejected_for_${RejectionReason.Unmapped.code}") shouldBe "1"
      metrics("raw_reads_considered") shouldBe "6"
    }
  }
}
