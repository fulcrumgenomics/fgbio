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
import com.fulcrumgenomics.umi.ConsensusCallerOptions._
import htsjdk.samtools.SAMFileHeader.SortOrder

/**
  * Tests for CallMolecularConsensusReads.
  *
  * This makes sure the tool runs end-to-end, and the majority of the tests that cover various options are covered
  * in [[ConsensusCallerTest]].
  */
class CallMolecularConsensusReadsTest extends UnitSpec {

  def newBam = makeTempFile("call_molecular_consensus_reads_test.", ".bam")

  "CallMolecularConsensusReads" should "run end-to-end" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.unsorted, baseQuality=30, readLength=100, readGroupId=Some("ABC"))
    val output  = newBam
    val rejects = newBam

    // Create 2000 paired end reads, where there are two pairs with the same coordinates and have the same group tag.
    Stream.range(0, 1000).foreach { idx =>
      val firstPair  = builder.addPair(s"READ:" + 2*idx,   0, 1+idx, 1000000+idx)
      val secondPair = builder.addPair(s"READ:" + 2*idx+1, 0, 1+idx, 1000000+idx)
      // set the read and add the unique molecule tag.
      (firstPair ++ secondPair).foreach { rec =>
        rec.setAttribute(DefaultTag, "GATTACA:" + idx)
        rec.setReadString("A" * rec.getReadLength)
      }
    }

    // Run the tool
    new CallMolecularConsensusReads(input=builder.toTempFile(), output=output, rejects=rejects, readGroupId="ABC").execute()

    // check we have no rejected records
    readBamRecs(rejects) shouldBe 'empty

    // we should have 1000 consensus paired end reads
    val records = readBamRecs(output)
    records.count { rec => rec.getFirstOfPairFlag } shouldBe 1000
    records.count { rec => rec.getSecondOfPairFlag } shouldBe 1000
    records.foreach { rec =>
      rec.getReadGroup.getId shouldBe "ABC"
      rec.getReadString shouldBe "A" * rec.getReadLength
      rec.getReadLength shouldBe 100
    }
  }
}
