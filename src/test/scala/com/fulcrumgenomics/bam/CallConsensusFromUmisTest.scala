/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.bam.ConsensusCallerOptions.DefaultAttribute
import htsjdk.samtools.util.CloserUtil
import htsjdk.samtools.{SAMRecordSetBuilder, SAMUtils}
import com.fulcrumgenomics.util.PhredValue.ZeroProbability

import scala.collection.JavaConversions._
import scala.collection.JavaConverters._

/**
  * Tests for CallConsensusFromUmis.
  */
class CallConsensusFromUmisTest extends UnitSpec {

  // There be dragons below!
  "CallConsensusFromUmis" should "should create two consensus for two UMI groups" in {
    val builder = new SAMRecordSetBuilder()
    builder.addFrag("READ1", 0, 1, false).setAttribute(DefaultAttribute, "GATTACA")
    builder.addFrag("READ2", 0, 1, false).setAttribute(DefaultAttribute, "GATTACA")
    builder.addFrag("READ3", 0, 1, false).setAttribute(DefaultAttribute, "ACATTAG")
    builder.addFrag("READ4", 0, 1, false).setAttribute(DefaultAttribute, "ACATTAG")
    builder.getRecords.foreach { rec =>
      rec.setReadString("A" * rec.getReadLength)
      rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new ConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new ConsensusCallerOptions(
        minReads=1,
        errorRatePreUmi=ZeroProbability,
        errorRatePostUmi=ZeroProbability
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
        rec.setAttribute(DefaultAttribute, "GATTACA")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new ConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new ConsensusCallerOptions(
        minReads=1,
        errorRatePreUmi=ZeroProbability,
        errorRatePostUmi=ZeroProbability
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
        rec.setAttribute(DefaultAttribute, "GATTACA")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    builder.getRecords.slice(2, 4).foreach {
      rec =>
        rec.setAttribute(DefaultAttribute, "ACATTAG")
        rec.setReadString("A" * rec.getReadLength)
        rec.setBaseQualityString(SAMUtils.phredToFastq(40).toString * rec.getReadLength)
    }
    val reader = builder.getSamReader
    val consensusCaller = new ConsensusCaller(
      input = reader.iterator().asScala,
      header = reader.getFileHeader,
      options = new ConsensusCallerOptions(
        minReads=1,
        errorRatePreUmi=ZeroProbability,
        errorRatePostUmi=ZeroProbability
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
