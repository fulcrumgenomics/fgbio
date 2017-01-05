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

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.Io
import dagr.commons.io.PathUtil
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.{SAMFileHeader, SAMFileWriterFactory, SAMRecord, SAMSequenceDictionary, SAMSequenceRecord, SamReaderFactory}

import scala.collection.JavaConversions.iterableAsScalaIterable
import scala.collection.mutable.ListBuffer

class FilterConsensusReadsTest extends UnitSpec {
  private val temp = makeTempFile("meh.", ".bam")

  /** Make a reference file that is 100 lines of 100 As. */
  private val ref = {
    val fa   = makeTempFile("ref.", ".fa")
    val dict = PathUtil.replaceExtension(fa, ".dict")

    val lines = new ListBuffer[String]()
    lines += ">chr1"
    Range.inclusive(1, 100).foreach(i => lines += ("A" * 100))
    Io.writeLines(fa, lines)

    val header = new SAMFileHeader
    header.setSequenceDictionary(new SAMSequenceDictionary(util.Arrays.asList(new SAMSequenceRecord("chr1", 100*100))))
    new SAMFileWriterFactory().makeSAMWriter(header, true, dict.toFile).close()
    dict.toFile.deleteOnExit()

    fa
  }

  /** Generates a FilterConsensusReads instance for testing non-end-to-end methods. */
  private def f(q: Int, d: Int, readErr: Double, baseErr: Double, nf: Double) : FilterConsensusReads =
    new FilterConsensusReads(input=temp, output=temp, ref=temp, reversePerBaseTags=false,
      minBaseQuality=q.toByte, minReads=d, maxReadErrorRate=readErr, maxBaseErrorRate=baseErr, maxNoCallFraction=nf)

  /** Tags up a SAMRecord for filtering. */
  private def tag(rec: SAMRecord, minDepth: Int, depth: Int, readErr: Float, depths: Array[Short]=null, errors: Array[Short]=null): SAMRecord = {
    rec.setReadString("A" * rec.getReadLength)
    rec.setAttribute(ConsensusTags.PerRead.RawReadCount, depth)
    rec.setAttribute(ConsensusTags.PerRead.MinRawReadCount, minDepth)
    rec.setAttribute(ConsensusTags.PerRead.RawReadErrorRate, readErr)
    rec.setAttribute(ConsensusTags.PerBase.RawReadCount, depths)
    rec.setAttribute(ConsensusTags.PerBase.RawReadErrors, errors)
    rec
  }

  /** Makes a SAMRecord with which to test filtering. */
  private def r(len: Int=10, q: Int=40, minDepth: Int, depth: Int, readErr: Float, depths: Array[Short]=null, errors: Array[Short]=null) : SAMRecord = {
    val builder = new SamRecordSetBuilder(readLength=len, baseQuality=q.toByte)
    tag(builder.addFrag(start=100).get, minDepth=minDepth, depth=depth, readErr=readErr, depths=depths, errors=errors)
  }

  /** Makes an array of Shorts of all the same value. */
  private def arr(value: Short, length: Int): Array[Short] = {
    val ss = new Array[Short](length)
    forloop (from=0, until=length) { i => ss(i) = value }
    ss
  }

  "FilterConsensusReads.filterRead" should "keep rec and not mask any bases" in {
    val filter = f(q=10, d=2, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=2, minDepth=2, readErr = 0.01f)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 0
  }

  it should "keep rec and not mask any bases when detail tags are present" in {
    val filter = f(q=10, d=2, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=2, minDepth=2, readErr=0.01f, depths=arr(2, 10), errors=arr(0, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 0
    rec.getReadString shouldBe ("A" * 10)
  }

  it should "filter out a read with depth < minDepth" in {
    val filter = f(q=10, d=5, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=4, minDepth=4, readErr=0.00f, depths=arr(4, 10), errors=arr(0, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.getReadString shouldBe ("A" * 10)
  }

  it should "filter out a read with too high an error rate" in {
    val filter = f(q=10, d=3, readErr=0.05, baseErr=0.5, nf=0.2)
    val rec    = r(depth=50, minDepth=50, readErr=0.20f, depths=arr(50, 10), errors=arr(10, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.getReadString shouldBe ("A" * 10)
  }

  it should "filter out bases that are below the required quality score" in {
    val filter = f(q=40, d=3, readErr=0.05, baseErr=0.1, nf=0.5)
    val rec    = r(depth=5, minDepth=5, readErr=0.00f, depths=arr(5, 10), errors=arr(0, 10))
    Seq(1,4,7).foreach(i => rec.getBaseQualities()(i) = 20)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 3
    rec.getReadString shouldBe "ANAANAANAA"
  }

  it should "filter out bases that have too much disagreement (errors/depth > baseErr)" in {
    val filter = f(q=40, d=10, readErr=0.05, baseErr=0.1, nf=0.5)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f, depths=arr(20, 10), errors=Array[Short](0,1,2,3,4,3,2,1,0,0))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 3
    rec.getReadString shouldBe "AAANNNAAAA"
  }

  it should "filter out the read if too many bases are masked in the read in the first place" in {
    val filter = f(q=20, d=5, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f)
    rec.setReadString("AANNNNNAAA")
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.getReadString shouldBe "AANNNNNAAA"
  }

  it should "filter out the read if too many bases are masked post-filtering" in {
    val filter = f(q=20, d=5, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f, depths=arr(20, 10), errors=Array[Short](3,3,3,0,0,0,0,3,3,3))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 6
    rec.getReadString shouldBe "NNNAAAANNN"
  }

  "FilterConsensusReads" should "not filter out any reads" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=45, sortOrder=SortOrder.queryname)
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q2", start1=100, start2=200).foreach(r => tag(r, minDepth=5, depth=5, readErr=0f, depths=arr(5, 10), errors=arr(0,10)))
    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=3, maxReadErrorRate=0.025, maxBaseErrorRate=0.1, maxNoCallFraction=0.1).execute()

    val recs = SamReaderFactory.make().open(out).toSeq
    recs.size shouldBe 4
    recs.exists(_.getReadString.contains("N")) shouldBe false
  }

  it should "filter out both reads of a pair whenever one fails" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=45, sortOrder=SortOrder.queryname)
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q2", start1=100, start2=200).foreach(r => tag(r, minDepth=5, depth=5, readErr=0f, depths=arr(5, 10), errors=arr(0,10)))
    builder.filter(_.getReadName == "q1").filter(_.getFirstOfPairFlag).foreach(_.setAttribute(ConsensusTags.PerRead.RawReadErrorRate, 0.2f))

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=3, maxReadErrorRate=0.025, maxBaseErrorRate=0.1, maxNoCallFraction=0.1).execute()

    val recs = SamReaderFactory.make().open(out).toSeq
    recs.size shouldBe 2
    recs.map(_.getReadName).distinct shouldBe Seq("q2")
  }

  it should "filter out all reads with the same query name when any primary read fails" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=45, sortOrder=SortOrder.queryname)
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q1", start1=100, start2=200).foreach { r =>
      tag(r, minDepth = 5, depth = 5, readErr = 0f, depths = arr(5, 10), errors = arr(0, 10))
      r.setSupplementaryAlignmentFlag(true)
    }

    builder.filterNot(_.getSupplementaryAlignmentFlag).filter(_.getFirstOfPairFlag).foreach(_.setAttribute(ConsensusTags.PerRead.RawReadErrorRate, 0.2f))

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=3, maxReadErrorRate=0.025, maxBaseErrorRate=0.1, maxNoCallFraction=0.1).execute()

    val recs = SamReaderFactory.make().open(out).toSeq
    recs.size shouldBe 0
  }

  it should "not filter out primary alignments or other supplementals when a supplemental or secondary fails" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=45, sortOrder=SortOrder.unsorted)
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q1", start1=110, start2=190).foreach { r =>
      tag(r, minDepth = 5, depth = 5, readErr = 0f, depths = arr(5, 10), errors = arr(0, 10))
      r.setSupplementaryAlignmentFlag(true)
    }

    builder.filter(_.getSupplementaryAlignmentFlag).filter(_.getFirstOfPairFlag).foreach(_.setAttribute(ConsensusTags.PerRead.RawReadErrorRate, 0.2f))

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=3, maxReadErrorRate=0.025, maxBaseErrorRate=0.1, maxNoCallFraction=0.1).execute()

    val recs = SamReaderFactory.make().open(out).toSeq
    recs.size shouldBe 3
    recs.count(!_.isSecondaryOrSupplementary) shouldBe 2
    recs.count(_.isSecondaryOrSupplementary) shouldBe 1
  }

  it should "reverse the per-base tags only when requested" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=45, sortOrder=SortOrder.unsorted)
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=7, depth=9, readErr=0f, depths=Array[Short](7,7,7,8,8,8,9,9,9,9), errors=arr(0,10)))

    Seq(true, false) foreach { reverseTags =>
      val out = makeTempFile("filtered.", ".bam")
      new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=reverseTags,
        minBaseQuality=45.toByte, minReads=3, maxReadErrorRate=0.025, maxBaseErrorRate=0.1, maxNoCallFraction=0.1).execute()

      val recs = SamReaderFactory.make().open(out).toSeq
      recs.size shouldBe 2
      recs.exists(_.getReadNegativeStrandFlag) shouldBe true
      recs.foreach { rec =>
        val depths = rec.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadCount)
        if (reverseTags && rec.getReadNegativeStrandFlag) depths shouldBe Array[Short](9,9,9,9,8,8,8,7,7,7) else depths shouldBe Array[Short](7,7,7,8,8,8,9,9,9,9)
      }
    }
  }
}
