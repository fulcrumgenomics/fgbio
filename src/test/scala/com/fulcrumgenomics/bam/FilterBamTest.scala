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

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.util.{Interval, IntervalList}

class FilterBamTest extends UnitSpec {
  def newBam = makeTempFile("filter_bam_test.", ".bam")

  // Make a set of records to run multiple tests on
  val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
  builder.addPair("ok_1", 0, 100, 300)
  builder.addPair("ok_2", 0, 200, 400)
  builder.addPair("ok_3", 0, 300, 500)
  builder.addPair("ok_4", 0, 400, 600)
  builder.addPair("unmapped_1",  0, 500, 600, record1Unmapped=true, record2Unmapped=true)
  builder.addPair("secondary_1", 0, 600, 800).foreach(_.setNotPrimaryAlignmentFlag(true))
  builder.addPair("duplicate_1", 0, 700, 900).foreach(_.setDuplicateReadFlag(true))
  builder.addPair("lowmapq_1"  , 0, 900, 1100, mapq1=5, mapq2=5)
  builder.addPair("ok_ontarget_1", 1, 5000, 5500)
  builder.addPair("ok_ontarget_2", 1, 5200, 5600)
  builder.addPair("ok_ontarget_3", 1, 5500, 5800)

  val intervals = new IntervalList(builder.header)
  intervals.add(new Interval(builder.header.getSequence(1).getSequenceName, 5000, 6000))

  "FilterBam" should "run ok on an empty bam" in {
    val empty = new SamRecordSetBuilder().toTempFile()
    val out   = newBam
    new FilterBam(input=empty, output=out).execute()
    val reader = SamReaderFactory.make.open(out.toFile)
    reader.iterator().hasNext shouldBe false
  }

  it should "not remove any records" in {
    val out  = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeDuplicates=false, removeUnmappedReads=false, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.toSeq.map(_.getSAMString)
  }

  it should "remove all but the 'ok' records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, minMapQ=10).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filter(_.getReadName.startsWith("ok_")).toSeq.map(_.getSAMString)
  }

  it should "remove just the unmapped records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=true, removeDuplicates=false, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filterNot(_.getReadName.startsWith("unmapped_")).toSeq.map(_.getSAMString)
  }

  it should "remove just the secondary records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=false, minMapQ=0, removeSecondaryAlignments=true).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filterNot(_.getReadName.startsWith("secondary_")).toSeq.map(_.getSAMString)
  }

  it should "remove just the duplicate records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=true, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filterNot(_.getReadName.startsWith("duplicate_")).toSeq.map(_.getSAMString)
  }

  it should "remove just the low mapq records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=false, minMapQ=20, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filterNot(_.getReadName.startsWith("lowmapq_")).toSeq.map(_.getSAMString)
  }

  it should "only return matching reads within the intervals" in {
    val out = newBam
    val ilist = makeTempFile("filter_bam_test.", ".interval_list")
    intervals.write(ilist.toFile)

    new FilterBam(input=builder.toTempFile(), output=out, minMapQ=20, intervals=Some(ilist)).execute()
    val recs = readBamRecs(out)
    recs.map(_.getSAMString) should contain theSameElementsAs builder.filter(_.getReadName.startsWith("ok_ontarget")).toSeq.map(_.getSAMString)
  }
}
