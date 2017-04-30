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

import com.fulcrumgenomics.bam.api.{SamOrder, SamSource}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.util.{Interval, IntervalList}

class FilterBamTest extends UnitSpec {
  def newBam = makeTempFile("filter_bam_test.", ".bam")

  // Make a set of records to run multiple tests on
  val builder = new SamBuilder(sort=Some(SamOrder.Coordinate))
  builder.addPair(name="ok_1", contig=0, start1=100, start2=300)
  builder.addPair(name="ok_2", contig=0, start1=200, start2=400)
  builder.addPair(name="ok_3", contig=0, start1=300, start2=500)
  builder.addPair(name="ok_4", contig=0, start1=400, start2=600)
  builder.addPair(name="unmapped_1",  contig=0, start1=500, start2=600, unmapped1=true, unmapped2=true)
  builder.addPair(name="secondary_1", contig=0, start1=600, start2=800).foreach(_.secondary = true)
  builder.addPair(name="duplicate_1", contig=0, start1=700, start2=900).foreach(_.duplicate = true)
  builder.addPair(name="lowmapq_1"  , contig=0, start1=900, start2=1100, mapq1=5, mapq2=5)
  builder.addPair(name="ok_ontarget_1", contig=1, start1=5000, start2=5500)
  builder.addPair(name="ok_ontarget_2", contig=1, start1=5200, start2=5600)
  builder.addPair(name="ok_ontarget_3", contig=1, start1=5500, start2=5800)

  val intervals = new IntervalList(builder.header)
  intervals.add(new Interval(builder.header.getSequence(1).getSequenceName, 5000, 6000))

  "FilterBam" should "run ok on an empty bam" in {
    val empty = new SamBuilder().toTempFile()
    val out   = newBam
    new FilterBam(input=empty, output=out).execute()
    val reader = SamSource(out)
    reader.iterator.hasNext shouldBe false
  }

  it should "not remove any records" in {
    val out  = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeDuplicates=false, removeUnmappedReads=false, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.toSeq.map(_.name)
  }

  it should "remove all but the 'ok' records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, minMapQ=10).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filter(_.name.startsWith("ok_")).toSeq.map(_.name)
  }

  it should "remove just the unmapped records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=true, removeDuplicates=false, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filterNot(_.name.startsWith("unmapped_")).toSeq.map(_.name)
  }

  it should "remove just the secondary records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=false, minMapQ=0, removeSecondaryAlignments=true).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filterNot(_.name.startsWith("secondary_")).toSeq.map(_.name)
  }

  it should "remove just the duplicate records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=true, minMapQ=0, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filterNot(_.name.startsWith("duplicate_")).toSeq.map(_.name)
  }

  it should "remove just the low mapq records" in {
    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, removeDuplicates=false, minMapQ=20, removeSecondaryAlignments=false).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filterNot(_.name.startsWith("lowmapq_")).toSeq.map(_.name)
  }

  it should "only return matching reads within the intervals" in {
    val out = newBam
    val ilist = makeTempFile("filter_bam_test.", ".interval_list")
    intervals.write(ilist.toFile)

    new FilterBam(input=builder.toTempFile(), output=out, minMapQ=20, intervals=Some(ilist)).execute()
    val recs = readBamRecs(out)
    recs.map(_.name) should contain theSameElementsAs builder.filter(_.name.startsWith("ok_ontarget")).toSeq.map(_.name)
  }
}
