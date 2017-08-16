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
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}

class FilterBamTest extends UnitSpec {
  def newBam: PathToBam = makeTempFile("filter_bam_test.", ".bam")

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

  it should "remove reads that are not part of a mapped pair" in {
    val builder = new SamBuilder(readLength=50)
    builder.addPair(name="q1", start1=100, start2=200)
    builder.addFrag(name="q2", start=150)
    builder.addFrag(name="q3", unmapped=true)
    builder.addPair(name="q4", start1=200, unmapped2=true)
    builder.addPair(name="q5", unmapped1=true, start2=300)
    builder.addPair(name="q6", unmapped1=true, unmapped2=true)

    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=true, minMapQ=0, removeSingleEndMappings=true).execute()
    val recs = readBamRecs(out)
    recs should have size 2
    recs.map(_.name).toSet shouldBe Set("q1")
  }

  it should "remove reads below the specified insert size" in {
    val builder = new SamBuilder(readLength=10)
    builder.addPair("q1", start1=100, start2=589) // isize=499
    builder.addPair("q2", start1=100, start2=590) // isize=500
    builder.addPair("q3", start1=100, start2=591) // isize=501
    builder.addPair("q4", start1=100, start2=120) // isize=30
    builder.addPair("q5", start1=200, start2=100, strand1=Minus, strand2=Plus) // isize=110
    builder.addFrag("q6", start=500)
    builder.addPair("q7", start1=100, start2=500).foreach(r => r.refIndex = if (r.firstOfPair) 1 else 2) // isize=0/undefined

    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, minMapQ=0, minInsertSize=Some(500)).execute()
    val recs = readBamRecs(out)
    recs should have size 4
    recs.map(_.name).toSet shouldBe Set("q2", "q3")
  }

  it should "remove reads above the specified insert size" in {
    val builder = new SamBuilder(readLength=10)
    builder.addPair("q1", start1=100, start2=589) // isize=499
    builder.addPair("q2", start1=100, start2=590) // isize=500
    builder.addPair("q3", start1=100, start2=591) // isize=501
    builder.addPair("q4", start1=100, start2=120) // isize=30
    builder.addPair("q5", start1=200, start2=100, strand1=Minus, strand2=Plus) // isize=110
    builder.addFrag("q6", start=500)

    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, minMapQ=0, maxInsertSize=Some(499)).execute()
    val recs = readBamRecs(out)
    recs should have size 7
    recs.map(_.name).toSet shouldBe Set("q1", "q4", "q5", "q6")
  }

  it should "remove reads with fewer than 50 mapped bases" in {
    val builder = new SamBuilder(readLength=100)
    builder.addPair("q01", start1=100, start2=400, cigar1="100M",      cigar2="100M")      // keep both
    builder.addPair("q02", start1=100, start2=400, cigar1="100=",      cigar2="49=2X49=")  // keep both
    builder.addPair("q03", start1=100, start2=400, cigar1="40M40I20M", cigar2="40M10D60M") // keep both
    builder.addPair("q04", start1=100, start2=400, cigar1="80S20M",    cigar2="100M")      // keep R2 only
    builder.addPair("q05", start1=100, start2=400, cigar1="40S20M40S", cigar2="50S25M25S") // keep neither
    builder.addPair("q06", start1=100, unmapped2=true, cigar1="100M")                      // keep R1 ony
    builder.addPair("q07", unmapped1=true, unmapped2=true)                                 // keep neither
    builder.addFrag("q08", unmapped=true)                                                  // discard fragment
    builder.addFrag("q09", start=100, cigar="80S20M")                                      // discard fragment
    builder.addFrag("q10", start=100, cigar="51M49S")                                      // keep fragment

    val out = newBam
    new FilterBam(input=builder.toTempFile(), output=out, removeUnmappedReads=false, minMapQ=0, minMappedBases=Some(50)).execute()
    val recs = readBamRecs(out)
    recs.map(_.id) should contain theSameElementsAs Seq("q01/1", "q01/2", "q02/1", "q02/2", "q03/1", "q03/2", "q04/2", "q06/1", "q10/1")
  }
}
