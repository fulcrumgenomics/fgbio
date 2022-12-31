/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.{PathToBam, PathToIntervals}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource}
import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Io, SortOrder}
import htsjdk.samtools.util.{Interval, IntervalList}
import org.scalatest.matchers.must.Matchers.convertToAnyMustWrapper

import scala.jdk.CollectionConverters.SeqHasAsJava

/** Unit test for [[NormalizeCoverage]]. */
class NormalizeCoverageTest extends UnitSpec {

  /* test will have a region of interest consisting of initially overlapping regions

  chr1:100-250
  chr1:200-400
  chr1:410-450

  The templates will consist of reads that
   - overlap each other (mates and supplementals)
   - map to different contigs
   - are unneeded for various reasons
   - overlap multiple regions of interest
   - target coverage = 2 (so that we can test that overlaps count as 1)
   - minMq=30
*/

  private val dict = SequenceDictionary(
    SequenceMetadata(name = "chr1", length = 10000),
    SequenceMetadata(name = "chr2", length = 50000),
    SequenceMetadata(name = "chr3", length = 10000)
  ).asSam
  val chr1Index: Int = dict.getSequenceIndex("chr1")
  val chr2Index: Int = dict.getSequenceIndex("chr2")
  val chr3Index: Int = dict.getSequenceIndex("chr3")


  val i1: Interval = new Interval("chr1", 100, 250)
  val i2: Interval = new Interval("chr1", 200, 400)
  val i3: Interval = new Interval("chr1", 410, 450)

  val il: IntervalList = new IntervalList(dict)
  il.addall(List(i1, i2, i3).asJava)

  val builder = new SamBuilder(sort = Some(SamOrder.Coordinate), readLength = 100)
  var expectedReadCount: Int = 0

  builder.header.setSequenceDictionary(dict)

  //  template: chr1:1-100, chr1:200-299 (mapq0) (needed 1)
  builder.addPair("needed_1", contig = chr1Index, start1 = 1, start2 = 200, mapq2 = 0)
  expectedReadCount += 2

  //  template: chr1:1-99,chr1:451-550 (not needed 1)
  builder.addPair("not_needed_1", contig = chr1Index, start1 = 1, start2 = 451, bases1 = "A" * 99, quals1 = "F" * 99, cigar1 = "99M")

  //  template: chr1:150-249, unmapped (needed 2)
  builder.addPair("needed_2", contig = chr1Index, start1 = 150, unmapped2 = true)
  expectedReadCount += 2

  //  template: chr1:150-249 (mapq=10), unmapped (not needed)
  builder.addPair("not_needed_2", contig = chr1Index, start1 = 150, mapq1 = 10, unmapped2 = true)

  //  template: chr1:120-219, chr1:130-229 and supplementals in chr2 (needed)
  builder.addPair("needed_3", contig = chr1Index, start1 = 130, start2 = 140)
  val supplementalReads: Seq[SamRecord] = builder.addPair("needed_3", contig = chr3Index, start1 = 1200, start2 = 1300)
  supplementalReads.foreach(r => r.supplementary = true)
  expectedReadCount += 4

  //  template : chr1:120-219, chr1:130-229 duplicate (not needed 3)
  val dupReads: Seq[SamRecord] = builder.addPair("not_needed_3", contig = chr1Index, start1 = 120, start2 = 130)
  dupReads.foreach(r => r.duplicate = true)

  //  template: chr1:120-130, chr1:120-130 (mapq=10) (two but no more are needed for the overlap region)
  for (i <- 1 to 3) {
    val shortReadsA: Seq[SamRecord] = builder.addPair(s"needed_perhaps_$i", contig = chr1Index, start1 = 120, start2 = 120, mapq2 = 10)
    shortReadsA.foreach(r => {
      r.bases = "A" * 11
      r.quals = "F" * 11
      r.cigar = "11M"
    })
  }
  expectedReadCount += 4

  //  template: chr1:120-130, chr1:120-130 secondary (not needed 4)
  val secondaryReads: Seq[SamRecord] = builder.addPair("not_needed_4", contig = chr1Index, start1 = 120, start2 = 120)
  secondaryReads.foreach(r => r.secondary = true)

  //  template: chr1:380-479, chr2:350-449 (needed 4)
  builder.addPair("needed_4", contig = chr1Index, contig2 = Some(chr2Index), start1 = 380, start2 = 350)
  expectedReadCount += 2

  //  template: chr1:380-479, chr2:350-449 non-pf (not needed 5)
  val nonPfReads: Seq[SamRecord] = builder.addPair("not_needed_5", contig = chr1Index, contig2 = Some(chr2Index), start1 = 120, start2 = 120)
  nonPfReads.foreach(r => r.pf = false)

  //  template: chr2:100-199, chr3:400-499 (not needed 6)
  builder.addPair("not_needed_6", contig = chr2Index, contig2 = Some(chr3Index), start1 = 100, start2 = 400)

  //  fragment 1: chr1:10-109 (needed)
  builder.addFrag("needed_5", contig = chr1Index, start = 10)
  expectedReadCount += 1

  //  fragment 2: chr1:451-550 (not needed
  builder.addFrag("not_needed_7", contig = chr1Index, start = 451)

  val intervalListFile: PathToIntervals = Io.makeTempFile("NormalizeCoverageTest", ".list")
  intervalListFile.toFile.deleteOnExit()
  il.write(intervalListFile)

  val bamFileIn: PathToBam = Io.makeTempFile("NormalizeCoverageTest", "_In.bam")
  bamFileIn.toFile.deleteOnExit()
  builder.write(bamFileIn)

  val bamFileOut: PathToBam = Io.makeTempFile("NormalizeCoverageTest", "_Out.bam")
  bamFileOut.toFile.deleteOnExit()

  val normalizeCoverage = new NormalizeCoverage(bamFileIn, bamFileOut, Some(intervalListFile), coverageTarget = 2, sortOrder = SamOrder.Coordinate)

  normalizeCoverage.execute()

  "normalizeCoverage" should "correctly include reads in a simple test" in {
    val reads = SamSource(bamFileOut)

    var optionalCount: Int = 0
    var readCount: Int = 0
    reads.foreach(r => {
      readCount += 1
      r.name must startWith("needed")
      r.name mustNot startWith("not_needed")
      if (r.name.startsWith("needed_perhaps")) {
        optionalCount += 1
      }
    })
    optionalCount mustEqual 4 // each optional template has 2 overlapping reads, and we need two templates, so 4 reads
    expectedReadCount must be > 10 // just to be clear that this isn't empty...
    readCount mustEqual expectedReadCount
    reads.close()
  }
}


