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
 *
 */

package com.fulcrumgenomics.bam

import java.nio.file.Files

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.{SAMReadGroupRecord, SAMRecordSetBuilder}
import com.fulcrumgenomics.bam.SplitType._
import com.fulcrumgenomics.bam.api.SamOrder

class SplitBamTest extends UnitSpec {

  private def newOutput = Files.createTempDirectory("split_bam").resolve("prefix")

  private val defaultRgId: String = new SamBuilder().header.getReadGroups.map(_.getId).toSeq.head
  private val otherRgId: String = s"$defaultRgId.2"

  // Make a set of records to run multiple tests on
  private val builder = {
    val b = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some(defaultRgId))
    b.header.addReadGroup(new SAMReadGroupRecord(otherRgId))
    b
  }
  builder.addPair(name="rg_default_1", contig=0, start1=100, start2=300, attrs=Map("RG" -> defaultRgId))
  builder.addPair(name="rg_default_2", contig=0, start1=100, start2=300, attrs=Map("RG" -> defaultRgId))
  builder.addPair(name="rg_other",     contig=0, start1=100, start2=300, attrs=Map("RG" -> otherRgId))
  builder.addPair(name="rg_unknown_1", contig=0, start1=100, start2=300, attrs=Map("RG" -> null))
  builder.addPair(name="rg_unknown_2", contig=0, start1=100, start2=300, attrs=Map("RG" -> null))
  builder.addPair(name="rg_unknown_3", contig=0, start1=100, start2=300, attrs=Map("RG" -> null))
  builder.header.getReadGroups.foreach { rg => rg.setLibrary("some-library"); rg.setSample("Sample") }

  "SplitBamByReadGroup" should "split reads into multiple outputs by read group" in {
    val input = builder.toTempFile()
    val output = newOutput
    val splitter = new SplitBam(input=input, output=output, splitBy=ReadGroup)
    splitter.execute()

    builder.header.getReadGroups.foreach { rg =>
      val path = splitter.toOutput(rg.getId)
      val bams = readBamRecs(path)
      val n    = if (rg.getId == defaultRgId) 4 else 2
      bams.length shouldBe n
      bams.foreach { bam => bam.readGroup.getId shouldBe rg.getId }
    }

    // unknown
    {
      val path = splitter.toOutput(splitter.unknown)
      val bams = readBamRecs(path)
      bams.length shouldBe 6
      bams.foreach { bam => bam.readGroup shouldBe null }
    }
  }

  it should "split reads into multiple outputs by library" in {
    val input = builder.toTempFile()
    val output = newOutput
    val splitter = new SplitBam(input=input, output=output, splitBy=Library)
    splitter.execute()

    {
      val library = "some-library"
      val path = splitter.toOutput(library)
      val bams = readBamRecs(path)
      bams.length shouldBe 6
      bams.map(_.readGroup.getLibrary).foreach { library => library shouldBe "some-library" }
    }

    // unknown
    {
      val path = splitter.toOutput(splitter.unknown)
      val bams = readBamRecs(path)
      bams.length shouldBe 6
      bams.foreach { bam => bam.readGroup shouldBe null }
    }
  }
}
