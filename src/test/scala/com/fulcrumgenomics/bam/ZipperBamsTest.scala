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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.bam.api.{SamOrder, SamSource}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMFileHeader.GroupOrder
import htsjdk.samtools.{SAMProgramRecord, SAMSequenceDictionary}
import com.fulcrumgenomics.FgBioDef._

class ZipperBamsTest extends UnitSpec {
  private val dict = new SamBuilder().dict

  // Read Group ID to be used in unmapped BAMs
  private val rgId = "abc"

  // Program group to be inserted into unmapped BAM
  private val unmappedPg   = {
    val tmp = new SAMProgramRecord("fqToBam")
    tmp.setCommandLine("FastqToBam do-some-stuff")
    tmp.setProgramName("FastqToBam")
    tmp
  }

  // Program group to be inserted into mapped BAM
  private val mappedPg   = {
    val tmp = new SAMProgramRecord("aligner")
    tmp.setCommandLine("aligner -t 8 -x 4 in.fastq ref.fa")
    tmp.setProgramName("SomeAligner")
    tmp
  }

  // Comments to be inserted into unmapped BAMs
  private val comments = Seq("This is a SAM Header comment.", "This is too.")

  /** Generate a new SamBuilder for unmapped BAMs */
  private def uBuilder(grouped: Boolean = true): SamBuilder = {
    val builder = new SamBuilder(readGroupId=Some(rgId))
    builder.header.setSequenceDictionary(new SAMSequenceDictionary()) // wipe out the sequence dictionary
    builder.header.addProgramRecord(unmappedPg)
    comments.foreach(builder.header.addComment)
    if (grouped) builder.header.setGroupOrder(GroupOrder.query)
    builder
  }

  /** Generate a builder that will mimic what an aligner might produce. */
  private def mBuilder(grouped: Boolean = false): SamBuilder = {
    val builder = new SamBuilder(readGroupId=None)
    builder.header.setReadGroups(Iterator().toJavaList)
    builder.header.addProgramRecord(mappedPg)
    builder
  }

  /** Run ZipperBams from two sam builders and return the output as a SamSource. */
  private def run(unmapped: SamBuilder,
                  mapped: SamBuilder,
                  remove: Iterable[String] = Nil,
                  reverse: Iterable[String] = Nil,
                  revcomp: Iterable[String] = Nil,
                  sort: Option[SamOrder] = None,
                  buffer: Int = 5000
                 ): SamSource = {
    val dir      = Io.makeTempDir("zipper-bams-test")
    val uBam     = dir.resolve("unmapped.bam")
    val mBam     = dir.resolve("mapped.bam")
    val zBam     = dir.resolve("zippered.bam")
    val dictPath = dir.resolve("ref.dict")

    Seq(uBam, mBam, zBam, dictPath, dir).foreach(_.toFile.deleteOnExit())

    dict.write(dictPath)
    unmapped.write(uBam)

    // Clean out any RG tags before writing the mapped BAM
    mapped.foreach(r => r.remove("RG"))
    mapped.write(mBam)

    val zipper = new ZipperBams(
      input         = mBam,
      unmapped      = uBam,
      ref           = dictPath,
      output        = zBam,
      tagsToRemove  = remove.toIndexedSeq,
      tagsToReverse = reverse.toIndexedSeq,
      tagsToRevcomp = revcomp.toIndexedSeq,
      sort          = sort,
      buffer        = buffer
    )

    executeFgbioTool(zipper)
    SamSource(zBam)
  }

  "ZipperBams" should "run on small BAM and copy over information from the unmapped BAM" in {
    val unmapped = uBuilder()
    val mapped   = mBuilder()

    unmapped.addPair(name="q1", unmapped1=true, unmapped2=true, attrs=Map("RX" -> "ACGT", "xy" -> 1234))
    unmapped.addPair(name="q2", unmapped1=true, unmapped2=true, attrs=Map("RX" -> "GGTA", "xy" -> 4567))

    mapped.addPair(name="q1", start1=100, start2=200, attrs=Map("PG" -> mappedPg.getId, "AS" -> 77))
    mapped.addPair(name="q2", start1=500, start2=700, attrs=Map("PG" -> mappedPg.getId, "AS" -> 16))

    val zippered = run(unmapped, mapped)

    // Check the header first
    zippered.comments should contain theSameElementsInOrderAs this.comments
    zippered.readGroups should contain theSameElementsAs  unmapped.header.getReadGroups.toSeq
    zippered.programGroups should contain theSameElementsAs Seq(unmappedPg, mappedPg)

    // Then check the records
    zippered.foreach { z =>
      // Find the matching input records
      val u = unmapped.find(u => u.name == z.name && u.firstOfPair == z.firstOfPair).value
      val m = mapped.find(m => m.name == z.name && m.flags == z.flags).value

      // All the attributes from the unmapped BAM should be there
      u.attributes.foreach { case (tag, value) => z(tag) == value shouldBe true }

      // And all the attributes from the mapped BAM
      m.attributes.foreach { case (tag, value) => z(tag) == value shouldBe true }

      // And we shouldn't have broken any of the alignment info
      z.basesString shouldBe m.basesString
      z.qualsString shouldBe m.qualsString
      z.refIndex shouldBe m.refIndex
      z.start shouldBe m.start
      z.cigar shouldBe m.cigar

      // And mate cigar should be present
      z.get("MQ").isDefined shouldBe true
    }
  }
}
