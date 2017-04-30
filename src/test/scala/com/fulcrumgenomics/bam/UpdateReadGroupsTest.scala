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

package com.fulcrumgenomics.bam

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._

import scala.collection.JavaConverters._

/**
  * Tests for UpdateReadGroups.
  */
class UpdateReadGroupsTest extends UnitSpec {
  def newBam = makeTempFile("update_read_group_test.", ".bam")

  private val builderOne   = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some("IDA"))
  private val builderTwo   = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some("IDB"))

  Stream(builderOne, builderTwo).foreach { builder =>
    builder.addPair("ok_1" + builder.readGroupId.get, contig=0, start1=100, start2=300)
    builder.addPair("ok_2" + builder.readGroupId.get, contig=0, start1=200, start2=400)
    builder.header.getReadGroups.foreach { rg => rg.setDescription("Description") }
  }

  private val headerAtoB              = makeHeader(("IDA", "IDB"))
  private val headerAandBtoCandD      = makeHeader(("IDA", "IDC"), ("IDB", "IDD"))
  private val headerBtoC              = makeHeader(("IDB", "IDC"))
  private val headerAtoBWithNewSample = makeHeaderWithNewSample(("IDA", "IDB", "NoName"))

  /** toFr are tuples, where each triple is old read group ID ("FR"), the new read group ID, and new sample name ("SM"). */
  def makeHeaderWithNewSample(fromAndTo: (String, String, String)*): SAMFileHeader = {
    val header = makeHeader(fromAndTo.map(x => (x._1, x._2)):_*)
    fromAndTo.foreach { case ((fr, to, sn)) => header.getReadGroup(to).setSample(sn) }
    header
  }

  /** toFr are tuples, where each tuple is old read group ID ("FR") and the new read group ID. */
  def makeHeader(fromAndTo: (String, String)*): SAMFileHeader = {
    val header = new SAMFileHeader()
    fromAndTo.foreach { case (fr, to) =>
      val rg = new SAMReadGroupRecord(to)
      rg.setSample("Sample")
      rg.setAttribute("FR", fr)
      header.addReadGroup(rg)
    }
    header
  }

  /** Writes the given SAM file header to a BAM file. */
  def headerToTempFile(header: SAMFileHeader): PathToBam = {
    val path = Files.createTempFile("SamRecordSet.", ".bam")
    path.toFile.deleteOnExit()
    SamWriter(path, header).close()
    path
  }

  /** Creates a temp SAM file after merging the header and records from multiple builders. */
  def mergeToTempFile(builder: SamBuilder*): PathToBam = {
    val path = Files.createTempFile("SamRecordSet.", ".bam")
    path.toFile.deleteOnExit()
    val merger = new SamFileHeaderMerger(SortOrder.coordinate, builder.map(_.header).asJava, false)
    val header = merger.getMergedHeader
    merger.hasReadGroupCollisions shouldBe false
    val writer = SamWriter(path, header, sort=Some(SamOrder.Coordinate))
    builder.foreach { b => writer ++= b }
    writer.close()
    path
  }

  "UpdateReadGroups" should "change the read group ID from IDA to IDB with an extra read group in the input file" in {
    val out = newBam
    new UpdateReadGroups(input=builderOne.toTempFile(), output=out, readGroupsFile=headerToTempFile(headerAtoB), ignoreMissingReadGroups=true).execute()
    val reader = SamSource(out)
    reader.foreach { rec => rec.readGroup.getId shouldBe headerAtoB.getReadGroups.iterator().next().getId }
    reader.close()
  }

  it should "change the read group ID from IDA to IDB" in {
    // remove the read group with ID "1" from the builder :/
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some("IDA"))
    builder.addPair("ok_1" + builder.readGroupId.get, contig=0, start1=100, start2=300)
    builder.addPair("ok_2" + builder.readGroupId.get, contig=0, start1=200, start2=400)
    builder.header.setReadGroups(builder.header.getReadGroups.filter { rg => rg.getId == "IDA" }.toJavaList)
    val out = newBam
    new UpdateReadGroups(input=builder.toTempFile(), output=out, readGroupsFile=headerToTempFile(headerAtoB), ignoreMissingReadGroups=false).execute()
    val reader = SamSource(out)
    reader.foreach { rec => rec.readGroup.getId shouldBe headerAtoB.getReadGroups.iterator().next().getId }
    reader.close()
  }

  it should "change two read groups (IDA, IDB) to (IDC, IDD)" in {
    val in  = mergeToTempFile(builderOne, builderTwo)
    val out = newBam
    new UpdateReadGroups(input=in, output=out, readGroupsFile=headerToTempFile(headerAandBtoCandD), ignoreMissingReadGroups=true).execute()
    val reader = SamSource(out)
    val builderOneReadNames = builderOne.iterator.map(_.name).toSet
    reader.foreach { rec =>
      if (builderOneReadNames.contains(rec.name)) rec.readGroup.getId shouldBe "IDC"
      else rec.readGroup.getId shouldBe "IDD"
    }
    reader.close()
  }

  it should "fail if missing read groups when ignoreMissingReadGroups is false" in {
    val out = newBam
    val in  = newBam
    val header = builderOne.header.clone()
    val extraRg = new SAMReadGroupRecord("X")
    extraRg.setKeySequence("XX")
    header.addReadGroup(extraRg)
    val writer = SamWriter(in, header, sort=Some(SamOrder.Coordinate))
    writer ++= builderOne.iterator
    writer.close()

    an[Exception] should be thrownBy new UpdateReadGroups(input=in, output=out, readGroupsFile=headerToTempFile(headerAtoB), ignoreMissingReadGroups=false).execute()
  }

  it should "allow missing read groups when ignoreMissingReadGroups is true" in {
    val out = newBam
    // headerAtoB is missing read group "ID:1".
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some("IDA"))
    val rg1 = new SAMReadGroupRecord("1")
    rg1.setSample("SomeOtherSample")
    builder.header.addReadGroup(rg1)

    builder.addPair("ok_1", contig=0, start1=100, start2=300).foreach { rec => rec(SAMTag.RG.name) = "IDA" }
    builder.addPair("ok_2", contig=0, start1=200, start2=400).foreach { rec => rec(SAMTag.RG.name) = "1" }
    new UpdateReadGroups(input=builder.toTempFile(), output=out, readGroupsFile=headerToTempFile(headerAtoB), ignoreMissingReadGroups=true).execute()
    val reader = SamSource(out)
    reader.iterator.foreach { rec =>
      if (rec.name == "ok_1") rec.readGroup.getId shouldBe "IDB"
      else rec.readGroup.getId shouldBe "1"
    }
    reader.close()
  }

  it should "not keep any read group attributes if keepReadGroupAttributes is false" in {
    val out = newBam
    new UpdateReadGroups(
      input                   = builderOne.toTempFile(),
      output                  = out,
      readGroupsFile          = headerToTempFile(headerAtoBWithNewSample),
      ignoreMissingReadGroups = true,
      keepReadGroupAttributes = false
    ).execute()
    val reader = SamSource(out)
    reader.foreach { rec =>
      rec.readGroup.getId shouldBe headerAtoBWithNewSample.getReadGroups.iterator().next().getId
      rec.readGroup.getSample shouldBe "NoName"
      rec.readGroup.getDescription shouldBe null
    }
    reader.close()
  }

  it should "keep read group attributes if keepReadGroupAttributes is true" in {
    val out = newBam
    new UpdateReadGroups(
      input                   = builderOne.toTempFile(),
      output                  = out,
      readGroupsFile          = headerToTempFile(headerAtoBWithNewSample),
      ignoreMissingReadGroups = true,
      keepReadGroupAttributes = true
    ).execute()
    val reader = SamSource(out)
    reader.foreach { rec =>
      rec.readGroup.getId shouldBe headerAtoBWithNewSample.getReadGroups.iterator().next().getId
      rec.readGroup.getSample shouldBe "NoName" // overwritten
      rec.readGroup.getDescription shouldBe "Description" // kept!
    }
    reader.close()
  }

  it should "fail if the new SAM header file has a read group without a \"FR\" attribute" in {
    val out = newBam
    an[Exception] should be thrownBy new UpdateReadGroups(
      input                   = builderOne.toTempFile(),
      output                  = out,
      readGroupsFile          = headerToTempFile(builderOne.header),
      ignoreMissingReadGroups = true
    ).execute()
  }

  it should "fail if the new SAM header file has a \"FR\" attribute that is not found in the input SAM file" in {
    val out = newBam
    an[Exception] should be thrownBy new UpdateReadGroups(
      input                   = builderOne.toTempFile(),
      output                  = out,
      readGroupsFile          = headerToTempFile(headerBtoC),
      ignoreMissingReadGroups = true
    ).execute()
  }

  it should "fail if a record in the input file has no read group" in {
    val out = newBam
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readGroupId=Some("IDA"))
    builder.addPair("ok_1" + builder.readGroupId.get, contig=0, start1=100, start2=300).foreach { rec => rec("RG") = null }
    an[Exception] should be thrownBy new UpdateReadGroups(
      input                   = builder.toTempFile(),
      output                  = out,
      readGroupsFile          = headerToTempFile(headerAtoB),
      ignoreMissingReadGroups = true
    ).execute()
  }
}
