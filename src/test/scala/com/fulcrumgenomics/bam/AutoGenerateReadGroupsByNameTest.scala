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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.Iso8601Date

class AutoGenerateReadGroupsByNameTest extends UnitSpec {

  "RunInfo" should "throw an exception if a read name is malformed" in {
    // Fail
    an[Exception] should be thrownBy RunInfoFromRead("")
    an[Exception] should be thrownBy RunInfoFromRead("field")
    an[Exception] should be thrownBy RunInfoFromRead("1:2:3:4:5:6")
    an[Exception] should be thrownBy RunInfoFromRead("1:2:3:4:5:6:7:8")
    an[Exception] should be thrownBy RunInfoFromRead("1:2:3:a:5:6:7")

    // Ok
    RunInfoFromRead("1:2:3:4:5:6:7")
    RunInfoFromRead("EAS139:136:FC706VJ:2:5:1000:12850")
  }

  private def getHeader(bam: PathToBam): SAMFileHeader = {
    val in = SamSource(bam)
    yieldAndThen(in.header) { in.safelyClose() }
  }

  private def runAddReadGroups(name: String*): PathToBam = {
    val builder = new SamBuilder()
    name.foreach { n => builder.addFrag(name=n, unmapped=true) }
    val in = builder.toTempFile()
    val out = makeTempFile("AddReadGroupsByNameTest", ".bam")
    new AutoGenerateReadGroupsByName(input=in, output=out, sample="sample", library="library").execute()
    out
  }

  "AutoGenerateReadGroupsByName" should "add a single read group for reads from one lane of a flowcell" in {
    val out = runAddReadGroups(
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:3:3:4",
      "instrument:run-number:flowcell-id:1:2:4:4",
      "instrument:run-number:flowcell-id:1:2:3:5"
    )
    val readGroups = getHeader(out).getReadGroups.toSeq
    readGroups should have size 1
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id.1")
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
  }

  it should "add a two read groups for reads from two lanes in one flowcell" in {
    val out = runAddReadGroups(
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:2:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId).toSeq should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit).toSeq should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.2")
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
  }

  it should "add two read groups for reads from the first lane for two flowcells" in {
    val out = runAddReadGroups(
      "instrument:run-number:flowcell-id-1:1:2:3:4",
      "instrument:run-number:flowcell-id-2:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId).toSeq should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit).toSeq should contain theSameElementsAs Seq("flowcell-id-1.1", "flowcell-id-2.1")
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
  }

  it should "add two read groups for reads from two different instruments" in {
    val out = runAddReadGroups(
      "instrument-1:run-number:flowcell-id:1:2:3:4",
      "instrument-2:run-number:flowcell-id:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId).toSeq should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit).toSeq should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.1")
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
  }

  it should "add two read groups for reads from two different run numbers" in {
    val out = runAddReadGroups(
      "instrument:run-number-1:flowcell-id:1:2:3:4",
      "instrument:run-number-2:flowcell-id:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId).toSeq should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit).toSeq should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.1")
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
  }

  it should "accept all optional parameters" in {
    val names = Seq(
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:3:3:4",
      "instrument:run-number:flowcell-id:1:2:4:4",
      "instrument:run-number:flowcell-id:1:2:3:5"
    )
    val builder = new SamBuilder()
    names.foreach { name => builder.addFrag(name=name, unmapped=true) }
    val in = builder.toTempFile()
    val out = makeTempFile("AddReadGroupsByNameTest", ".bam")
    new AutoGenerateReadGroupsByName(
      input=in,
      output=out,
      sample="sample",
      library="library",
      sequencingCenter = Some("sequencingCenter"),
      predictedInsertSize = Some(255),
      programGroup = Some("programGroup"),
      platformModel = Some("platformModel"),
      description = Some("description"),
      runDate = Some(new Iso8601Date("2001-01-01")),
      comments = List("comment-1", "comment-2")
    ).execute()

    val header     = getHeader(out)
    val readGroups = header.getReadGroups
    readGroups should have size 1
    readGroups.foreach { rg => rg.getId shouldBe "1" }
    readGroups.foreach { rg => rg.getPlatformUnit shouldBe "flowcell-id.1" }
    readGroups.foreach { rg => rg.getSample shouldBe "sample" }
    readGroups.foreach { rg => rg.getLibrary shouldBe "library" }
    readGroups.foreach { rg => rg.getSequencingCenter shouldBe "sequencingCenter" }
    readGroups.foreach { rg => rg.getPredictedMedianInsertSize shouldBe 255 }
    readGroups.foreach { rg => rg.getProgramGroup shouldBe "programGroup" }
    readGroups.foreach { rg => rg.getPlatformModel shouldBe "platformModel" }
    readGroups.foreach { rg => rg.getDescription shouldBe "description" }
    readGroups.foreach { rg => rg.getRunDate shouldBe new Iso8601Date("2001-01-01") }
    readGroups.foreach { rg => rg.getPlatform shouldBe "ILLUMINA" }

    header.getComments should contain theSameElementsAs Seq("comment-1", "comment-2").map(co => "@CO\t" + co)
  }
}
