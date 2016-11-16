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

import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import htsjdk.samtools.{SAMFileHeader, SamReaderFactory}

class AddReadGroupsByNameTest extends UnitSpec {

  "RecordInfo" should "throw an exception if a read name is malformed" in {
    // Fail
    an[Exception] should be thrownBy RunInfo("")
    an[Exception] should be thrownBy RunInfo("field")
    an[Exception] should be thrownBy RunInfo("1:2:3:4:5:6")
    an[Exception] should be thrownBy RunInfo("1:2:3:4:5:6:7:8")
    an[Exception] should be thrownBy RunInfo("1:2:3:a:5:6:7")

    // Ok
    RunInfo("1:2:3:4:5:6:7")
  }

  private def getHeader(bam: PathToBam): SAMFileHeader = {
    val in = SamReaderFactory.make().open(bam.toFile)
    val header = in.getFileHeader
    in.close()
    header
  }

  private def getAddReadGroupsByNameOutput(name: String*): PathToBam = {
    val builder = new SamRecordSetBuilder()
    name.foreach { n => builder.addFrag(name=n, unmapped=true) }
    val in = builder.toTempFile()
    val out = makeTempFile("AddReadGroupsByNameTest", ".bam")
    new AddReadGroupsByName(input=in, output=out, sample="sample", library="library").execute()
    out
  }

  import scala.collection.JavaConversions.asScalaBuffer

  "AddReadGroupsByName" should "add a single read group for reads from one lane of a flowcell" in {
    val out = getAddReadGroupsByNameOutput(
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:1:3:3:4",
      "instrument:run-number:flowcell-id:1:2:4:4",
      "instrument:run-number:flowcell-id:1:2:3:5"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 1
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id.1")
  }

  it should "add a two read groups for reads from two lanes in one flowcell" in {
    val out = getAddReadGroupsByNameOutput(
      "instrument:run-number:flowcell-id:1:2:3:4",
      "instrument:run-number:flowcell-id:2:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.2")
  }

  it should "add two read groups for reads from the first lane for two flowcells" in {
    val out = getAddReadGroupsByNameOutput(
      "instrument:run-number:flowcell-id-1:1:2:3:4",
      "instrument:run-number:flowcell-id-2:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id-1.1", "flowcell-id-2.1")
  }

  it should "add two read groups for reads from two different instruments" in {
    val out = getAddReadGroupsByNameOutput(
      "instrument-1:run-number:flowcell-id:1:2:3:4",
      "instrument-2:run-number:flowcell-id:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.1")
    readGroups.map(_.getPlatformModel) should contain theSameElementsAs Seq("instrument-1.run-number", "instrument-2.run-number")
  }

  it should "add two read groups for reads from two different run numbers" in {
    val out = getAddReadGroupsByNameOutput(
      "instrument:run-number-1:flowcell-id:1:2:3:4",
      "instrument:run-number-2:flowcell-id:1:2:3:4"
    )
    val readGroups = getHeader(out).getReadGroups
    readGroups should have size 2
    readGroups.map(_.getId) should contain theSameElementsAs Seq("1", "2")
    readGroups.map(_.getPlatformUnit) should contain theSameElementsAs Seq("flowcell-id.1", "flowcell-id.1")
    readGroups.map(_.getPlatformModel) should contain theSameElementsAs Seq("instrument.run-number-1", "instrument.run-number-2")
  }
}
