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

package com.fulcrumgenomics.illumina

import com.fulcrumgenomics.FgBioDef.FilePath
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.{Io, ReadStructure}
import htsjdk.samtools.util.Iso8601Date

object RunInfoTest extends UnitSpec {
  /** Creates the RunInfo.xml text and writes it to a file.
    * Ignores the ImageDimensions and ImageChannels elements under Run, as well as TileSet in FlowcellLayout. */
  def runInfo(date: String, readStructure: ReadStructure): FilePath = {
    // *** WARNING ***
    // This code is used by other tests (ex. ExtractIlluminaRunInfoTest, so please do not modify.
    val pre = s"""<?xml version="1.0"?>
                 |<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="4">
                 |  <Run Id="${date}_NS123456_0011_ABCDEFGHIJ" Number="11">
                 |    <Flowcell>NS123456</Flowcell>
                 |    <Instrument>BCDEFGHIJ</Instrument>
                 |    <Date>${date}</Date>
                 |    <Reads>
                 |""".stripMargin

    val reads = readStructure.zipWithIndex.map { case (segment, idx) =>
      val isIndexedRead = if (segment.symbol == 'B') "Y" else "N"
      s"""      <Read Number="${idx+1}" NumCycles="${segment.length}" IsIndexedRead="${isIndexedRead}" />"""
    }.mkString("\n")

    val post =
      s"""
         |    </Reads>
         |    <FlowcellLayout LaneCount="4" SurfaceCount="2" SwathCount="1" TileCount="12" SectionPerLane="3" LanePerSection="2">
         |    </FlowcellLayout>
         |  </Run>
         |</RunInfo>""".stripMargin

    val text = Seq(pre, reads, post).mkString("").split("\n")
    val out  = makeTempFile(prefix="RunInfo", suffix=".xml")
    Io.writeLines(out, text)
    out
  }
}

class RunInfoTest extends UnitSpec {
  import RunInfoTest._

  "RunInfo" should "parse a RunInfo.xml with a six character date" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T")))
    info.run_date shouldBe new Iso8601Date("2017-02-04")
  }

  it should "parse a RunInfo.xml with a eight character date" in {
    val info = RunInfo(runInfo(date="20170204", readStructure=ReadStructure("8B150T")))
    info.run_date shouldBe new Iso8601Date("2017-02-04")
  }

  it should "handle a non-indexed run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("150T")))
    info.read_structure.toString shouldBe "150T"
  }

  it should "handle a single index run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T")))
    info.read_structure.toString shouldBe "8B150T"
  }

  it should "handle a dual indexed run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T150T8B")))
    info.read_structure.toString shouldBe "8B150T150T8B"
  }

  it should "handle a complicated read structure" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B10T8B10T10T8B")))
    info.run_barcode             shouldBe "BCDEFGHIJ_NS123456"
    info.flowcell_barcode        shouldBe "NS123456"
    info.instrument_name         shouldBe "BCDEFGHIJ"
    info.run_date                shouldBe  new Iso8601Date("2017-02-04")
    info.read_structure.toString shouldBe "8B10T8B10T10T8B"
    info.num_lanes               shouldBe 4
  }
}
