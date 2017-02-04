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

package com.fulcrumgenomics.basecalling

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.{DelimitedDataParser, Io, ReadStructure}
import com.fulcrumgenomics.FgBioDef.FilePath
import htsjdk.samtools.util.Iso8601Date

class ExtractIlluminaRunInfoTest extends UnitSpec {

  /** Creates the RunInfo.xml text and writes it to a file.
    * Ignores the ImageDimensions and ImageChannels elements under Run, as well as TileSet in FlowcellLayout. */
  private def runInfo(date: String, readStructure: ReadStructure): FilePath = {
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

  "RunInfo" should "parse a RunInfo.xml with a six character date" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T")))
    info.runDate shouldBe new Iso8601Date("2017-02-04")
  }

  it should "parse a RunInfo.xml with a eight character date" in {
    val info = RunInfo(runInfo(date="20170204", readStructure=ReadStructure("8B150T")))
    info.runDate shouldBe new Iso8601Date("2017-02-04")
  }

  it should "handle a non-indexed run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("150T")))
    info.readStructure.toString shouldBe "150T"
  }

  it should "handle a single index run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T")))
    info.readStructure.toString shouldBe "8B150T"
  }

  it should "handle a dual indexed run" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B150T150T8B")))
    info.readStructure.toString shouldBe "8B150T150T8B"
  }

  it should "handle a complicated read structure" in {
    val info = RunInfo(runInfo(date="170204", readStructure=ReadStructure("8B10T8B10T10T8B")))
    info.runBarcode             shouldBe "BCDEFGHIJ_NS123456"
    info.flowcellBarcode        shouldBe "NS123456"
    info.runDate                shouldBe  new Iso8601Date("2017-02-04")
    info.readStructure.toString shouldBe "8B10T8B10T10T8B"
    info.numLanes               shouldBe 4
  }

  "ExtractRunInfo" should "run end-to-end" in {
    val in = runInfo(date="20170204", readStructure=ReadStructure("8B150T"))
    val out = makeTempFile("out", ".csv")
    new ExtractIlluminaRunInfo(input=in, output=out).execute()
    val parser = DelimitedDataParser(out, '\t')
    parser.headers should contain theSameElementsInOrderAs ExtractIlluminaRunInfo.HeaderColumns
    parser.hasNext shouldBe true
    val row = parser.next()
    row[String]("run_barcode")      shouldBe "BCDEFGHIJ_NS123456"
    row[String]("flowcell_barcode") shouldBe "NS123456"
    row[String]("run_date")         shouldBe new Iso8601Date("2017-02-04").toString
    row[String]("read_structure")   shouldBe "8B150T"
    row[String]("num_lanes")        shouldBe "4"
  }
}
