/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.cmdline

import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.SAMFileHeader

/** Some basic test for the CLP classes. */
class FgBioToolTest extends UnitSpec {
  private val toolInfo = FgBioToolInfo(
    name                    = "foo",
    args                    = "--these are --some args".split(" ").toList,
    commandLineWithDefaults = "--these are --some args --with defaults",
    description             = "a description",
    version                 = "0.1.2"
  )

  private def validateProgramGroup(header: SAMFileHeader, id: String) = {
    val pg = header.getProgramRecord(id)
    pg should not be null
    pg.getProgramName shouldBe toolInfo.name
    pg.getCommandLine shouldBe toolInfo.commandLineWithDefaults
    pg.getProgramVersion shouldBe toolInfo.version
    pg.getPreviousProgramGroupId shouldBe null
  }

  "FgBioToolInfo.addProgramGroupTo" should "add a program groups to an empty header" in {
    val header = new SAMFileHeader()

    // no ID given
    toolInfo.addProgramGroupTo(header=header).getId shouldBe "1"
    header.getProgramRecords.size() shouldBe 1
    validateProgramGroup(header=header, id="1")

    // ID given
    toolInfo.addProgramGroupTo(header=header, id=Some("id")).getId shouldBe "id"
    header.getProgramRecords.size() shouldBe 2
    validateProgramGroup(header=header, id="id")
  }

  it should "add a program group to a header with existing program groups with ID collisions" in {
    val header = new SAMFileHeader()

    // no ID given, so ID should be 1
    toolInfo.addProgramGroupTo(header=header).getId shouldBe "1"
    header.getProgramRecords.size() shouldBe 1
    validateProgramGroup(header=header, id="1")

    // no ID given, collision, so ID should be 2
    toolInfo.addProgramGroupTo(header=header).getId shouldBe "2"
    header.getProgramRecords.size() shouldBe 2
    validateProgramGroup(header=header, id="2")
  }
}