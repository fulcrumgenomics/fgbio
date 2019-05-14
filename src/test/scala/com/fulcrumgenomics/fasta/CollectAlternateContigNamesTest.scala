/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.fasta.AssemblyReportColumn._
import com.fulcrumgenomics.fasta.SequenceRole._

class CollectAlternateContigNamesTest extends UnitSpec {

  private val dir        = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/fasta")
  private val reportHg19 = dir.resolve("GRCh37.p13.assembly_report.txt")
  private val reportHg38 = dir.resolve("GRCh38.p12.assembly_report.txt")

  "CollectAlternateContigNames" should "read get the UCSC-style-names for alternates in GRCh37.p13" in {
    val output = makeTempFile("test.", ".txt")
    val tool = new CollectAlternateContigNames(
      input         = reportHg19,
      output        = output,
      primary       = RefSeqAccession,
      alternates    = Seq(UcscName),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold)
    )
    executeFgbioTool(tool)

    val lines = Io.readLines(output).toList
    lines should have length 84
    lines.head shouldBe "NC_000001.10\tchr1"
    lines.last shouldBe "NT_113948.1\tchr19_gl000208_random"
  }

  it should "read get the UCSC-style-names for alternates in GRCh38.p12" in {
    val output = makeTempFile("test.", ".txt")
    val tool = new CollectAlternateContigNames(
      input         = reportHg38,
      output        = output,
      primary       = RefSeqAccession,
      alternates    = Seq(UcscName),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold)
    )
    executeFgbioTool(tool)

    val lines = Io.readLines(output).toList
    lines should have length 193
    lines.head shouldBe "NC_000001.11\tchr1"
    lines.last shouldBe "NT_187405.1\tchrUn_KI270312v1"
  }

  it should "read get the UCSC-style-names and GenBank accessions for alternates in GRCh37.p13 and output a single alternate per line" in {
    val output = makeTempFile("test.", ".txt")
    val tool = new CollectAlternateContigNames(
      input         = reportHg19,
      output        = output,
      primary       = RefSeqAccession,
      alternates    = Seq(UcscName, GenBankAccession),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold),
      singleLine = false
    )
    executeFgbioTool(tool)

    val lines = Io.readLines(output).toList
    lines should have length 168
    lines.head shouldBe "NC_000001.10\tchr1"
    lines.last shouldBe "NT_113948.1\tGL000208.1"
  }

  it should "read get the UCSC-style-names and GenBank accessions for alternates in GRCh37.p13 and output all alternates on a single line" in {
    val output = makeTempFile("test.", ".txt")
    val tool = new CollectAlternateContigNames(
      input         = reportHg19,
      output        = output,
      primary       = RefSeqAccession,
      alternates    = Seq(UcscName, GenBankAccession),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold),
      singleLine = true
    )
    executeFgbioTool(tool)

    val lines = Io.readLines(output).toList
    lines should have length 84
    lines.head shouldBe "NC_000001.10\tchr1\tCM000663.1"
    lines.last shouldBe "NT_113948.1\tchr19_gl000208_random\tGL000208.1"
  }
}
