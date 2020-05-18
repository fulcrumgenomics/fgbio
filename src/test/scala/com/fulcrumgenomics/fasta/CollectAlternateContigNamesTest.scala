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
import com.fulcrumgenomics.fasta.{AssemblyReportColumn => Column}
import com.fulcrumgenomics.fasta.SequenceRole._
import com.fulcrumgenomics.testing.UnitSpec

class CollectAlternateContigNamesTest extends UnitSpec {

  private val dir        = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/fasta")
  private val reportHg19 = dir.resolve("GRCh37.p13.assembly_report.txt")
  private val reportHg38 = dir.resolve("GRCh38.p12.assembly_report.txt")

  "CollectAlternateContigNames" should "read the UCSC-style-names for alternates in GRCh37.p13" in {
    val output = makeTempFile("test.", ".dict")
    val tool = new CollectAlternateContigNames(
      input         = reportHg19,
      output        = output,
      primary       = Column.RefSeqAccession,
      alternates    = Seq(Column.UcscName),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold)
    )
    executeFgbioTool(tool)

    val dict = SequenceDictionary(output)
    dict should have length 84
    dict.head.name shouldBe "NC_000001.10"
    dict.head.aliases should contain theSameElementsInOrderAs Seq("chr1")
    dict.last.name shouldBe "NC_012920.1"
    dict.last.aliases should contain theSameElementsInOrderAs Seq("chrM")
  }

  it should "read the UCSC-style-names for alternates in GRCh38.p12" in {
    // Note: skips molecules without a ref-seq accession
    val output = makeTempFile("test.", ".dict")
    val tool = new CollectAlternateContigNames(
      input         = reportHg38,
      output        = output,
      primary       = Column.RefSeqAccession,
      alternates    = Seq(Column.UcscName),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold)
    )
    executeFgbioTool(tool)

    val dict = SequenceDictionary(output)
    dict should have length 193
    dict.head.name shouldBe "NC_000001.11"
    dict.head.aliases should contain theSameElementsInOrderAs Seq("chr1")
    dict.last.name shouldBe "NC_012920.1"
    dict.last.aliases should contain theSameElementsInOrderAs Seq("chrM")
  }

  it should "read the assigned-molecules for alternates in GRCh38.p12" in {
    val output = makeTempFile("test.", ".dict")
    val tool = new CollectAlternateContigNames(
      input                = reportHg38,
      output                = output,
      primary               = Column.RefSeqAccession,
      alternates            = Seq(Column.AssignedMolecule),
      sequenceRoles         = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold),
      skipMissingAlternates = false
    )
    executeFgbioTool(tool)

    val dict = SequenceDictionary(output)
    dict should have length 193
    dict.head.name shouldBe "NC_000001.11"
    dict.head.aliases should contain theSameElementsInOrderAs Seq("1")
    dict.last.name shouldBe "NC_012920.1"
    dict.last.aliases should contain theSameElementsInOrderAs Seq("MT")
    // make sure only assembled-molecules have aliases
    dict.foreach { metadata: SequenceMetadata =>
      metadata.aliases.nonEmpty shouldBe (SequenceRole(metadata) == AssembledMolecule)
    }
  }

  it should "read the assigned-molecules for alternates in GRCh38.p12 but only those with aliases" in {
    val output = makeTempFile("test.", ".dict")
    val tool = new CollectAlternateContigNames(
      input                 = reportHg38,
      output                = output,
      primary               = Column.RefSeqAccession,
      alternates            = Seq(Column.AssignedMolecule),
      sequenceRoles         = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold),
      skipMissingAlternates = true
    )
    executeFgbioTool(tool)

    val dict = SequenceDictionary(output)
    dict should have length 25
    dict.head.name shouldBe "NC_000001.11"
    dict.head.aliases should contain theSameElementsInOrderAs Seq("1")
    dict.last.name shouldBe "NC_012920.1"
    dict.last.aliases should contain theSameElementsInOrderAs Seq("MT")
  }

  it should "update an existing sequence dictionary" in {
    val firstOutput  = makeTempFile("test.", ".1.dict")
    val secondOutput = makeTempFile("test.", ".2.dict")
    val firstTool = new CollectAlternateContigNames(
      input         = reportHg38,
      output        = firstOutput,
      primary       = Column.RefSeqAccession,
      alternates    = Seq(Column.GenBankAccession),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold)
    )
    executeFgbioTool(firstTool)

    val secondTool = new CollectAlternateContigNames(
      input         = reportHg38,
      output        = secondOutput,
      primary       = Column.RefSeqAccession,
      alternates    = Seq(Column.UcscName),
      sequenceRoles = Seq(AssembledMolecule, UnlocalizedScaffold, UnplacedScaffold),
      existing      = Some(firstOutput)
    )
    executeFgbioTool(secondTool)

    val dict = SequenceDictionary(secondOutput)
    dict should have length 193
    dict.head.name shouldBe "NC_000001.11"
    dict.head.aliases should contain theSameElementsInOrderAs Seq("CM000663.2", "chr1")
    dict.last.name shouldBe "NC_012920.1"
    dict.last.aliases should contain theSameElementsInOrderAs Seq("J01415.2", "chrM")
  }
}
