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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.{VcfContigHeader, VcfSource}

class UpdateVcfContigNamesTest extends UnitSpec {

  private val testDir        = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf/update_vcf_contig_names/")
  private val sourceDictPath = testDir.resolve("GRCh38.p12.dict")
  private val targetDictPath = testDir.resolve("hg38.dict")
  private val sourceDict     = SequenceDictionary.extract(this.sourceDictPath)
  private val targetDict     = SequenceDictionary.extract(this.targetDictPath)

  /** A VCF builder whose header has the contigs of the source sequence dictionary. */
  private def builder(samples: String*): VcfBuilder = {
    val contigs = VcfContigHeader.fromSequenceDictionary(sourceDict)
    VcfBuilder(VcfBuilder.DefaultHeader.copy(contigs=contigs, samples=samples.toIndexedSeq))
  }

  "UpdateVcfContigNames" should "update the contig names" in {
    val vcf = builder("s1")
      .add(chrom=sourceDict(0).name, pos=1, id="rs1", alleles=Seq("A", "C"), qual=50, info=Map("DP" -> 20),
        filters=Seq("LowQD"), gts=Seq(Gt(sample="s1", gt="0|1", ad=Seq(10, 10)))) // chr1
      .add(chrom=sourceDict(10).name, pos=10, alleles=Seq("A", "C")) // chr2
      .toTempFile()

    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = vcf, output = output, dict = this.targetDictPath, skipMissing = false)
    executeFgbioTool(tool)

    val reader = VcfSource(output)

    // Check the sequence dictionary, including that the contig lines are those of the target dictionary
    reader.header.dict.sameAs(targetDict) shouldBe true
    reader.header.contigs should contain theSameElementsInOrderAs VcfContigHeader.fromSequenceDictionary(targetDict)
    reader.header.contigs.head.md5 shouldBe targetDict(0).md5

    // Check the variants
    val variants = reader.toIndexedSeq
    variants should have size 2
    val first = variants.head
    first.chrom shouldBe "chr1"
    first.pos shouldBe 1
    first.id shouldBe Some("rs1")
    first.qual shouldBe Some(50)
    first.filters shouldBe Set("LowQD")
    first[Int]("DP") shouldBe 20
    first.gt("s1").gtVcfStyle shouldBe "0|1"
    first.gt("s1")[IndexedSeq[Int]]("AD") should contain theSameElementsInOrderAs Seq(10, 10)
    val second = variants.last
    second.chrom shouldBe "chr2"
    second.pos shouldBe 10

    reader.close()
  }

  it should "fail when a contig cannot be updated" in {
    val vcf    = builder().add(chrom=sourceDict(594).name, pos=1, alleles=Seq("A", "C")).toTempFile() // dummy
    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = vcf, output = output, dict = this.targetDictPath, skipMissing = false)

    an[Exception] should be thrownBy executeFgbioTool(tool)
  }

  it should "skip contigs that cannot be updated when --skip-missing is used" in {
    val vcf    = builder().add(chrom=sourceDict(594).name, pos=1, alleles=Seq("A", "C")).toTempFile() // dummy
    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = vcf, output = output, dict = this.targetDictPath, skipMissing = true)
    executeFgbioTool(tool)

    val reader = VcfSource(output)

    // Check the sequence dictionary
    reader.header.dict.sameAs(targetDict) shouldBe true

    // Check the variants
    reader.toIndexedSeq should have size 0
    reader.close()
  }

  it should "copy records through unchanged other than the contig name, and not need an index" in {
    val record = Seq("AC=.,3;AF=0.5,.;MQ=60.0", "GT:AD:FT:XF", "0|1:10,.:lowq:1.50,.")
    val input  = makeTempFile("input", ".vcf")
    Io.writeLines(input, Seq(
      "##fileformat=VCFv4.2",
      s"##contig=<ID=${sourceDict(0).name},length=${sourceDict(0).length}>",
      "##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">",
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">",
      "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"mq\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">",
      "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">",
      "##FORMAT=<ID=XF,Number=2,Type=Float,Description=\"xf\">",
      Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1").mkString("\t"),
      (Seq(sourceDict(0).name, "100", ".", "A", "C,G", "50", "PASS") ++ record).mkString("\t")
    ))

    val output = makeTempFile("output", ".vcf")
    executeFgbioTool(new UpdateVcfContigNames(input=input, output=output, dict=this.targetDictPath))

    Io.readLines(output).filterNot(_.startsWith("#")).toSeq shouldBe Seq((Seq("chr1", "100", ".", "A", "C,G", "50", "PASS") ++ record).mkString("\t"))
  }
}
