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

import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder}
import htsjdk.variant.vcf.VCFFileReader

class UpdateVcfContigNamesTest extends UnitSpec {

  private val testDir        = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf/update_vcf_contig_names/")
  private val sourceDictPath = testDir.resolve("GRCh38.p12.dict")
  private val targetDictPath = testDir.resolve("hg38.dict")
  private val sourceDict     = SequenceDictionary.extract(this.sourceDictPath)
  private val targetDict     = SequenceDictionary.extract(this.targetDictPath)

  "UpdateVcfContigNames" should "update the contig names" in {


    val builder = new VariantContextSetBuilder()
    builder.setSequenceDictionary(sourceDict)
    builder.addVariant(refIdx = 0, start = 1, variantAlleles = List("A", "C")) // chr1
    builder.addVariant(refIdx = 10, start = 10, variantAlleles = List("A", "C")) // chr2

    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = builder.toTempFile(), output = output, dict = this.targetDictPath, skipMissing = false)
    executeFgbioTool(tool)

    val reader = new VCFFileReader(output)

    // Check the sequence dictionary
    val outputDict = reader.getFileHeader.getSequenceDictionary.fromSam
    outputDict.sameAs(targetDict) shouldBe true

    // Check the variants
    val variants = reader.toIndexedSeq
    variants should have size 2
    val first = variants.head
    first.getContig shouldBe "chr1"
    first.getStart shouldBe 1
    val second = variants.last
    second.getContig shouldBe "chr2"
    second.getStart shouldBe 10

    reader.close()
  }

  it should "fail when a contig cannot be updated" in {
    val builder = new VariantContextSetBuilder()
    builder.setSequenceDictionary(sourceDict)
    builder.addVariant(refIdx = 594, start = 1, variantAlleles = List("A", "C")) // dummy

    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = builder.toTempFile(), output = output, dict = this.targetDictPath, skipMissing = false)

    an[Exception] should be thrownBy executeFgbioTool(tool)

  }

  it should "skip contigs that cannot be updated when --skip-missing is used" in {
    val builder = new VariantContextSetBuilder()
    builder.setSequenceDictionary(sourceDict)
    builder.addVariant(refIdx = 594, start = 1, variantAlleles = List("A", "C")) // dummy

    val output = makeTempFile("output", ".vcf.gz")
    val tool   = new UpdateVcfContigNames(input = builder.toTempFile(), output = output, dict = this.targetDictPath, skipMissing = true)
    executeFgbioTool(tool)

    val reader = new VCFFileReader(output)

    // Check the sequence dictionary
    val outputDict = reader.getFileHeader.getSequenceDictionary.fromSam
    outputDict.sameAs(targetDict) shouldBe true

    // Check the variants
    val variants = reader.toIndexedSeq
    variants should have size 0
  }
}
