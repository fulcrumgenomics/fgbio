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

package com.fulcrumgenomics.testing

import com.fulcrumgenomics.FgBioDef._

/**
  * Tests for VariantContextSetBuilder.
  */
class VariantContextSetBuilderTest extends  UnitSpec {
  "VariantContextSetBuilder" should "add genotypes for two samples" in {
    val sampleNames = Seq("Sample1", "Sample2")
    val builder     = new VariantContextSetBuilder(sampleNames=sampleNames)
      .addVariant(refIdx=0, start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A"), sampleName=Some(sampleNames.head))
      .addVariant(refIdx=0, start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("C"), sampleName=Some(sampleNames.last))
    builder.iterator.next().getGenotypes.size shouldBe 2
    builder.iterator.next().getGenotypesOrderedByName.map(_.getSampleName).toSeq should contain theSameElementsInOrderAs sampleNames
    builder.size shouldBe 1
  }

  it should "add a no call genotype" in {
    val builder     = new VariantContextSetBuilder()
      .addVariant(refIdx=0, start=1, variantAlleles=List("A", "C"))
    val vcf = builder.toTempFile()
    val variants = readVcfRecs(vcf).toIndexedSeq
    variants.length shouldBe 1
    variants.head.getGenotype(0).isNoCall shouldBe true
  }

  it should "add genotype attributes" in {
    val builder     = new VariantContextSetBuilder()
      .addVariant(
        refIdx=0,
        start=1,
        variantAlleles=List("A", "C"),
        genotypeAttributes = Map(
          "GQ" -> 2,
          "DP" -> 4,
          "AD" -> Seq(1, 3),
          "PL" -> Seq(1, 2, 3),
          "ZZ" -> 3
        )
      )
    val vcf = builder.toTempFile()
    val variants = readVcfRecs(vcf).toIndexedSeq
    variants.length shouldBe 1
    val genotype = variants.head.getGenotype(0)
    genotype.isNoCall shouldBe true
    genotype.getGQ shouldBe 2
    genotype.getDP shouldBe 4
    genotype.getAD should contain theSameElementsInOrderAs Seq(1, 3)
    genotype.getPL should contain theSameElementsInOrderAs Seq(1, 2, 3)
    genotype.getAnyAttribute("ZZ") shouldBe "3"
  }
}
