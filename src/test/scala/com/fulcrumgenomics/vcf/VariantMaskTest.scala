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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, UnitSpec, VariantContextSetBuilder}
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

class VariantMaskTest extends UnitSpec {
  val ref = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr0").add("NNNNNNNNNN", times=100) // length 1000
    builder.add("chr1").add("AAAAAAAAAA", times=100) // length 1000
    builder.add("chr2").add("CCCCCCCCCC", times=100) // length 1000
    builder.add("chr3").add("GGGGGGGGGG", times=100) // length 1000
    builder.add("chr4").add("TTTTTTTTTT", times=100) // length 1000
    builder.add("chr5").add("ACGTACGTAC", times=100) // length 1000
    builder.toTempFile()
  }

  val dict = SAMSequenceDictionaryExtractor.extractDictionary(ref.toFile)

  "VariantMask" should "mask SNPs as individual bases" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","T"))
    builder.addVariant(refIdx=1, start=102, variantAlleles=List("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1, 99)  shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe false
    mask.isVariant(1, 102) shouldBe true
    mask.isVariant(1, 103) shouldBe false
  }

  it should "mask all deleted bases for deletions, plus the upstream base" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("AA","A"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1,  99) shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe true
    mask.isVariant(1, 102) shouldBe false
  }

  it should "mask just the upstream base for insertions" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","AA"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1,  99) shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe false
    mask.isVariant(1, 102) shouldBe false
  }

  it should "allow querying be sequence name as well as ref index" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant("chr0", 100) shouldBe false
    mask.isVariant("chr1", 100) shouldBe true
    mask.isVariant("chr2", 100) shouldBe false
    mask.isVariant("chr3", 100) shouldBe false
    mask.isVariant("chr4", 100) shouldBe false
  }

  it should "construct a mask ok from a VCF path" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","C"), genotypeAlleles=List("A", "C"))
    val mask = VariantMask(builder.toTempFile())
    mask.isVariant(1, 100) shouldBe true
  }

  it should "throw an exception if a VCF doesn't have a sequence dictionary in it" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","C"), genotypeAlleles=List("A", "C"))
    builder.setSequenceDictionary(new SAMSequenceDictionary())
    val vcf = builder.toTempFile()
    an[Exception] shouldBe thrownBy { VariantMask(vcf) }
  }

  it should "throw an exception if requested to traverse backwards to an earlier reference" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","T"))
    builder.addVariant(refIdx=2, start=200, variantAlleles=List("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(2, 200) shouldBe true
    an[Exception] shouldBe thrownBy { mask.isVariant(1, 100) }
  }

  it should "throw an exception if invalid reference sequences or positions are requested" in {
    val builder = new VariantContextSetBuilder().setSequenceDictionary(dict)
    builder.addVariant(refIdx=1, start=100, variantAlleles=List("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    an[Exception] shouldBe thrownBy { mask.isVariant(-1, 100) }        // invalid index (low)
    an[Exception] shouldBe thrownBy { mask.isVariant("chrNope", 100) } // invalid ref name
    an[Exception] shouldBe thrownBy { mask.isVariant(0, -1) }          // invalid position (low)
    an[Exception] shouldBe thrownBy { mask.isVariant(0, 1e6.toInt) }   // invalid position (high)
    an[Exception] shouldBe thrownBy { mask.isVariant(9, 100) }         // invalid index (high)
  }
}
