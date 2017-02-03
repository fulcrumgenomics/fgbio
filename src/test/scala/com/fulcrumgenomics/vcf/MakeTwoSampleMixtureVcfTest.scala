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
 */

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder}
import htsjdk.variant.variantcontext.Allele
import com.fulcrumgenomics.FgBioDef._
import htsjdk.variant.vcf.VCFFileReader

class MakeTwoSampleMixtureVcfTest extends UnitSpec {
  private val builder = new VariantContextSetBuilder(sampleNames = List("s1", "s2"))
  private val (_A, _C, _G, _T, _N) = ("A", "C", "G", "T", Allele.NO_CALL_STRING)

  def addVariant(pos: Int, refAllele: String, s1Allele1: String, s1Allele2: String, s2Allele1: String, s2Allele2: String) = {
    val alleles = List(refAllele) ++ Set(s1Allele1, s1Allele2, s2Allele1, s2Allele2).filterNot(_ == refAllele).filterNot(_ == _N)
    builder.addVariant(start=pos, variantAlleles=alleles, sampleName=Some("s1"), genotypeAlleles=List(s1Allele1, s1Allele2))
    builder.addVariant(start=pos, variantAlleles=alleles, sampleName=Some("s2"), genotypeAlleles=List(s2Allele1, s2Allele2))
  }

  addVariant(10, _A, _A, _A, _A, _A)  // Monomorphic, should not come out
  addVariant(20, _A, _A, _C, _A, _C)
  addVariant(30, _A, _C, _C, _A, _A)
  addVariant(40, _A, _A, _A, _C, _C)
  addVariant(50, _A, _C, _C, _C, _A)
  addVariant(60, _A, _C, _C, _C, _C)

  addVariant(100, _A, _C, _T, _C, _T) // should be no-called for tri-allelic (including ref)
  addVariant(200, _A, _C, _T, _A, _C) // should be no-called for tri-allelic
  addVariant(300, _A, _N, _N, _A, _C) // should be no-called for missing_gt

  "MakeTwoSampleMixtureVcf" should "make a tumor only VCF" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.6, tumorOnly=true, noCallIsHomRef=false).execute()

    val reader = new VCFFileReader(out.toFile, false)
    val variants = reader.toIndexedSeq

    variants.size shouldBe 8
    reader.getFileHeader.getSampleNamesInOrder.size() shouldBe 1

    variants(0).getStart shouldBe 20
    variants(0).isFiltered shouldBe false
    variants(0).getReference.getBaseString shouldBe _A
    variants(0).getGenotype(0).getGenotypeString shouldBe "A/C"
    variants(0).getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 0.5

    variants(1).getStart shouldBe 30
    variants(1).isFiltered shouldBe false
    variants(1).getReference.getBaseString shouldBe _A
    variants(1).getGenotype(0).getGenotypeString shouldBe "A/C"
    variants(1).getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 0.6

    variants(2).getStart shouldBe 40
    variants(2).isFiltered shouldBe false
    variants(2).getReference.getBaseString shouldBe _A
    variants(2).getGenotype(0).getGenotypeString shouldBe "A/C"
    variants(2).getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 0.4

    variants(3).getStart shouldBe 50
    variants(3).isFiltered shouldBe false
    variants(3).getReference.getBaseString shouldBe _A
    variants(3).getGenotype(0).getGenotypeString shouldBe "A/C"
    variants(3).getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 0.8

    variants(4).getStart shouldBe 60
    variants(4).isFiltered shouldBe false
    variants(4).getReference.getBaseString shouldBe _A
    variants(4).getGenotype(0).getGenotypeString shouldBe "C/C"
    variants(4).getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 1.0

    variants(5).getStart shouldBe 100
    variants(5).isFiltered shouldBe true
    variants(5).getFilters.contains(MakeTwoSampleMixtureVcf.MultiAllelicFilter) shouldBe true

    variants(6).getStart shouldBe 200
    variants(6).isFiltered shouldBe true
    variants(6).getFilters.contains(MakeTwoSampleMixtureVcf.MultiAllelicFilter) shouldBe true

    variants(7).getStart shouldBe 300
    variants(7).isFiltered shouldBe true
    variants(7).getFilters.contains(MakeTwoSampleMixtureVcf.UnknownGtFilter) shouldBe true
  }

  it should "treat no-calls as hom-ref" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.5, tumorOnly=true, noCallIsHomRef=true).execute()
    val reader = new VCFFileReader(out.toFile, false)
    val variants = reader.toIndexedSeq
    val variant = variants.find(_.getStart == 300).getOrElse(fail("Could not find expected variant"))

    variant.isFiltered shouldBe false
    variant.getGenotype(0).getExtendedAttribute("AF").asInstanceOf[String].toDouble shouldBe 0.25
  }

  it should "create a tumor/normal VCF" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.5, tumorOnly=false, noCallIsHomRef=true).execute()
    val reader = new VCFFileReader(out.toFile, false)
    reader.getFileHeader.getGenotypeSamples.toSeq shouldBe Seq("tumor", "normal")
    val variants = reader.toIndexedSeq

    // The variants at the following positions have the alt allele in the normal and should get flagged as such in t/n mode
    Seq(20, 40, 50, 60, 300).map(pos => variants.find(_.getStart == pos)).foreach {
      case None => fail("Missing variant.")
      case Some(v) =>
        v.isFiltered shouldBe true
        v.getFilters.contains(MakeTwoSampleMixtureVcf.GermlineFilter) shouldBe true
    }
  }
}
