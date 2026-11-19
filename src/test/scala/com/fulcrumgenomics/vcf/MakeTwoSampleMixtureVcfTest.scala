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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.{Variant, VcfCount, VcfFieldType, VcfInfoHeader, VcfSource}
import htsjdk.samtools.util.{Interval, IntervalList}

import scala.annotation.unused

class MakeTwoSampleMixtureVcfTest extends UnitSpec {
  private val builder = VcfBuilder(Seq("s1", "s2"))
  @unused private val (_A, _C,  _G, _T, _N) = ("A", "C", "G", "T", ".")

  def addVariant(pos: Int, refAllele: String, s1Allele1: String, s1Allele2: String, s2Allele1: String, s2Allele2: String) = {
    val alleles = List(refAllele) ++ Set(s1Allele1, s1Allele2, s2Allele1, s2Allele2).filterNot(_ == refAllele).filterNot(_ == _N)
    builder.add(pos=pos, alleles=alleles, gts=Seq(Gt(sample="s1", gt=s"$s1Allele1/$s1Allele2"), Gt(sample="s2", gt=s"$s2Allele1/$s2Allele2")))
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

  /** Gets the single allele fraction of the tumor genotype. */
  private def tumorAf(variant: Variant): Double = variant.gt("tumor")[Seq[Float]]("AF").head.toString.toDouble

  "MakeTwoSampleMixtureVcf" should "make a tumor only VCF" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.6, tumorOnly=true, noCallIsHomRef=false).execute()

    val reader = VcfSource(out)
    val variants = reader.toIndexedSeq

    variants.size shouldBe 8
    reader.header.samples.size shouldBe 1

    variants(0).pos shouldBe 20
    variants(0).isFiltered shouldBe false
    variants(0).alleles.ref.value shouldBe _A
    variants(0).gt("tumor").gtWithBases shouldBe "A/C"
    tumorAf(variants(0)) shouldBe 0.5

    variants(1).pos shouldBe 30
    variants(1).isFiltered shouldBe false
    variants(1).alleles.ref.value shouldBe _A
    variants(1).gt("tumor").gtWithBases shouldBe "A/C"
    tumorAf(variants(1)) shouldBe 0.6

    variants(2).pos shouldBe 40
    variants(2).isFiltered shouldBe false
    variants(2).alleles.ref.value shouldBe _A
    variants(2).gt("tumor").gtWithBases shouldBe "A/C"
    tumorAf(variants(2)) shouldBe 0.4

    variants(3).pos shouldBe 50
    variants(3).isFiltered shouldBe false
    variants(3).alleles.ref.value shouldBe _A
    variants(3).gt("tumor").gtWithBases shouldBe "A/C"
    tumorAf(variants(3)) shouldBe 0.8

    variants(4).pos shouldBe 60
    variants(4).isFiltered shouldBe false
    variants(4).alleles.ref.value shouldBe _A
    variants(4).gt("tumor").gtWithBases shouldBe "C/C"
    tumorAf(variants(4)) shouldBe 1.0

    variants(5).pos shouldBe 100
    variants(5).isFiltered shouldBe true
    variants(5).filters.contains(MakeTwoSampleMixtureVcf.MultiAllelicFilter) shouldBe true

    variants(6).pos shouldBe 200
    variants(6).isFiltered shouldBe true
    variants(6).filters.contains(MakeTwoSampleMixtureVcf.MultiAllelicFilter) shouldBe true

    variants(7).pos shouldBe 300
    variants(7).isFiltered shouldBe true
    variants(7).filters.contains(MakeTwoSampleMixtureVcf.UnknownGtFilter) shouldBe true
  }

  it should "treat no-calls as hom-ref" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.5, tumorOnly=true, noCallIsHomRef=true).execute()
    val variants = VcfSource(out).toIndexedSeq
    val variant = variants.find(_.pos == 300).getOrElse(fail("Could not find expected variant"))

    variant.isFiltered shouldBe false
    tumorAf(variant) shouldBe 0.25
  }

  it should "create a tumor/normal VCF" in {
    val in = this.builder.toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.5, tumorOnly=false, noCallIsHomRef=true).execute()
    val reader = VcfSource(out)
    reader.header.samples shouldBe Seq("tumor", "normal")
    val variants = reader.toIndexedSeq

    // The variants at the following positions have the alt allele in the normal and should get flagged as such in t/n mode
    Seq(20, 40, 50, 60, 300).map(pos => variants.find(_.pos == pos)).foreach {
      case None => fail("Missing variant.")
      case Some(v) =>
        v.isFiltered shouldBe true
        v.filters.contains(MakeTwoSampleMixtureVcf.GermlineFilter) shouldBe true
    }
  }

  it should "only mix variants that start in or enclose the given intervals" in {
    val in = VcfBuilder(Seq("s1", "s2"))
      .withInfoHeaders(VcfInfoHeader(id="END", count=VcfCount.Fixed(1), kind=VcfFieldType.Integer, description="End position."))
      .add(chrom="chr1", pos=10,  alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .add(chrom="chr1", pos=95,  alleles=Seq("ACGTAC", "A"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0"))) // ends in 100-200
      .add(chrom="chr1", pos=100, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .add(chrom="chr1", pos=150, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .add(chrom="chr1", pos=250, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .add(chrom="chr1", pos=290, alleles=Seq("A", "<DEL>"), info=Map("END" -> 320), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0"))) // encloses 300-310
      .add(chrom="chr2", pos=100, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .toTempFile()

    val intervals = new IntervalList(new SamBuilder().header)
    intervals.add(new Interval("chr1", 100, 200))
    intervals.add(new Interval("chr1", 300, 310))
    intervals.add(new Interval("chr2", 1, 1000))
    val intervalsPath = makeTempFile("intervals.", ".interval_list")
    intervals.write(intervalsPath)

    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", intervals=Some(intervalsPath)).execute()

    val source = VcfSource(out)
    source.map(v => (v.chrom, v.pos)).toSeq should contain theSameElementsInOrderAs Seq(("chr1", 100), ("chr1", 150), ("chr1", 290), ("chr2", 100))
    source.safelyClose()
  }

  it should "fail with a clear message if the input VCF is not indexed and intervals are given" in {
    val in        = VcfBuilder(Seq("s1", "s2")).add(pos=10, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0"))).toTempFile()
    val unindexed = makeTempFile("mixture.", ".vcf.gz")
    java.nio.file.Files.copy(in, unindexed, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    val intervals = new IntervalList(new SamBuilder().header)
    intervals.add(new Interval("chr1", 1, 1000))
    val intervalsPath = makeTempFile("intervals.", ".interval_list")
    intervals.write(intervalsPath)

    val tool = new MakeTwoSampleMixtureVcf(input=unindexed, output=makeTempFile("mixture.", ".vcf.gz"), tumor="s1", normal="s2", intervals=Some(intervalsPath))
    val ex   = intercept[FailureException] { tool.execute() }
    ex.message.value should include ("--input VCF must be indexed")
    ex.message.value should include (unindexed.toString)
  }

  it should "set the tumor and normal genotypes, allele fraction, ID and filters" in {
    val in = VcfBuilder(Seq("s1", "s2"))
      .add(pos=10, id="rs10", alleles=Seq("A", "C"), qual=40, filters=Seq("LowQD"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0")))
      .add(pos=20, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="1/1"), Gt(sample="s2", gt="0|1")))
      .add(pos=30, alleles=Seq("A", "C", "G"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/2")))
      .add(pos=40, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="./."), Gt(sample="s2", gt="0/1")))
      .toTempFile()

    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.3, noCallIsHomRef=false).execute()

    val source = VcfSource(out)
    source.header.samples should contain theSameElementsInOrderAs Seq("tumor", "normal")
    val Seq(v10, v20, v30, v40) = source.toIndexedSeq
    source.safelyClose()

    v10.id shouldBe Some("rs10")
    v10.qual shouldBe None
    v10.filters shouldBe Set("PASS")
    v10.gt("tumor").gtWithBases shouldBe "A/C"
    v10.gt("tumor")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.15f)
    v10.gt("normal").gtWithBases shouldBe "A/A"

    v20.filters shouldBe Set(MakeTwoSampleMixtureVcf.GermlineFilter)
    v20.gt("tumor").gtWithBases shouldBe "A/C"
    v20.gt("tumor")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.65f)
    v20.gt("normal").gtWithBases shouldBe "A/C"
    v20.gt("normal").phased shouldBe false

    v30.filters shouldBe Set(MakeTwoSampleMixtureVcf.GermlineFilter, MakeTwoSampleMixtureVcf.MultiAllelicFilter, MakeTwoSampleMixtureVcf.UnknownGtFilter)
    v30.gt("tumor").isNoCall shouldBe true
    v30.gt("normal").isNoCall shouldBe true

    v40.filters shouldBe Set(MakeTwoSampleMixtureVcf.GermlineFilter, MakeTwoSampleMixtureVcf.UnknownGtFilter)
    v40.gt("tumor").isNoCall shouldBe true
    v40.gt("normal").gtWithBases shouldBe "A/C"
  }

  it should "mix the alt allele the two samples carry, even when it is not the variant's first alt allele" in {
    val in = VcfBuilder(Seq("s1", "s2"))
      .add(pos=10, alleles=Seq("A", "C", "G"), gts=Seq(Gt(sample="s1", gt="0/2"), Gt(sample="s2", gt="0/0")))
      .add(pos=20, alleles=Seq("A", "C", "G"), gts=Seq(Gt(sample="s1", gt="2/2"), Gt(sample="s2", gt="2/2")))
      .toTempFile()
    val out = makeTempFile("mixture.", ".vcf.gz")
    new MakeTwoSampleMixtureVcf(input=in, output=out, tumor="s1", normal="s2", tumorFraction=0.5).execute()

    val Seq(v10, v20) = VcfSource(out).toIndexedSeq
    v10.filters shouldBe Set("PASS")
    v10.gt("tumor").gtWithBases shouldBe "A/G"
    tumorAf(v10) shouldBe 0.25
    v10.gt("normal").gtWithBases shouldBe "A/A"

    v20.filters shouldBe Set(MakeTwoSampleMixtureVcf.GermlineFilter)
    v20.gt("tumor").gtWithBases shouldBe "G/G"
    tumorAf(v20) shouldBe 1.0
  }
}
