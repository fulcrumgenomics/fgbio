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

package com.fulcrumgenomics.vcf.filtration

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec, VariantContextSetBuilder}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.variant.vcf.VCFFileReader

class FilterSomaticVcfTest extends UnitSpec {
  // A pair of (paths to) VCFs for use in testing below
  private lazy val Seq(tumorOnlyVcf, tumorNormalVcf) = {
    Seq(Seq("tumor"), Seq("tumor", "normal")).map { samples =>
      val builder = new VariantContextSetBuilder(samples)
      val includeNormal = samples.contains("normal")

      // C>A SNP at 100
      builder.addVariant(start=100, sampleName=Some("tumor"), variantAlleles=List("C", "A"), genotypeAlleles=List("C", "A"))
      if (includeNormal) builder.addVariant(start=100, sampleName=Some("normal"), variantAlleles=List("C", "A"), genotypeAlleles=List("C", "C"))

      // A>G SNP at 200
      builder.addVariant(start=200, sampleName=Some("tumor"), variantAlleles=List("A", "G"), genotypeAlleles=List("A", "G"))
      if (includeNormal) builder.addVariant(start=200, sampleName=Some("normal"), variantAlleles=List("A", "G"), genotypeAlleles=List("A", "A"))

      // Indel at 300
      builder.addVariant(start=300, sampleName=Some("tumor"), variantAlleles=List("AAA", "A"), genotypeAlleles=List("AAA", "A"))
      if (includeNormal) builder.addVariant(start=300, sampleName=Some("normal"), variantAlleles=List("AAA", "A"), genotypeAlleles=List("AAA", "AAA"))

      // A>T SNP at 400
      builder.addVariant(start=400, sampleName=Some("tumor"), variantAlleles=List("A", "T"), genotypeAlleles=List("A", "T"))
      if (includeNormal) builder.addVariant(start=400, sampleName=Some("normal"), variantAlleles=List("A", "T"), genotypeAlleles=List("A", "A"))

      builder.toTempFile()
    }
  }

  // BAM file that has reads over each of the variant sites
  private lazy val bam = {
    val rlen    = 40
    val builder = new SamRecordSetBuilder(readLength=rlen, baseQuality=40, sortOrder=SortOrder.coordinate)

    // C>A SNP at 100 with signal for end repair artifact
    for (pos <- 61 to 100; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen).foreach(_.setReadString("C"*rlen))
      builder.addPair(start1=pos-rlen, start2=pos     ).foreach(_.setReadString("C"*rlen))
    }
    builder.addPair(start1=20, start2=62).foreach(_.setReadString("A"*rlen))
    builder.addPair(start1=21, start2=61).foreach(_.setReadString("A"*rlen))
    builder.addPair(start1=16, start2=63).foreach(_.setReadString("A"*rlen))

    // A>G SNP at 200 with low MAF but even reads
    for (pos <- 161 to 200; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen).foreach(_.setReadString("A"*rlen))
      builder.addPair(start1=pos-rlen, start2=pos     ).foreach(_.setReadString("A"*rlen))
      if (i == 1 && pos % 4 == 0) {
        builder.addPair(start1=pos,      start2=pos+rlen).foreach(_.setReadString("G"*rlen))
        builder.addPair(start1=pos-rlen, start2=pos     ).foreach(_.setReadString("G"*rlen))
      }
    }

    // No reads at the indel position because we won't look at those anyway

    // A>T SNP at 400 with signal for end repair artifact
    for (pos <- 361 to 400; i <- 1 to 3) {
      builder.addPair(start1=pos,      start2=pos+rlen).foreach(_.setReadString("A"*rlen))
      builder.addPair(start1=pos-rlen, start2=pos     ).foreach(_.setReadString("A"*rlen))
    }
    builder.addPair(start1=398, start2=430).foreach(_.setReadString("T"*rlen))
    builder.addPair(start1=398, start2=440).foreach(_.setReadString("T"*rlen))
    builder.addPair(start1=399, start2=437).foreach(_.setReadString("T"*rlen))
    builder.addPair(start1=400, start2=449).foreach(_.setReadString("T"*rlen))
    builder.addPair(start1=400, start2=441).foreach(_.setReadString("T"*rlen))

    builder.toTempFile()
  }

  private val EndRepairInfoKey   = new EndRepairArtifactLikelihoodFilter().InfoKey
  private val EndRepairFilterKey = new EndRepairArtifactLikelihoodFilter().FilterKey

  "FilterSomaticVcf" should "work on an empty VCF" in {
    val emptyVcf    = new VariantContextSetBuilder(Seq("tumor")).toTempFile()
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=emptyVcf, output=filteredVcf, bam=bam).execute()
    val reader = new VCFFileReader(filteredVcf.toFile, false)
    reader.getFileHeader.getInfoHeaderLines.exists(_.getID == EndRepairInfoKey) shouldBe true
    reader.getFileHeader.getFilterLines.exists(_.getID == EndRepairFilterKey) shouldBe true
    reader.safelyClose()
  }

  it should "work on a single sample VCF" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 4
    variants(0).hasAttribute(EndRepairInfoKey) shouldBe true
    variants(1).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(2).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(3).hasAttribute(EndRepairInfoKey) shouldBe true
    variants.exists(_.getFilters.contains(EndRepairFilterKey)) shouldBe false // no threshold == no filtering
  }

  it should "fail on a single-sample VCF if an invalid sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam, sample=Some("WhoDis")).execute() }
  }

  it should "fail on a multi-sample VCF if no sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam).execute() }
  }

  it should "fail on a multi-sample VCF if an invalid sample name is provided" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    an[Exception] shouldBe thrownBy { new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam, sample=Some("WhoDis")).execute() }
  }

  it should "work on a multi-sample VCF if a sample name is given" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorNormalVcf, output=filteredVcf, bam=bam, sample=Some("tumor")).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 4
    variants(0).hasAttribute(EndRepairInfoKey) shouldBe true
    variants(1).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(2).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(3).hasAttribute(EndRepairInfoKey) shouldBe true
    variants.exists(_.getFilters.contains(EndRepairFilterKey)) shouldBe false // no threshold == no filtering
  }

  it should "apply filters if the end repair p-value threshold is supplied" in {
    val filteredVcf = makeTempFile("filtered.", ".vcf")
    new FilterSomaticVcf(input=tumorOnlyVcf, output=filteredVcf, bam=bam, sample=Some("tumor"), endRepairDistance=4, endRepairPValue=Some(0.001)).execute()
    val variants = readVcfRecs(filteredVcf)
    variants should have size 4
    variants(0).hasAttribute(EndRepairInfoKey) shouldBe true
    variants(1).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(2).hasAttribute(EndRepairInfoKey) shouldBe false
    variants(3).hasAttribute(EndRepairInfoKey) shouldBe true

    variants(0).getFilters.contains(EndRepairFilterKey) shouldBe true
    variants(1).getFilters.contains(EndRepairFilterKey) shouldBe false
    variants(2).getFilters.contains(EndRepairFilterKey) shouldBe false
    variants(3).getFilters.contains(EndRepairFilterKey) shouldBe true
  }
}
