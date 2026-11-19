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
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.MakeMixtureVcf.Sample
import com.fulcrumgenomics.vcf.api.{Allele, Variant, VcfContigHeader, VcfCount, VcfFieldType, VcfFilterHeader, VcfFormatHeader, VcfGeneralHeader, VcfInfoHeader, VcfSource}

import scala.collection.mutable

class MakeMixtureVcfTest extends UnitSpec {
  private val dummyVcf = makeTempFile("dummy.", ".vcf")

  /** A VCF builder for the given samples with a FORMAT header line for the allele fraction (AF). */
  private def builder(samples: String*): VcfBuilder = {
    VcfBuilder(samples).withFormatHeaders(MakeMixtureVcf.HeaderLines.collect { case f: VcfFormatHeader => f }:_*)
  }

  /** Builds a map of allele fractions, all zero, for the alleles of the variant. */
  private def emptyFractions(variant: Variant): mutable.Map[Allele, Double] = {
    mutable.Map[Allele, Double](variant.alleles.toSeq.map(a => a -> 0d):_*)
  }

  "MakeMixtureVcf.determineSamples" should "auto-populate samples from the header" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix")
    val ss      = mixer.determineSamples(builder(samples:_*).header)
    ss shouldBe samples.map(s => Sample(s, 0.25))
  }

  it should "work when only a subset of samples are specified" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s2", "s3"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder("s1", "s2", "s3", "s4").header)
    ss shouldBe Seq(Sample("s2", 0.5), Sample("s3", 0.5))
  }

  it should "work when specific sample proportions are specified" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s2@0.9", "s3@0.1"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder("s1", "s2", "s3", "s4").header)
    ss shouldBe Seq(Sample("s2", 0.9), Sample("s3", 0.1))
  }

  it should "divide the remaining proportion equally between unannotated samples" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.5", "s2", "s3"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder("s1", "s2", "s3", "s4").header)
    ss shouldBe Seq(Sample("s1", 0.5), Sample("s2", 0.25), Sample("s3", 0.25))
  }

  it should "throw an exception if a sample is specified that's not in the VCF" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1", "s2", "s7"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder("s1", "s2", "s3", "s4").header) }
  }

  it should "throw an exception if the proportions add up to less than 1" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.25", "s2@0.25"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder("s1", "s2", "s3", "s4").header) }
  }

  it should "throw an exception if the proportions add up to more than 1" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.75", "s2@0.75"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder("s1", "s2", "s3", "s4").header) }
  }

  it should "throw an exception if a sample proportion is negative" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.75", "s2@0.75", "s3@-0.5"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder("s1", "s2", "s3").header) }
  }

  it should "throw an exception if the same sample is specified more than once" in {
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.5", "s1@0.5"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder("s1", "s2", "s3", "s4").header) }
  }

  "MakeMixtureVcf.updateAlleleFractionsForSample" should "extract simple fractions when no AF field is used" in {
    val variant = builder("s1", "s2")
      .add(pos=10, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A/C")))
      .iterator.next()

    val mixer     = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix")
    val fractions = emptyFractions(variant)

    mixer.updateAlleleFractionsForSample(variant, Sample("s1", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.5
    fractions(Allele("C")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(variant, Sample("s2", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.75
    fractions(Allele("C")) shouldBe 0.25
  }

  it should "extract a single-valued AF from a field" in {
    val vcf = builder("s1", "s2")
      .add(pos=10, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A/C", attrs=Map("AF" -> 0.4))))
      .toTempFile()
    val variant   = VcfSource(vcf).iterator.next()
    val mixer     = new MakeMixtureVcf(input=vcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = emptyFractions(variant)

    mixer.updateAlleleFractionsForSample(variant, Sample("s1", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.5
    fractions(Allele("C")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(variant, Sample("s2", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.8
    fractions(Allele("C")) shouldBe 0.2
  }

  it should "extract a multi-valued AF from a field" in {
    val vcf = builder("s1", "s2")
      .add(pos=10, alleles=Seq("A", "C", "T"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A/C/T", attrs=Map("AF" -> Seq(0.4, 0.1)))))
      .toTempFile()
    val variant   = VcfSource(vcf).iterator.next()
    val mixer     = new MakeMixtureVcf(input=vcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = emptyFractions(variant)

    mixer.updateAlleleFractionsForSample(variant, Sample("s1", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.5
    fractions(Allele("C")) shouldBe 0.0
    fractions(Allele("T")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(variant, Sample("s2", 0.5), fractions)
    fractions(Allele("A")) shouldBe 0.75
    fractions(Allele("C")) shouldBe 0.2
    fractions(Allele("T")) shouldBe 0.05
  }

  it should "fail if a sample with non-hom-ref genotypes is missing it's AF attribute" in {
    val variant = builder("s1", "s2")
      .add(pos=10, alleles=Seq("A", "C", "T"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A/C")))
      .iterator.next()
    val mixer     = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = emptyFractions(variant)

    an[Exception] should be thrownBy { mixer.updateAlleleFractionsForSample(variant, Sample("s2", 0.5), fractions) }
  }

  it should "treat a genotype with a spanning deletion as heterozygous when given a single AF value" in {
    val variant = builder("s1")
      .add(pos=10, alleles=Seq("A", "*"), gts=Seq(Gt(sample="s1", gt="A/*", attrs=Map("AF" -> 0.4))))
      .iterator.next()
    val mixer     = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = emptyFractions(variant)

    mixer.updateAlleleFractionsForSample(variant, Sample("s1", 1.0), fractions)
    fractions(Allele("A")) shouldBe 0.6 +- 1e-9
    fractions(Allele("*")) shouldBe 0.4 +- 1e-9
  }

  Seq(true, false).foreach { noCallsAreHomRef =>
    "MakeMixtureVcf" should s"run end to end and create a valid output VCF with no-calls-are-hom-ref=${noCallsAreHomRef}" in {
      val samples = Seq("s1", "s2", "s3") // s1=0.5, s2=0.25, s3=0.25
      val in = builder(samples:_*)
        // First variant should come back as A/C with AF = 0.1
        .add(pos=10, alleles=Seq("A","C"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A"), Gt(sample="s3", gt="A/C", attrs=Map("AF" -> 0.4))))
        // Second variant should come back as C/C with AF = 1
        .add(pos=20, alleles=Seq("A","C"), gts=Seq(Gt(sample="s1", gt="C", attrs=Map("AF" -> 1.0)), Gt(sample="s2", gt="C", attrs=Map("AF" -> 1.0)), Gt(sample="s3", gt="C", attrs=Map("AF" -> 1.0))))
        // Third variant should come back as A/C/T with AFs = 0.0625, 0.0625
        .add(pos=30, alleles=Seq("A","C","T"), gts=Seq(Gt(sample="s1", gt="A"), Gt(sample="s2", gt="A/T", attrs=Map("AF" -> 0.25)), Gt(sample="s3", gt="A/C", attrs=Map("AF" -> 0.25))))
        // Fourth variant should come back as A/C/T with AFs = 0.0625, 0.0625 when noCallsAreHomRef, and a no-call otherwise
        .add(pos=40, alleles=Seq("A","C","T"), gts=Seq(Gt(sample="s1", gt="./."), Gt(sample="s2", gt="A/T", attrs=Map("AF" -> 0.25)), Gt(sample="s3", gt="A/C", attrs=Map("AF" -> 0.25))))
        .toTempFile()
      val out = makeTempFile("mixture.", ".vcf")

      val mixer = new MakeMixtureVcf(input=in, output=out, samples=Seq("s1@0.5", "s2@0.25", "s3@0.25"),
        outputSampleName="mix", noCallIsHomRef=noCallsAreHomRef, alleleFractionField=Some("AF"))
      mixer.execute()

      val Seq(v1, v2, v3, v4) = VcfSource(out).toSeq

      v1.pos shouldBe 10
      v1.gt("mix").gtWithBases shouldBe "A/C"
      v1.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.1f)

      v2.pos shouldBe 20
      v2.gt("mix").gtWithBases shouldBe "C"
      v2.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(1.0f)

      v3.pos shouldBe 30
      v3.gt("mix").gtWithBases shouldBe "A/C/T"
      v3.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.0625f, 0.0625f)

      v4.pos shouldBe 40
      if (noCallsAreHomRef) {
        v4.gt("mix").gtWithBases shouldBe "A/C/T"
        v4.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.0625f, 0.0625f)
      }
      else {
        v4.gt("mix").isNoCall shouldBe true
        v4.filters should contain (MakeMixtureVcf.UnknownGtFilter)
      }
    }
  }

  "MakeMixtureVcf" should "keep only the alleles called in the mixed samples, the ID and the input filters" in {
    val in = VcfBuilder(Seq("s1", "s2", "s3"))
      .add(pos=10, id="rs10", alleles=Seq("A", "C", "T"), qual=50, info=Map("DP" -> 30), filters=Seq("PASS"),
        gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0"), Gt(sample="s3", gt="0/2")))
      .add(pos=20, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/0"), Gt(sample="s2", gt="0/0"), Gt(sample="s3", gt="1/1")))
      .add(pos=30, alleles=Seq("A", "C"), filters=Seq("LowQD"), gts=Seq(Gt(sample="s1", gt="1/1"), Gt(sample="s2", gt="0|1"), Gt(sample="s3", gt="0/0")))
      .add(pos=40, alleles=Seq("A", "C", "G"), gts=Seq(Gt(sample="s1", gt="1/2"), Gt(sample="s2", gt="2/2"), Gt(sample="s3", gt="0/0")))
      .add(pos=50, alleles=Seq("A", "C"), filters=Seq("LowQD"), gts=Seq(Gt(sample="s1", gt="./."), Gt(sample="s2", gt="0/1"), Gt(sample="s3", gt="0/0")))
      .toTempFile()

    Seq(true, false).foreach { noCallIsHomRef =>
      val out = makeTempFile("mixture.", ".vcf.gz")
      new MakeMixtureVcf(input=in, output=out, samples=Seq("s1", "s2"), outputSampleName="mix", noCallIsHomRef=noCallIsHomRef).execute()

      val source = VcfSource(out)
      source.header.samples should contain theSameElementsInOrderAs Seq("mix")
      val Seq(v10, v20, v30, v40, v50) = source.toIndexedSeq
      source.safelyClose()

      // Only the alleles called in s1 and s2 are kept; the QUAL and INFO are dropped
      v10.id shouldBe Some("rs10")
      v10.alleles.toSeq.map(_.value) should contain theSameElementsInOrderAs Seq("A", "C")
      v10.qual shouldBe None
      v10.attrs shouldBe empty
      v10.filters shouldBe Set("PASS")
      v10.gt("mix").gtWithBases shouldBe "A/C"
      v10.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.25f)

      // Both mixed samples are hom-ref, so the mixture is a haploid hom-ref with no alt alleles
      v20.alleles.alts shouldBe empty
      v20.gt("mix").gtWithBases shouldBe "A"
      v20.gt("mix").attrs.contains("AF") shouldBe false
      v20.filters shouldBe Set("PASS")

      v30.filters shouldBe Set("LowQD")
      v30.gt("mix").gtWithBases shouldBe "A/C"
      v30.gt("mix").phased shouldBe false
      v30.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.75f)

      v40.alleles.toSeq.map(_.value) should contain theSameElementsInOrderAs Seq("A", "C", "G")
      v40.gt("mix").gtWithBases shouldBe "C/G"
      v40.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.25f, 0.75f)

      if (noCallIsHomRef) {
        v50.filters shouldBe Set("LowQD")
        v50.gt("mix").gtWithBases shouldBe "A/C"
        v50.gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.25f)
      }
      else {
        v50.filters shouldBe Set("LowQD", MakeMixtureVcf.UnknownGtFilter)
        v50.gt("mix").isNoCall shouldBe true
      }
    }
  }

  it should "replace the input's definitions of the header lines it adds" in {
    val in = VcfBuilder(Seq("s1", "s2"))
      .withFormatHeaders(VcfFormatHeader(MakeMixtureVcf.AlleleFractionField, VcfCount.Fixed(1), VcfFieldType.Float, "A different AF."))
      .withFilterHeaders(VcfFilterHeader(MakeMixtureVcf.UnknownGtFilter, "A different filter."))
      .add(pos=10, alleles=Seq("A", "C", "G"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/2")))
      .toTempFile()
    val out = makeTempFile("mixture.", ".vcf")
    new MakeMixtureVcf(input=in, output=out, outputSampleName="mix").execute()

    val source = VcfSource(out)
    source.header.format(MakeMixtureVcf.AlleleFractionField).count shouldBe VcfCount.OnePerAltAllele
    source.header.filter(MakeMixtureVcf.UnknownGtFilter).description shouldBe MakeMixtureVcf.HeaderLines.collect { case f: VcfFilterHeader => f.description }.head
    source.iterator.next().gt("mix")[Seq[Float]]("AF") should contain theSameElementsInOrderAs Seq(0.25f, 0.25f)
    source.safelyClose()
  }

  it should "write the allele fractions with the requested precision" in {
    val in = VcfBuilder(Seq("s1", "s2", "s3"))
      .add(pos=10, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="0/0"), Gt(sample="s3", gt="0/0")))
      .toTempFile()
    Seq(2 -> "0.17", 5 -> "0.16667").foreach { case (precision, expected) =>
      val out = makeTempFile("mixture.", ".vcf")
      new MakeMixtureVcf(input=in, output=out, outputSampleName="mix", precision=precision).execute()
      val record = Io.readLines(out).filterNot(_.startsWith("#")).next().split('\t')
      record(8) shouldBe "GT:AF"
      record(9) shouldBe s"0/1:$expected"
    }
  }

  "MakeMixtureVcf.makeWriter" should "add the given INFO, FORMAT and FILTER lines, replacing any of the same ID in the input" in {
    val input = VcfBuilder(Seq("s1"))
      .withInfoHeaders(VcfInfoHeader("XX", VcfCount.Fixed(1), VcfFieldType.String, "The input's XX."))
      .header
    val info  = VcfInfoHeader("XX", VcfCount.Fixed(1), VcfFieldType.Integer, "Our XX.")
    val out   = makeTempFile("mixture.", ".vcf")
    MakeMixtureVcf.makeWriter(out, input, MakeMixtureVcf.HeaderLines :+ info, "mix").close()

    val source = VcfSource(out)
    source.header.info("XX") shouldBe info
    source.header.format(MakeMixtureVcf.AlleleFractionField).description shouldBe "Alt allele fraction."
    source.header.filter.contains(MakeMixtureVcf.UnknownGtFilter) shouldBe true
    source.header.samples should contain theSameElementsInOrderAs Seq("mix")
    source.safelyClose()
  }

  it should "fail if given a header line other than an INFO, FORMAT or FILTER line" in {
    val input = VcfBuilder(Seq("s1")).header
    Seq(VcfContigHeader(index=0, name="chrX"), VcfGeneralHeader("ALT", "DEL")).foreach { line =>
      val ex = intercept[IllegalArgumentException] { MakeMixtureVcf.makeWriter(makeTempFile("mixture.", ".vcf"), input, Seq(line), "mix") }
      ex.getMessage should include ("Only INFO, FORMAT and FILTER entries")
    }
  }
}
