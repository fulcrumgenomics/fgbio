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
import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder}
import com.fulcrumgenomics.vcf.MakeMixtureVcf.Sample
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.mutable
import scala.collection.JavaConversions.iterableAsScalaIterable

class MakeMixtureVcfTest extends UnitSpec {
  private val dummyVcf = makeTempFile("dummy.", ".vcf")

  "MakeMixtureVcf.determineSamples" should "auto-populate samples from the header" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix")
    val ss      = mixer.determineSamples(builder.header)
    ss shouldBe samples.map(s => Sample(s, 0.25))
  }

  it should "work when only a subset of samples are specified" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s2", "s3"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder.header)
    ss shouldBe Seq(Sample("s2", 0.5), Sample("s3", 0.5))
  }

  it should "work when specific sample proportions are specified" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s2@0.9", "s3@0.1"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder.header)
    ss shouldBe Seq(Sample("s2", 0.9), Sample("s3", 0.1))
  }

  it should "divide the remaining proportion equally between unannotated samples" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.5", "s2", "s3"), outputSampleName="mix")
    val ss      = mixer.determineSamples(builder.header)
    ss shouldBe Seq(Sample("s1", 0.5), Sample("s2", 0.25), Sample("s3", 0.25))
  }

  it should "throw an exception if a sample is specified that's not in the VCF" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1", "s2", "s7"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder.header) }
  }

  it should "throw an exception if the proportions add up to less than 1" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.25", "s2@0.25"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder.header) }
  }

  it should "throw an exception if the proportions add up to more than 1" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.75", "s2@0.75"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder.header) }
  }

  it should "throw an exception if a sample proportion is negative" in {
    val samples = Seq("s1", "s2", "s3")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.75", "s2@0.75", "s3@-0.5"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder.header) }
  }

  it should "throw an exception if the same sample is specified more than once" in {
    val samples = Seq("s1", "s2", "s3", "s4")
    val builder = new VariantContextSetBuilder(samples)
    val mixer   = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, samples=Seq("s1@0.5", "s1@0.5"), outputSampleName="mix")
    an[Exception] should be thrownBy { mixer.determineSamples(builder.header) }
  }

  "MakeMixtureVcf.updateAlleleFractionsForSample" should "extract simple fractions when no AF field is used" in {
    val samples = Seq("s1", "s2")
    val builder = new VariantContextSetBuilder(samples)
    builder.addVariant(start=10, sampleName=Some("s1"), variantAlleles=List("A", "C"), genotypeAlleles=List("A"))
    builder.addVariant(start=10, sampleName=Some("s2"), variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"))
    val ctx = builder.iterator.next()

    val mixer     = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix")
    val fractions = mutable.Map[Allele,Double](ctx.getAlleles.map(a => a -> 0d).toSeq:_*)

    mixer.updateAlleleFractionsForSample(ctx, Sample("s1", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.5
    fractions(ctx.getAllele("C")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(ctx, Sample("s2", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.75
    fractions(ctx.getAllele("C")) shouldBe 0.25
  }

  it should "extract a single-valued AF from a field" in {
    val samples = Seq("s1", "s2")
    val builder = new VariantContextSetBuilder(samples)
    MakeMixtureVcf.HeaderLines.foreach(builder.addMetaDataLine)
    builder.addVariant(start=10, sampleName=Some("s1"), variantAlleles=List("A", "C"), genotypeAlleles=List("A"))
    builder.addVariant(start=10, sampleName=Some("s2"), variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), genotypeAttributes=Map("AF" -> 0.4))
    val vcf   = builder.toTempFile()
    val in    = new VCFFileReader(vcf.toFile, false)
    val ctx   = in.iterator().next()
    val mixer = new MakeMixtureVcf(input=vcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = mutable.Map[Allele,Double](ctx.getAlleles.map(a => a -> 0d).toSeq:_*)

    mixer.updateAlleleFractionsForSample(ctx, Sample("s1", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.5
    fractions(ctx.getAllele("C")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(ctx, Sample("s2", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.8
    fractions(ctx.getAllele("C")) shouldBe 0.2
  }

  it should "extract a multi-valued AF from a field" in {
    val samples = Seq("s1", "s2")
    val builder = new VariantContextSetBuilder(samples)
    MakeMixtureVcf.HeaderLines.foreach(builder.addMetaDataLine)
    builder.addVariant(start=10, sampleName=Some("s1"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A"))
    builder.addVariant(start=10, sampleName=Some("s2"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A", "C", "T"),
      genotypeAttributes=Map("AF" -> Array[Double](0.4, 0.1)))

    val vcf   = builder.toTempFile()
    val in    = new VCFFileReader(vcf.toFile, false)
    val ctx   = in.iterator().next()
    val mixer = new MakeMixtureVcf(input=vcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = mutable.Map[Allele,Double](ctx.getAlleles.map(a => a -> 0d).toSeq:_*)

    mixer.updateAlleleFractionsForSample(ctx, Sample("s1", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.5
    fractions(ctx.getAllele("C")) shouldBe 0.0
    fractions(ctx.getAllele("T")) shouldBe 0.0

    mixer.updateAlleleFractionsForSample(ctx, Sample("s2", 0.5), fractions)
    fractions(ctx.getAllele("A")) shouldBe 0.75
    fractions(ctx.getAllele("C")) shouldBe 0.2
    fractions(ctx.getAllele("T")) shouldBe 0.05
  }

  it should "fail if a sample with non-hom-ref genotypes is missing it's AF attribute" in {
    val samples = Seq("s1", "s2")
    val builder = new VariantContextSetBuilder(samples)
    MakeMixtureVcf.HeaderLines.foreach(builder.addMetaDataLine)
    builder.addVariant(start=10, sampleName=Some("s1"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A"))
    builder.addVariant(start=10, sampleName=Some("s2"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A", "C"))

    val ctx   = builder.iterator.next()
    val mixer = new MakeMixtureVcf(input=dummyVcf, output=dummyVcf, outputSampleName="mix", alleleFractionField=Some("AF"))
    val fractions = mutable.Map[Allele,Double](ctx.getAlleles.map(a => a -> 0d).toSeq:_*)

    an[Exception] should be thrownBy { mixer.updateAlleleFractionsForSample(ctx, Sample("s2", 0.5), fractions) }
  }

  Seq(true, false).foreach { noCallsAreHomRef =>
    "MakeMixtureVcf" should s"run end to end and create a valid output VCF with no-calls-are-hom-ref=${noCallsAreHomRef}" in {
      val samples = Seq("s1", "s2", "s3") // s1=0.5, s2=0.25, s3=0.25
      val builder = new VariantContextSetBuilder(samples)
      MakeMixtureVcf.HeaderLines.foreach(builder.addMetaDataLine)

      // First variant should come back as A/C with AF = 0.1
      builder.addVariant(start=10, sampleName=Some("s1"), variantAlleles=List("A","C"), genotypeAlleles=List("A"))
      builder.addVariant(start=10, sampleName=Some("s2"), variantAlleles=List("A","C"), genotypeAlleles=List("A"))
      builder.addVariant(start=10, sampleName=Some("s3"), variantAlleles=List("A","C"), genotypeAlleles=List("A", "C"),
        genotypeAttributes=Map("AF" -> 0.4))

      // Second variant should come back as C/C with AF = 1
      builder.addVariant(start=20, sampleName=Some("s1"), variantAlleles=List("A","C"), genotypeAlleles=List("C"), genotypeAttributes=Map("AF" -> 1.0))
      builder.addVariant(start=20, sampleName=Some("s2"), variantAlleles=List("A","C"), genotypeAlleles=List("C"), genotypeAttributes=Map("AF" -> 1.0))
      builder.addVariant(start=20, sampleName=Some("s3"), variantAlleles=List("A","C"), genotypeAlleles=List("C"), genotypeAttributes=Map("AF" -> 1.0))

      // Third variant should come back as A/C/T with AFs = 0.0625, 0.0625
      builder.addVariant(start=30, sampleName=Some("s1"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A"))
      builder.addVariant(start=30, sampleName=Some("s2"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A","T"), genotypeAttributes=Map("AF" -> 0.25))
      builder.addVariant(start=30, sampleName=Some("s3"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A","C"), genotypeAttributes=Map("AF" -> 0.25))

      // Fourth variant should come back as A/C/T with AFs = 0.0625, 0.0625 when noCallsAreHomRef, and a no-call otherwise
      builder.addVariant(start=40, sampleName=Some("s1"), variantAlleles=List("A","C","T"), genotypeAlleles=List(".", "."))
      builder.addVariant(start=40, sampleName=Some("s2"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A","T"), genotypeAttributes=Map("AF" -> 0.25))
      builder.addVariant(start=40, sampleName=Some("s3"), variantAlleles=List("A","C","T"), genotypeAlleles=List("A","C"), genotypeAttributes=Map("AF" -> 0.25))

      val in  = builder.toTempFile()
      val out = makeTempFile("mixture.", ".vcf")

      val mixer = new MakeMixtureVcf(input=in, output=out, samples=Seq("s1@0.5", "s2@0.25", "s3@0.25"),
        outputSampleName="mix", noCallIsHomRef=noCallsAreHomRef, alleleFractionField=Some("AF"))
      mixer.execute()

      val Seq(ctx1, ctx2, ctx3, ctx4) = new VCFFileReader(out.toFile, false).toSeq

      ctx1.getStart shouldBe 10
      ctx1.getGenotype("mix").getGenotypeString(true) shouldBe "A/C"
      ctx1.getGenotype("mix").getExtendedAttribute("AF").toString.toDouble shouldBe 0.1

      ctx2.getStart shouldBe 20
      ctx2.getGenotype("mix").getGenotypeString(true) shouldBe "C"
      ctx2.getGenotype("mix").getExtendedAttribute("AF").toString.toDouble shouldBe 1.0

      ctx3.getStart shouldBe 30
      ctx3.getGenotype("mix").getGenotypeString(true) shouldBe "A/C/T"
      ctx3.getGenotype("mix").getExtendedAttribute("AF").toString.split(',').map(_.toDouble) shouldBe Array(0.0625, 0.0625)

      ctx4.getStart shouldBe 40
      if (noCallsAreHomRef) {
        ctx4.getGenotype("mix").getGenotypeString(true) shouldBe "A/C/T"
        ctx4.getGenotype("mix").getExtendedAttribute("AF").toString.split(',').map(_.toDouble) shouldBe Array(0.0625, 0.0625)
      }
      else {
        ctx4.getGenotype("mix").isNoCall shouldBe true
        ctx4.isFiltered shouldBe true
      }
    }
  }
}
