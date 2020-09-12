/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.FixVcfPhaseSet.VcfPhaseSetUpdater
import com.fulcrumgenomics.vcf.FixVcfPhaseSet.VcfPhaseSetUpdater.Result._
import com.fulcrumgenomics.vcf.api.{Genotype, Variant, VcfHeader}
import org.scalatest.OptionValues

class FixVcfPhaseSetTest extends UnitSpec with OptionValues {

  private case class TestCase(variant: Variant, genotype: Genotype, header: VcfHeader)

  private def testCase(builder: VcfBuilder = VcfBuilder(samples=Seq("sample")),
                       phased: Boolean = false,
                       phaseSet: Option[Any] = None
                       ): TestCase = {
    val gtAttrs = phaseSet.map { value => "PS" -> value }.toMap
    val gt = Gt(
      sample = "sample",
      gt     = if (phased) "0|1" else "0/1",
      attrs  = gtAttrs
    )

    builder.add(
      chrom   = "chr1",
      pos     = 1,
      alleles = Seq("A", "C"),
      gts     = Seq(gt)
    )

    val variant = builder.head
    TestCase(
      variant  = variant,
      genotype = variant.genotypes.get("sample").value,
      header   = builder.header
    )
  }

  "VcfPhaseSetUpdater" should "ignore unphased variants with no phase set value when phaseGenotypesWithPhaseSet=false" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=false, phaseSet=None)
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[Int]("PS").isEmpty shouldBe true
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe NotPhasedNoPhaseSetValue(genotype)
    updater.update(variant=variant) shouldBe variant
  }

  it should "ignore unphased variants with a phase set value if phaseGenotypesWithPhaseSet=false" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=false, phaseSet=Some("1"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").value shouldBe "1"
    val updatedGenotype = genotype.copy(attrs=Map.empty) // PS should be removed
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe NotPhasedWithPhaseSetValue(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  it should "update the phase set value for unphased variants if phaseGenotypesWithPhaseSet=true" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=false, phaseSet=Some("1"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=true)
    genotype.get[String]("PS").value shouldBe "1"
    val updatedGenotype = genotype.copy(phased=true, attrs=Map("PS" -> 1)) // note the type is changed!
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe UpdatedPhaseOnly(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  it should "ignore phased variants with no phase set" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=true, phaseSet=None)
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").isEmpty shouldBe true
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe PhasedMissingPhaseSetValue(genotype)
    updater.update(variant=variant) shouldBe variant
  }

  it should "do nothing to a phased variant with the correct the phase set value" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=true, phaseSet=Some("1"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").value shouldBe "1"
    val updatedGenotype = genotype.copy(phased=true, attrs=Map("PS" -> 1)) // note the type is changed!
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe Valid(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  it should "update the phase set if it is the wrong type/kind" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=true, phaseSet=Some("ABC"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").value shouldBe "ABC"
    val updatedGenotype = genotype.copy(phased=true, attrs=Map("PS" -> 1)) // note the type is changed!
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe UpdatedPhaseSetValue(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  it should "update the phase set if it is the wrong value" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=true, phaseSet=Some("2"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=false, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").value shouldBe "2"
    val updatedGenotype = genotype.copy(phased=true, attrs=Map("PS" -> 1)) // note the type is changed!
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe UpdatedPhaseSetValue(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  it should "keep the original value when the phase set is updated" in {
    val TestCase(variant: Variant, genotype: Genotype, header: VcfHeader) = testCase(phased=true, phaseSet=Some("123"))
    val updater = new VcfPhaseSetUpdater(header=header, keepOriginal=true, phaseGenotypesWithPhaseSet=false)
    genotype.get[String]("PS").value shouldBe "123"
    val updatedGenotype = genotype.copy(phased=true, attrs=Map("PS" -> 1, "OPS" -> "123")) // note the type is changed!
    updater.updateGenotype(variant=variant, genotype=genotype) shouldBe UpdatedPhaseSetValue(updatedGenotype)
    updater.update(variant=variant) shouldBe variant.copy(genotypes=Map("sample" -> updatedGenotype))
  }

  // Compares genotypes individually to make it a bit easier to debug later
  private def compareVariant(actual: Variant, expected: Variant): Unit = withClue(f"${actual.chrom}:${actual.pos}" ) {
    actual.genotypes.size shouldBe expected.genotypes.size
    actual.genotypes.zip(expected.genotypes).foreach { case (left, right) => left shouldBe right }
    actual shouldBe expected
  }

  private def compareVariants(actual: Seq[Variant], expected: Seq[Variant]): Unit = {
    actual.length shouldBe expected.length
    actual.zip(expected).foreach { case (left, right) => compareVariant(actual=left, expected=right) }
  }

  def build(): (VcfBuilder, Seq[Variant]) = {
    // Developer notes:
    // - sample s1 is phased and wrong phase set values
    // - sample s2 is unphased
    // - sample s3 has correct phase set values, but is not phased, so will have it updated
    // - some PS values are string, some are integers, and these should be handled gracefully
    val builder = VcfBuilder(samples=Seq("s1", "s2", "s3"))

    // S3 unphased
    builder.add(
      chrom   = "chr1",
      pos     = 1,
      alleles = Seq("A", "C"),
      gts     = Seq(Gt("s1", "0|1", attrs=Map("PS" -> "ABC")), Gt("s2", "0/1"), Gt("s3", "0/1"))
    )

    builder.add(
      chrom   = "chr1",
      pos     = 2,
      alleles = Seq("A", "C"),
      gts     = Seq(Gt("s1", "1|1", attrs=Map("PS" -> "ABC")), Gt("s2", "0/1"), Gt("s3", "0/1", attrs=Map("PS" -> 2)))
    )

    // S1 is unphased
    builder.add(
      chrom   = "chr1",
      pos     = 3,
      alleles = Seq("A", "C"),
      gts     = Seq(Gt("s1", "0/1"), Gt("s2", "0/1"), Gt("s3", "0/1", attrs=Map("PS" -> 2)))
    )

    val Seq(v1: Variant, v2: Variant, v3: Variant) = builder.toSeq
    val expected = IndexedSeq(
      v1.copy(
        genotypes = Map(
          "s1" -> v1.genotypes("s1").copy(attrs=Map("PS" -> 1)),
          "s2" -> v1.genotypes("s2"),
          "s3" -> v1.genotypes("s3"),
        )
      ),
      v2.copy(
        genotypes = Map(
          "s1" -> v2.genotypes("s1").copy(attrs=Map("PS" -> 1)),
          "s2" -> v2.genotypes("s2"),
          "s3" -> v2.genotypes("s3").copy(phased=true, attrs=Map("PS" -> 2)),
        )
      ),
      v3.copy(
        genotypes = Map(
          "s1" -> v3.genotypes("s1"),
          "s2" -> v3.genotypes("s2"),
          "s3" -> v3.genotypes("s3").copy(phased=true, attrs=Map("PS" -> 2)),
        )
      )
    )

    (builder, expected)
  }

  it should "update the phase set value to the position of the first variant in the set" in {
    val (builder: VcfBuilder, expected: Seq[Variant]) = build()
    // convert all phase set values to string, as VcfPhaseSetUpdater expects all genotypes to be String values.
    val actual: Seq[Variant] = builder.toIndexedSeq.map { variant =>
      val genotypes = variant.genotypes.map { case (sample, genotype) =>
        val maybePhaseSet = genotype.get[Any]("PS").map {value => ("PS", value.toString) }
        val attrs         = genotype.attrs.filterNot(_._1 == "PS") ++ maybePhaseSet
        sample -> genotype.copy(attrs=attrs)
      }
      variant.copy(genotypes=genotypes)
    }
    val updater = new VcfPhaseSetUpdater(header=builder.header, keepOriginal=false, phaseGenotypesWithPhaseSet=true)
    compareVariants(
      actual   = actual.map { variant => updater.update(variant=variant)},
      expected = expected
    )
  }

  "FixVcfPhaseSet" should "run end-to-end" in {
    val (builder: VcfBuilder, expected: Seq[Variant]) = build()
    val output = makeTempFile("out.", ".vcf.gz")
    val tool   = new FixVcfPhaseSet(input=builder.toTempFile(), output=output, keepOriginal=false, phaseGenotypesWithPhaseSet=true)
    executeFgbioTool(tool)

    val actual = readVcfRecs(vcf=output)
    compareVariants(
      actual   = actual,
      expected = expected
    )
  }
}
