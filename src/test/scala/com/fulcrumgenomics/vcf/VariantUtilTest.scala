/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api._
import htsjdk.variant.vcf.{VCFHeaderLineCount, VCFHeaderLineType, VCFStandardHeaderLines}

class VariantUtilTest extends UnitSpec {
  /** Builds a single variant at chr1:10 with the given alleles and genotypes. */
  private def variant(alleles: Seq[String], gts: Gt*): Variant = {
    VcfBuilder(gts.map(_.sample)).add(pos=10, alleles=alleles, gts=gts).iterator.next()
  }

  "VariantUtil.StandardFormatHeaderLines" should "match HTSJDK's definitions of the standard FORMAT fields" in {
    VariantUtil.StandardFormatHeaderLines.foreach { line =>
      val htsjdk = VCFStandardHeaderLines.getFormatLine(line.id)
      val count  = htsjdk.getCountType match {
        case VCFHeaderLineCount.INTEGER   => VcfCount.Fixed(htsjdk.getCount)
        case VCFHeaderLineCount.R         => VcfCount.OnePerAllele
        case VCFHeaderLineCount.A         => VcfCount.OnePerAltAllele
        case VCFHeaderLineCount.G         => VcfCount.OnePerGenotype
        case VCFHeaderLineCount.UNBOUNDED => VcfCount.Unknown
        case other                        => fail(s"Unexpected count type: $other")
      }
      val kind = htsjdk.getType match {
        case VCFHeaderLineType.Integer => VcfFieldType.Integer
        case VCFHeaderLineType.String  => VcfFieldType.String
        case other                     => fail(s"Unexpected type: $other")
      }
      line shouldBe VcfFormatHeader(htsjdk.getID, count, kind, htsjdk.getDescription)
    }
    VariantUtil.StandardFormatHeaderLines.map(_.id) should contain theSameElementsAs Seq("GT", "GQ", "DP", "AD", "PL", "FT")
  }

  "VariantUtil.subsetToSamples" should "keep only the genotypes of the samples, and the alleles called in them" in {
    val v = variant(Seq("A", "C", "G", "T"), Gt(sample="s1", gt="0/3"), Gt(sample="s2", gt="./."), Gt(sample="s3", gt="1/2"))
    val subset = VariantUtil.subsetToSamples(v, Set("s1", "s2"))
    subset.alleles.toSeq.map(_.value) should contain theSameElementsInOrderAs Seq("A", "T")
    subset.genotypes.keySet shouldBe Set("s1", "s2")
    subset.gt("s1").gtWithBases shouldBe "A/T"
    subset.gt("s2").isNoCall shouldBe true
  }

  it should "keep only the reference allele if the samples call no alt allele" in {
    val v = variant(Seq("A", "C"), Gt(sample="s1", gt="0/0"), Gt(sample="s2", gt="1/1"))
    VariantUtil.subsetToSamples(v, Set("s1")).alleles.toSeq shouldBe Seq(Allele("A"))
  }

  it should "keep a spanning deletion allele called in the samples" in {
    val v = variant(Seq("A", "C", "*"), Gt(sample="s1", gt="0|2"), Gt(sample="s2", gt="1/1"))
    val subset = VariantUtil.subsetToSamples(v, Set("s1"))
    subset.alleles.toSeq.map(_.value) should contain theSameElementsInOrderAs Seq("A", "*")
    subset.isSnp shouldBe true
  }

  it should "return the variant itself if it has no other samples and all of its alt alleles are called" in {
    val v = variant(Seq("A", "C", "G"), Gt(sample="s1", gt="0/1"), Gt(sample="s2", gt="1/2"))
    VariantUtil.subsetToSamples(v, Set("s1", "s2", "s3")) should be theSameInstanceAs v
  }

  "VariantUtil.toDoubles" should "convert numeric values to Doubles" in {
    VariantUtil.toDoubles(0.05f, "AF") shouldBe Seq(0.05)
    VariantUtil.toDoubles(0.25, "AF") shouldBe Seq(0.25)
    VariantUtil.toDoubles(3, "AF") shouldBe Seq(3.0)
    VariantUtil.toDoubles("0.1,0.2", "AF") shouldBe Seq(0.1, 0.2)
    VariantUtil.toDoubles(ArrayAttr(Seq(0.1f, 0.2f)), "AF") shouldBe Seq(0.1, 0.2)
    VariantUtil.toDoubles(Seq("0.3"), "AF") shouldBe Seq(0.3)
  }

  it should "throw an exception naming the field for missing or non-numeric values" in {
    Seq[Any](Variant.MissingFloat, Seq(0.1f, Variant.MissingFloat), 'c', "0.1,x").foreach { value =>
      val ex = intercept[IllegalArgumentException] { VariantUtil.toDoubles(value, "AF of sample s1 at chr1:10") }
      ex.getMessage should include ("AF of sample s1 at chr1:10")
    }
  }

  "VariantUtil.toInt" should "convert a single integer value to an Int" in {
    VariantUtil.toInt(5, "PS") shouldBe 5
    VariantUtil.toInt("5", "PS") shouldBe 5
    VariantUtil.toInt(ArrayAttr(Seq("5")), "PS") shouldBe 5
    VariantUtil.toInt(Seq(5), "PS") shouldBe 5
  }

  it should "throw an exception naming the field for missing, non-integer or multiple values" in {
    Seq[Any](Variant.MissingInt, ".", "x", "1.5", 1.5f, "5,6", ArrayAttr(Seq("5", "6")), Seq.empty[Int]).foreach { value =>
      val ex = intercept[IllegalArgumentException] { VariantUtil.toInt(value, "PS at chr1:10") }
      ex.getMessage should include ("PS at chr1:10")
    }
  }
}
