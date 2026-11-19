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

import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, UnitSpec, VcfBuilder}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.api.VcfContigHeader

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

  val dict = SequenceDictionary.extract(ref)

  /** A builder for a single-sample VCF with the contigs in the given sequence dictionary. */
  private def newBuilder(dict: SequenceDictionary = this.dict): VcfBuilder = {
    val contigs = VcfContigHeader.fromSequenceDictionary(dict)
    VcfBuilder(VcfBuilder.DefaultHeader.copy(contigs=contigs, samples=IndexedSeq("Sample")))
  }

  "VariantMask" should "mask SNPs as individual bases" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","T"))
    builder.add(chrom=dict(1).name, pos=102, alleles=Seq("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1, 99)  shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe false
    mask.isVariant(1, 102) shouldBe true
    mask.isVariant(1, 103) shouldBe false
  }

  it should "mask all deleted bases for deletions, plus the upstream base" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("AA","A"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1,  99) shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe true
    mask.isVariant(1, 102) shouldBe false
  }

  it should "mask just the upstream base for insertions" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","AA"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1,  99) shouldBe false
    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 101) shouldBe false
    mask.isVariant(1, 102) shouldBe false
  }

  it should "allow querying be sequence name as well as ref index" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant("chr0", 100) shouldBe false
    mask.isVariant("chr1", 100) shouldBe true
    mask.isVariant("chr2", 100) shouldBe false
    mask.isVariant("chr3", 100) shouldBe false
    mask.isVariant("chr4", 100) shouldBe false
  }

  it should "construct a mask ok from a VCF path" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","C"), gts=Seq(Gt(sample="Sample", gt="A/C")))
    val mask = VariantMask(builder.toTempFile())
    mask.isVariant(1, 100) shouldBe true
  }

  it should "construct a mask from a VCF path with single breakend alleles" in {
    val vcf = makeTempFile("VariantMaskTest.", ".vcf")
    Io.writeLines(vcf, Seq("##fileformat=VCFv4.3") ++
      dict.iterator.map(s => s"##contig=<ID=${s.name},length=${s.length}>") ++
      Seq(
        Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO").mkString("\t"),
        Seq(dict(1).name, "100", "bnd1", "A", "A.", ".", "PASS", ".").mkString("\t"),
        Seq(dict(1).name, "200", "bnd2", "A", ".A", ".", "PASS", ".").mkString("\t"),
        Seq(dict(1).name, "300", ".", "A", "T", ".", "PASS", ".").mkString("\t")
      )
    )
    val mask = VariantMask(vcf)
    Seq(99, 101, 199, 201, 299, 301).foreach(pos => mask.isVariant(1, pos) shouldBe false)
    Seq(100, 200, 300).foreach(pos => mask.isVariant(1, pos) shouldBe true)
  }

  it should "throw an exception if a VCF doesn't have a sequence dictionary in it" in {
    val builder = newBuilder(SequenceDictionary())
    val vcf = builder.toTempFile()
    an[Exception] shouldBe thrownBy { VariantMask(vcf) }
  }

  it should "throw an exception if requested to traverse backwards to an earlier reference" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","T"))
    builder.add(chrom=dict(2).name, pos=200, alleles=Seq("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(2, 200) shouldBe true
    an[Exception] shouldBe thrownBy { mask.isVariant(1, 100) }
  }

  it should "throw an exception if invalid reference sequences or positions are requested" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","T"))
    val mask = VariantMask(builder.iterator, dict)

    an[Exception] shouldBe thrownBy { mask.isVariant(-1, 100) }        // invalid index (low)
    an[Exception] shouldBe thrownBy { mask.isVariant("chrNope", 100) } // invalid ref name
    an[Exception] shouldBe thrownBy { mask.isVariant(0, -1) }          // invalid position (low)
    an[Exception] shouldBe thrownBy { mask.isVariant(0, 1e6.toInt) }   // invalid position (high)
    an[Exception] shouldBe thrownBy { mask.isVariant(9, 100) }         // invalid index (high)
  }

  it should "not mask filtered variants" in {
    val builder = newBuilder()
    builder.add(chrom=dict(1).name, pos=100, alleles=Seq("A","T"), filters=Seq("PASS"))
    builder.add(chrom=dict(1).name, pos=102, alleles=Seq("A","T"), filters=Seq("LowQD"))
    builder.add(chrom=dict(1).name, pos=104, alleles=Seq("A","T"))
    val mask = VariantMask(builder.toTempFile())

    mask.isVariant(1, 100) shouldBe true
    mask.isVariant(1, 102) shouldBe false
    mask.isVariant(1, 104) shouldBe true
  }
}
