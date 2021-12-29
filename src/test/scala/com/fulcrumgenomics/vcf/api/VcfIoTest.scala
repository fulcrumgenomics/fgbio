/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import org.scalatest.OptionValues

import scala.collection.compat._
import scala.collection.immutable.ListMap

object VcfIoTest {
  val Header = VcfHeader(
    contigs = IndexedSeq(VcfContigHeader(index=0, name="chr1", length=Some(5000))),
    infos   = Seq(
      VcfInfoHeader(id="AC",  count=VcfCount.OnePerAltAllele, kind=VcfFieldType.Integer, description="Alternate allele counts in genotypes."),
      VcfInfoHeader(id="DP",  count=VcfCount.Fixed(1),        kind=VcfFieldType.Integer, description="Depth across all samples."),
      VcfInfoHeader(id="QD",  count=VcfCount.Fixed(1),        kind=VcfFieldType.Float,   description="Quality/depth."),
      VcfInfoHeader(id="STR", count=VcfCount.Fixed(0),        kind=VcfFieldType.Flag,    description="Quality/depth.")
    ),
    formats = Seq(
      VcfFormatHeader(id="GT", count=VcfCount.Fixed(1),       kind=VcfFieldType.String, description="Genotype string."),
      VcfFormatHeader(id="AD", count=VcfCount.OnePerAllele,   kind=VcfFieldType.Integer, description="Depth per allele."),
      VcfFormatHeader(id="GQ", count=VcfCount.Fixed(1),       kind=VcfFieldType.Integer, description="Genotype quality."),
      VcfFormatHeader(id="PL", count=VcfCount.OnePerGenotype, kind=VcfFieldType.Integer, description="Phred scaled genotype likelihoods.")
    ),
    filters = Seq(
      VcfFilterHeader(id="LowQD", description="Low Quality/Depth value"),
      VcfFilterHeader(id="LowAB", description="Low/poor allele balance.")
    ),
    others  = Seq(
      VcfGeneralHeader(headerType="simple",   id="is as simple does"),
      VcfGeneralHeader(headerType="compound", id="time-of-day", data=Map("hour" -> "13", "minute" -> "53")),
      VcfGeneralHeader(headerType="ALT",      id="NON_REF", data=Map("Description" -> "Represents any non-reference allele."))
    ),
    samples = IndexedSeq("s1")
  )
}

/** Tests for reading and writing VCFs using the scala API. */
class VcfIoTest extends UnitSpec with OptionValues {
  case class Result(header: VcfHeader, vs: IndexedSeq[Variant])

  /** Writes out the variants to a file and reads them back, returning the header and the variants. */
  @inline private def roundtrip(variants: IterableOnce[Variant], header: VcfHeader = VcfIoTest.Header): Result = {
    val vcf = makeTempFile("test.", ".vcf")
    val out = VcfWriter(vcf, header)
    out ++= variants
    out.close()
    val in = VcfSource(vcf)
    val results = Result(in.header, in.toIndexedSeq)
    in.close()

    // Check that the header didn't get mangled other than perhaps having lines resorted
    results.header.samples should contain theSameElementsInOrderAs header.samples
    results.header.contigs should contain theSameElementsInOrderAs header.contigs
    results.header.infos   should contain theSameElementsAs        header.infos
    results.header.formats should contain theSameElementsAs        header.formats
    results.header.filters should contain theSameElementsAs        header.filters
    results.header.others  should contain theSameElementsAs        header.others

    results
  }

  "VcfReader and VcfWriter" should "write a VCF and read back a VCF" in {
    val alleles = AlleleSet(ref=Allele("C"), alts=IndexedSeq(Allele("G"), Allele("T")))
    val variant = Variant(
      chrom     = "chr1",
      pos       = 1000,
      id        = Some("testvar"),
      alleles   = alleles,
      qual      = Some(1234.5),
      attrs     = ListMap("DP" -> 120, "AC" -> Seq(1, 1), "STR" -> "."),
      genotypes = Map("s1" -> Genotype(alleles, "s1", alleles.alts, phased=false, attrs=Map("AD" -> Seq(0, 50, 54))))
    )

    val Result(header, variants) = roundtrip(Seq(variant))
    header.samples should contain theSameElementsInOrderAs Seq("s1")
    header.dict(0).name shouldBe "chr1"

    variants should have size 1

    variants.head.chrom                            shouldBe "chr1"
    variants.head.pos                              shouldBe 1000
    variants.head.end                              shouldBe 1000
    variants.head.alleles.ref.toString             shouldBe "C"
    variants.head.alleles.alts.map(_.toString)     should contain theSameElementsInOrderAs Seq("G", "T")
    variants.head.qual.value                       shouldBe 1234.5
    variants.head.apply[Int]("DP")                 shouldBe 120
    variants.head.apply[ArrayAttr[Int]]("AC")      should contain theSameElementsInOrderAs Seq(1, 1)
    variants.head.attrs.contains("STR")            shouldBe true
    variants.head.attrs.contains("QD")             shouldBe false

    variants.head.genotypes.size           shouldBe 1
    val gt = variants.head.genotypes("s1")
    gt.gtWithBases                         shouldBe "G/T"
    gt.gtVcfStyle                          shouldBe "1/2"
    gt[ArrayAttr[Int]]("AD")               should contain theSameElementsInOrderAs Seq(0, 50, 54)
    gt.attrs.contains("PL")                shouldBe false
  }


  it should "write and read back a multi-sample VCF with some missing fields for the second sample" in {
    val alleles = AlleleSet(ref=Allele("C"), alts=IndexedSeq(Allele("G")))
    val gts = Map(
      "s1" -> Genotype(alleles, "s1", calls=alleles.toIndexedSeq, phased=false, attrs=Map("AD" -> Seq(50, 54), "GQ" -> 99, "PL" -> Seq(999, 41, 998))),
      "s2" -> Genotype(alleles, "s2", calls=alleles.alts ++ alleles.alts, phased=false, attrs=Map("AD" -> Seq(0, 101)))
    )
    val inVariants = Seq(Variant(chrom="chr1", pos=1200, alleles=alleles, attrs=ListMap("DP" -> 205, "AC" -> Seq(1, 3)), genotypes=gts))
    val inHeader   = VcfIoTest.Header.copy(samples=IndexedSeq("s1", "s2"))

    val result = roundtrip(variants=inVariants, header=inHeader)
    result.header.samples should contain theSameElementsInOrderAs Seq("s1", "s2")

    result.vs should have size 1
    result.vs.head.chrom                            shouldBe "chr1"
    result.vs.head.pos                              shouldBe 1200
    result.vs.head.end                              shouldBe 1200
    result.vs.head.alleles.ref.toString             shouldBe "C"
    result.vs.head.alleles.alts.map(_.toString)     should contain theSameElementsInOrderAs Seq("G")
    result.vs.head.qual                             shouldBe None
    result.vs.head.apply[Int]("DP")                 shouldBe 205
    result.vs.head.apply[ArrayAttr[Int]]("AC")      should contain theSameElementsInOrderAs Seq(1, 3)

    result.vs.head.genotypes.size shouldBe 2
    val Seq(s1, s2): Seq[Genotype] = Seq("s1", "s2").map(result.vs.head.genotypes.apply)

    s1.gtWithBases                shouldBe "C/G"
    s1.gtVcfStyle                 shouldBe "0/1"
    s1[ArrayAttr[Int]]("AD")      should contain theSameElementsInOrderAs Seq(50, 54)
    s1[ArrayAttr[Int]]("PL")      should contain theSameElementsInOrderAs Seq(999, 41, 998)
    s1[Int]("GQ")                 shouldBe 99

    s2.gtWithBases                shouldBe "G/G"
    s2.gtVcfStyle                 shouldBe "1/1"
    s2[ArrayAttr[Int]]("AD")      should contain theSameElementsInOrderAs Seq(0, 101)
    s2.attrs.contains("PL")       shouldBe false
    s2.attrs.contains("GQ")       shouldBe false
  }

  it should "write a VCF with an index and perform queries properly" in {
    val alleles = AlleleSet(ref=Allele("C"), alts=IndexedSeq(Allele("G"), Allele("T")))
    val variants = Range.inclusive(1000, 2000, step=10).map { s =>
      Variant(chrom="chr1", pos=s, alleles=alleles, genotypes=Map("s1" -> Genotype(alleles, "s1", alleles.alts)))
    }
    val vcf = makeTempFile("queryable.", ".vcf.gz")
    val out = VcfWriter(vcf, VcfIoTest.Header)
    out ++= variants
    out.close()

    val in = VcfSource(vcf)
    in.query("chr1", 1000, 1020) should have size 3
    in.query("chr1", 1500, 1599) should have size 10
    in.query("chr1", 999,  2001) should have size 101
    in.iterator should have size 101
  }

  it should "round-trip a VCF that has no-calls in the genotypes" in {
    val alleles = AlleleSet(ref=Allele("C"), alts=IndexedSeq(Allele("G"), Allele("T")))
    val variants = Seq(
      Variant(chrom="chr1", pos=10, alleles=alleles, genotypes=Map("s1" -> Genotype(alleles, "s1", IndexedSeq(NoCallAllele, NoCallAllele)))),
      Variant(chrom="chr1", pos=20, alleles=alleles, genotypes=Map("s1" -> Genotype(alleles, "s1", IndexedSeq(Allele("G"),  NoCallAllele))))
    )

    val out = roundtrip(variants)
    out.vs.size shouldBe 2
    out.vs.find(_.pos == 10).value.gt("s1").gtWithBases shouldBe "./."
    out.vs.find(_.pos == 20).value.gt("s1").gtWithBases shouldBe "G/."
  }

  it should "round-trip variants with a variety of filters" in {
    val vs = Seq(
      Variant(chrom="chr1", pos=10, alleles=AlleleSet("A", "C"), id=Some("PASS"),  filters=Variant.PassingFilters),
      Variant(chrom="chr1", pos=20, alleles=AlleleSet("A", "C"), id=None        ,  filters=Variant.EmptyFilters),
      Variant(chrom="chr1", pos=30, alleles=AlleleSet("A", "C"), id=Some("LowQD"), filters=Set("LowQD")),
    )
    val result = roundtrip(vs, header=VcfIoTest.Header.copy(samples=IndexedSeq.empty))
    result.vs.foreach { v => v.filters.headOption shouldBe v.id }
  }

  it should "round-trip variants with many samples and not mess up the genotypes" in {
    val samples = Seq("F", "E", "G", "C", "D", "0")
    val builder = VcfBuilder(samples=samples)
    builder.add(chrom="chr1", pos=100, alleles=Seq("A", "C", "G", "T"), gts=Seq(
      Gt(sample="F", gt="0/1"),
      Gt(sample="E", gt="0/2"),
      Gt(sample="G", gt="0/3"),
      Gt(sample="C", gt="0/0"),
      Gt(sample="D", gt="1/3"),
      Gt(sample="0", gt="2/3")
    ))

    val vcf = builder.toTempFile()
    val in  = VcfSource(vcf)
    in.header.samples should contain theSameElementsInOrderAs samples
    val v   = in.iterator.next()

    v.gt("F").callIndices.mkString("/") shouldBe "0/1"
    v.gt("E").callIndices.mkString("/") shouldBe "0/2"
    v.gt("G").callIndices.mkString("/") shouldBe "0/3"
    v.gt("C").callIndices.mkString("/") shouldBe "0/0"
    v.gt("D").callIndices.mkString("/") shouldBe "1/3"
    v.gt("0").callIndices.mkString("/") shouldBe "2/3"
  }

  it should "not attempt to index a VCF when streaming to a file handle or other kind of non-regular file" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples=samples)
    val writer  = VcfWriter(Io.DevNull, header = builder.header)
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample="sample", gt="0/1")))
    noException shouldBe thrownBy { writer.write(builder.toSeq) }
    writer.close()
  }
}
