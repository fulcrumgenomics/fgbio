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

import com.fulcrumgenomics.FgBioDef.PathToVcf
import com.fulcrumgenomics.commons.CommonsDef.SafelyClosable
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import htsjdk.samtools.util.FileExtensions.{BCF => BcfExtension, COMPRESSED_VCF => VcfGzExtension, CSI => CsiExtensions, TABIX_INDEX => TbiExtension, TRIBBLE_INDEX => IdxExtension, VCF => VcfExtension}
import org.scalatest.OptionValues

import java.nio.file.{Files, Path}
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
    val vcf = makeTempFile(getClass.getSimpleName, VcfExtension)
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
    val vcf = makeTempFile("queryable.", VcfGzExtension)
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

  it should "not allow a VCF writer to write to a directory path" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples = samples)
    val output  = Io.makeTempDir(getClass.getSimpleName)
    output.toFile.deleteOnExit()
    an[IllegalArgumentException] shouldBe thrownBy { VcfWriter(output, header = builder.header) }
  }

  it should "write a sibling index when writing to a plaintext VCF file" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples = samples)
    val output  = makeTempFile(getClass.getSimpleName, VcfExtension)
    val writer  = VcfWriter(output, header = builder.header)
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    writer.close()
    val source = VcfSource(output)
    source.isQueryable shouldBe true
    source.safelyClose()
  }

  it should "write a sibling index when writing to a compressed VCF file" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples = samples)
    val output  = makeTempFile(getClass.getSimpleName, VcfGzExtension)
    val writer  = VcfWriter(output, header = builder.header)
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    writer.close()
    val source = VcfSource(output)
    source.isQueryable shouldBe true
    source.safelyClose()
  }

  it should "write a sibling index when writing to a BCF file" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples = samples)
    val output  = makeTempFile(getClass.getSimpleName, BcfExtension)
    val writer  = VcfWriter(output, header = builder.header)
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    writer.close()
    val source = VcfSource(output)
    source.isQueryable shouldBe true
    source.safelyClose()
  }

  it should "not write a sibling index when asked not to index" in {
    val builder = VcfBuilder(samples = Seq("sample"))
    val output  = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, header = builder.header, index = false).close()
    val source = VcfSource(output)
    source.isQueryable shouldBe false
    source.safelyClose()
  }

  it should "write variants out of coordinate order when asked not to index" in {
    val builder = VcfBuilder(samples = Seq("sample"))
    builder.add(chrom = "chr1", pos = 200, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    val writer = VcfWriter(output, header = builder.header, index = false)
    writer.write(builder.toSeq.sortBy(-_.pos))
    writer.close()
    val source = VcfSource(output)
    source.map(_.pos).toList shouldBe List(200, 100)
    source.safelyClose()
  }

  it should "not attempt to write a sibling index when streaming to a named pipe like '/dev/null'" in {
    val samples = Seq("sample")
    val builder = VcfBuilder(samples = samples)
    val writer  = VcfWriter(Io.DevNull, header = builder.header)
    builder.add(chrom = "chr1", pos = 100, alleles = Seq("A", "C"), gts = Seq(Gt(sample = "sample", gt = "0/1")))
    noException shouldBe thrownBy { writer.write(builder.toSeq); writer.close() }
    Files.exists(PathUtil.pathTo(Io.DevNull.getFileName.toString + CsiExtensions)) shouldBe false
    Files.exists(PathUtil.pathTo(Io.DevNull.getFileName.toString + IdxExtension))  shouldBe false
    Files.exists(PathUtil.pathTo(Io.DevNull.getFileName.toString + TbiExtension))  shouldBe false
  }

  /** Writes a header-only VCF 4.5 whose FORMAT lines use each of the VCF 4.4 and 4.5 `Number=` codes. */
  private def writeVcf45WithNewNumberCodes(): Path = {
    val vcf = makeTempFile(getClass.getSimpleName, VcfExtension)
    Io.writeLines(vcf, Seq(
      "##fileformat=VCFv4.5",
      "##contig=<ID=chr1,length=5000>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##FORMAT=<ID=XP,Number=P,Type=Integer,Description=\"One per GT allele\">",
      "##FORMAT=<ID=XLA,Number=LA,Type=Integer,Description=\"One per local alt allele\">",
      "##FORMAT=<ID=XLR,Number=LR,Type=Integer,Description=\"One per local allele\">",
      "##FORMAT=<ID=XLG,Number=LG,Type=Integer,Description=\"One per local genotype\">",
      "##FORMAT=<ID=XM,Number=M,Type=Integer,Description=\"One per base modification\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"
    ))
    vcf
  }

  /** The counts of the FORMAT lines written by [[writeVcf45WithNewNumberCodes()]], keyed by ID. */
  private val NewNumberCodeCounts: Map[String, VcfCount] = Map(
    "XP"  -> VcfCount.OnePerGtAllele,
    "XLA" -> VcfCount.OnePerLocalAltAllele,
    "XLR" -> VcfCount.OnePerLocalAllele,
    "XLG" -> VcfCount.OnePerLocalGenotype,
    "XM"  -> VcfCount.OnePerBaseModification
  )

  it should "read FORMAT lines that use the VCF 4.4 and 4.5 Number codes" in {
    val source = VcfSource(writeVcf45WithNewNumberCodes())
    val counts = source.header.formats.map(f => f.id -> f.count).toMap
    source.safelyClose()
    NewNumberCodeCounts.foreach { case (id, count) => counts(id) shouldBe count }
  }

  it should "write back FORMAT lines that use the VCF 4.4 and 4.5 Number codes" in {
    val source = VcfSource(writeVcf45WithNewNumberCodes())
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, source.header).close()
    source.safelyClose()

    val reread = VcfSource(output)
    val counts = reread.header.formats.map(f => f.id -> f.count).toMap
    reread.safelyClose()
    NewNumberCodeCounts.foreach { case (id, count) => counts(id) shouldBe count }
  }

  it should "keep the version a VCF header declares when reading" in {
    val source = VcfSource(writeVcf45WithNewNumberCodes())
    source.header.version.value shouldBe VcfVersion.V4_5
    source.safelyClose()
  }

  it should "write a VCF at the version its header declares" in {
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, VcfIoTest.Header.copy(version = Some(VcfVersion.V4_3))).close()
    Io.readLines(output).next() shouldBe "##fileformat=VCFv4.3"
  }

  it should "write a VCF as VCFv4.2 when its header declares no version" in {
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, VcfIoTest.Header.copy(version = None)).close()
    Io.readLines(output).next() shouldBe "##fileformat=VCFv4.2"
  }

  it should "write a VCF at the version its lines need when its header declares no version" in {
    val local  = VcfFormatHeader(id="XLA", count=VcfCount.OnePerLocalAltAllele, kind=VcfFieldType.Integer, description="x")
    val header = VcfIoTest.Header.copy(formats = VcfIoTest.Header.formats :+ local, version = None)
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, header).close()
    Io.readLines(output).next() shouldBe "##fileformat=VCFv4.5"
  }

  it should "write a BCF at the version its header declares" in {
    val source = VcfSource(writeVcf45WithNewNumberCodes())
    val output = makeTempFile(getClass.getSimpleName, BcfExtension)
    VcfWriter(output, source.header).close()
    source.safelyClose()

    val reread = VcfSource(output)
    reread.header.version.value shouldBe VcfVersion.V4_5
    reread.header.formats.find(_.id == "XLA").value.count shouldBe VcfCount.OnePerLocalAltAllele
    reread.safelyClose()
  }

  /** The header lines of the VCFs written by [[writeRawVcf]]. */
  private val RawHeaderLines: Seq[String] = Seq(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=5000>",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">",
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"mq\">",
    "##INFO=<ID=S1,Number=1,Type=String,Description=\"s\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">",
    "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">",
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
    "##FORMAT=<ID=XF,Number=2,Type=Float,Description=\"xf\">",
    "##FORMAT=<ID=XI,Number=2,Type=Integer,Description=\"xi\">",
    Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1", "s2").mkString("\t")
  )

  /** Writes a two-sample VCF with the given records, each given as its tab-separated fields. */
  private def writeRawVcf(records: Seq[String]*): PathToVcf = {
    val vcf = makeTempFile(getClass.getSimpleName, VcfExtension)
    Io.writeLines(vcf, RawHeaderLines ++ records.map(_.mkString("\t")))
    vcf
  }

  /** Reads the VCF or BCF with a VcfSource and writes it back out to the given path with a VcfWriter. */
  private def copyVcf(input: PathToVcf, output: PathToVcf, writeFullFormatField: Boolean = false): Unit = {
    val source = VcfSource(input)
    val writer = VcfWriter(output, source.header, writeFullFormatField=writeFullFormatField)
    writer ++= source
    writer.close()
    source.safelyClose()
  }

  /** Reads the VCF or BCF with a VcfSource and writes it back out as a VCF with a VcfWriter, returning the records
    * written as text. */
  private def roundTripText(input: PathToVcf, writeFullFormatField: Boolean = false): Seq[String] = {
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    copyVcf(input, output, writeFullFormatField=writeFullFormatField)
    Io.readLines(output).filterNot(_.startsWith("#")).toSeq
  }

  "VcfSource and VcfWriter" should "keep missing values in multi-valued INFO and FORMAT fields" in {
    val vcf = writeRawVcf(
      Seq("chr1", "100", ".", "A", "C,G", "50", "PASS", "AC=.,3;AF=0.5,.", "GT:XF:XI", "0/1:1.5,.:.,7", "0/0:.,2.5:3,.")
    )

    val variant = VcfSource(vcf).iterator.next()
    variant[IndexedSeq[Int]]("AC") shouldBe Seq(Variant.MissingInt, 3)
    Variant.isMissingValue(variant[IndexedSeq[Float]]("AF").last) shouldBe true

    roundTripText(vcf) shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C,G", "50", "PASS", "AC=.,3;AF=0.5,.", "GT:XF:XI", "0/1:1.5,.:.,7", "0/0:.,2.5:3,.").mkString("\t")
    )
  }

  it should "write a missing Int or Float on its own as a missing value" in {
    val alleles = AlleleSet(ref=Allele("A"), alts=IndexedSeq(Allele("C")))
    val variant = Variant(chrom="chr1", pos=100, alleles=alleles, attrs=ListMap("DP" -> Variant.MissingInt, "QD" -> Variant.MissingFloat),
      genotypes=Map("s1" -> Genotype(alleles, "s1", alleles.toIndexedSeq, attrs=Map("GQ" -> Variant.MissingInt, "AD" -> Seq(1, 2)))))
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    val writer = VcfWriter(output, VcfIoTest.Header)
    writer.write(variant)
    writer.close()

    Io.readLines(output).filterNot(_.startsWith("#")).toSeq shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C", ".", ".", "DP=.;QD=.", "GT:AD", "0/1:1,2").mkString("\t")
    )
  }

  it should "read and write the genotype filters (FT)" in {
    val vcf = writeRawVcf(Seq("chr1", "100", ".", "A", "C", "50", "PASS", ".", "GT:DP:FT", "0/1:20:lowq", "0/1:30:a;b"))

    val variant = VcfSource(vcf).iterator.next()
    variant.gt("s1")[String]("FT") shouldBe "lowq"
    variant.gt("s2")[String]("FT") shouldBe "a;b"

    roundTripText(vcf) shouldBe Seq(Seq("chr1", "100", ".", "A", "C", "50", "PASS", ".", "GT:DP:FT", "0/1:20:lowq", "0/1:30:a;b").mkString("\t"))
  }

  it should "write genotype filters (FT) given as a sequence" in {
    val alleles = AlleleSet(ref=Allele("A"), alts=IndexedSeq(Allele("C")))
    val variant = Variant(chrom="chr1", pos=100, alleles=alleles,
      genotypes=Map("s1" -> Genotype(alleles, "s1", alleles.toIndexedSeq, attrs=Map("FT" -> Seq("lowq", "lowdp")))))
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    val writer = VcfWriter(output, VcfIoTest.Header)
    writer.write(variant)
    writer.close()

    Io.readLines(output).filterNot(_.startsWith("#")).toSeq shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C", ".", ".", ".", "GT:FT", "0/1:lowdp;lowq").mkString("\t")
    )
  }

  it should "keep the commas in a Number=1 String INFO value" in {
    val vcf = writeRawVcf(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "S1=a,b", "GT", "0/1", "0/0"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "S1=.,.", "GT", "0/1", "0/0")
    )
    VcfSource(vcf).toIndexedSeq.map(_[String]("S1")) shouldBe Seq("a,b", ".,.")
    roundTripText(vcf) shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "S1=a,b", "GT", "0/1", "0/0").mkString("\t"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "S1=.,.", "GT", "0/1", "0/0").mkString("\t")
    )
  }

  it should "write a BCF to a path whose extension is upper case" in {
    val vcf = writeRawVcf(Seq("chr1", "100", ".", "A", "C,G", "50", "PASS", "AF=0.5,0.25", "GT:XF", "0/1:1.5,.", "0/0:0.25,2"))
    val bcf = makeTempFile(getClass.getSimpleName, ".BCF")
    copyVcf(vcf, bcf)
    roundTripText(bcf) shouldBe roundTripText(vcf)
  }

  it should "write Float values as the shortest decimal that reads back as the same Float" in {
    val vcf = writeRawVcf(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "AF=0.500;MQ=60", "GT:XF", "0/1:0,1e-10", "0/0:12345678.9,-2.50"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "AF=0.0001;MQ=59.123456789", "GT", "0/1", "0/0")
    )
    roundTripText(vcf) shouldBe Seq(
      // NB: a Float holds only about 7 significant digits, so 12345678.9 and 59.123456789 cannot be kept exactly
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "AF=0.5;MQ=60", "GT:XF", "0/1:0,0.0000000001", "0/0:12345679,-2.5").mkString("\t"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "AF=0.0001;MQ=59.123455", "GT", "0/1", "0/0").mkString("\t")
    )
  }

  it should "read NaN and infinite Float values written in any case, with or without a sign" in {
    val vcf = writeRawVcf(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "MQ=nan", "GT:XF", "0/1:Inf,-INF", "0/0:+infinity,NaN"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "MQ=-Infinity", "GT", "0/1", "0/0")
    )
    val Seq(v1, v2) = VcfSource(vcf).toIndexedSeq
    v1[Float]("MQ").isNaN shouldBe true
    v1.gt("s1")[IndexedSeq[Float]]("XF") shouldBe Seq(Float.PositiveInfinity, Float.NegativeInfinity)
    v1.gt("s2")[IndexedSeq[Float]]("XF").head shouldBe Float.PositiveInfinity
    v2[Float]("MQ") shouldBe Float.NegativeInfinity

    roundTripText(vcf) shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "MQ=NaN", "GT:XF", "0/1:Inf,-Inf", "0/0:Inf,NaN").mkString("\t"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "MQ=-Inf", "GT", "0/1", "0/0").mkString("\t")
    )
  }

  it should "leave out trailing missing FORMAT fields unless asked to write them all" in {
    val vcf = writeRawVcf(Seq("chr1", "100", ".", "A", "C", "50", "PASS", ".", "GT:DP:GQ:PS", "0/1:10:.:.", "0|1:12:30:100"))
    roundTripText(vcf) shouldBe Seq(Seq("chr1", "100", ".", "A", "C", "50", "PASS", ".", "GT:DP:GQ:PS", "0/1:10", "0|1:12:30:100").mkString("\t"))
    roundTripText(vcf, writeFullFormatField=true) shouldBe Seq(
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", ".", "GT:DP:GQ:PS", "0/1:10:.:.", "0|1:12:30:100").mkString("\t")
    )
  }

  it should "round-trip INFO and FORMAT values, including missing values, through a BCF" in {
    val records = Seq(
      Seq("chr1", "100", ".", "A", "C,G", "50", "PASS", "AC=.,3;AF=0.5,.;MQ=59.5;S1=a,b", "GT:DP:FT:XF:XI", "0/1:20:lowq:1.5,.:.,7", "0/0:3:.:0.25,2:3,4"),
      Seq("chr1", "200", ".", "A", "C", "50", "PASS", "AC=1;AF=0.0001;MQ=60", "GT:DP", "0/1:10", "0/0:12")
    )
    val vcf = writeRawVcf(records:_*)
    val bcf = makeTempFile(getClass.getSimpleName, BcfExtension)
    copyVcf(vcf, bcf)

    val fromVcf = VcfSource(vcf).toIndexedSeq
    val fromBcf = VcfSource(bcf).toIndexedSeq
    fromBcf.map(_.attrs.map { case (k, v) => k -> v.toString }) shouldBe fromVcf.map(_.attrs.map { case (k, v) => k -> v.toString })
    fromBcf.map(_.gts.map(_.attrs.map { case (k, v) => k -> v.toString }).toSeq) shouldBe fromVcf.map(_.gts.map(_.attrs.map { case (k, v) => k -> v.toString }).toSeq)

    roundTripText(bcf) shouldBe roundTripText(vcf)
  }

  it should "read and write single breakend alleles, through a VCF and a BCF" in {
    val records = Seq(
      Seq("chr1", "100", "bnd1", "C", "C.", "50", "PASS", ".", "GT", "0/1", "0/0"),
      Seq("chr1", "200", "bnd2", "T", ".T,G]chr1:300]", "50", "PASS", ".", "GT", "1/2", "0/1")
    )
    val vcf = writeRawVcf(records:_*)

    val variants = VcfSource(vcf).toIndexedSeq
    variants.map(_.alleles.alts.map(_.value)) shouldBe Seq(Seq("C."), Seq(".T", "G]chr1:300]"))
    variants.flatMap(_.alleles.alts).foreach(_.isInstanceOf[Allele.BreakendAllele] shouldBe true)
    variants.map(_.gt("s1").gtWithBases) shouldBe Seq("C/C.", ".T/G]chr1:300]")
    variants.map(_.end) shouldBe Seq(100, 200)

    // VCF -> VCF
    roundTripText(vcf) shouldBe records.map(_.mkString("\t"))

    // VCF -> BCF -> VCF
    val bcf = makeTempFile(getClass.getSimpleName, BcfExtension)
    copyVcf(vcf, bcf)
    VcfSource(bcf).toIndexedSeq.map(_.alleles.toSeq.map(_.value)) shouldBe variants.map(_.alleles.toSeq.map(_.value))
    roundTripText(bcf) shouldBe records.map(_.mkString("\t"))
  }

  it should "fail with a message naming the VCF, record, field, type and value when a typed value is invalid" in {
    case class Case(info: String, format: String, s1: String, expected: String)
    Seq(
      Case("AC=3000000000", "GT",    "0/1",          "Invalid Integer value '3000000000' for INFO field AC at chr1:100"),
      Case("AC=1,1.5",      "GT",    "0/1",          "Invalid Integer value '1.5' for INFO field AC at chr1:100"),
      Case("MQ=high",       "GT",    "0/1",          "Invalid Float value 'high' for INFO field MQ at chr1:100"),
      Case(".",             "GT:PS", "0/1:1.5",      "Invalid Integer value '1.5' for FORMAT field PS of sample s1 at chr1:100"),
      Case(".",             "GT:XI", "0/1:1,99999999999", "Invalid Integer value '99999999999' for FORMAT field XI of sample s1 at chr1:100")
    ).foreach { c =>
      val vcf = writeRawVcf(Seq("chr1", "100", ".", "A", "C", "50", "PASS", c.info, c.format, c.s1, "0/0" + ":." * (c.format.count(_ == ':'))))
      val source = VcfSource(vcf)
      val ex = intercept[IllegalArgumentException] { source.toIndexedSeq }
      source.safelyClose()
      ex.getMessage shouldBe s"${c.expected} in $vcf"
    }
  }

  "VcfConversions.formatFloat" should "format Floats as plain, shortest round-trip decimals" in {
    VcfConversions.formatFloat(60f) shouldBe "60"
    VcfConversions.formatFloat(0.5f) shouldBe "0.5"
    VcfConversions.formatFloat(0.1f) shouldBe "0.1"
    VcfConversions.formatFloat(1e-10f) shouldBe "0.0000000001"
    VcfConversions.formatFloat(1.2345679e7f) shouldBe "12345679"
    VcfConversions.formatFloat(-2.5f) shouldBe "-2.5"
    VcfConversions.formatFloat(0f) shouldBe "0"
    VcfConversions.formatFloat(-0f) shouldBe "-0"
    VcfConversions.formatFloat(Float.NaN) shouldBe "NaN"
    VcfConversions.formatFloat(Float.PositiveInfinity) shouldBe "Inf"
    VcfConversions.formatFloat(Float.NegativeInfinity) shouldBe "-Inf"
    Seq(0.1f, 1/3f, 123.456f, 3.4028235e38f, 1.4e-45f).foreach { f => VcfConversions.formatFloat(f).toFloat shouldBe f }
  }

  "VcfSource and VcfWriter" should "round-trip the length, assembly, MD5, species and URL of contig lines" in {
    val contig = VcfContigHeader(index=0, name="chr1", length=Some(5000), assembly=Some("hg38"),
      md5=Some("6aef897c3d6ff0c78aff06ac189178dd"), species=Some("Homo sapiens"), url=Some("file:/ref.fa"))
    val header = VcfIoTest.Header.copy(contigs=IndexedSeq(contig))
    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    VcfWriter(output, header).close()

    Io.readLines(output).filter(_.startsWith("##contig")).toSeq shouldBe Seq(
      "##contig=<ID=chr1,length=5000,assembly=hg38,md5=6aef897c3d6ff0c78aff06ac189178dd,species=\"Homo sapiens\",URL=file:/ref.fa>"
    )

    val source = VcfSource(output)
    source.header.contigs should contain theSameElementsInOrderAs Seq(contig)
    source.header.dict(0).md5 shouldBe contig.md5
    source.header.dict(0).assembly shouldBe contig.assembly
    source.header.dict(0).species shouldBe contig.species
    source.header.dict(0).uri shouldBe contig.url
    source.safelyClose()
  }

  "VcfContigHeader.fromSequenceMetadata" should "keep the length, assembly, MD5, species and URL of the sequence" in {
    val metadata = SequenceMetadata(name="chr1", length=5000, index=0, assembly=Some("hg38"), md5=Some("abc"),
      species=Some("Homo sapiens"), uri=Some("file:/ref.fa"))
    VcfContigHeader.fromSequenceMetadata(metadata) shouldBe VcfContigHeader(index=0, name="chr1", length=Some(5000),
      assembly=Some("hg38"), md5=Some("abc"), species=Some("Homo sapiens"), url=Some("file:/ref.fa"))
  }

  it should "leave out the length of a sequence of unknown length" in {
    VcfContigHeader.fromSequenceMetadata(SequenceMetadata(name="chr1", index=0)).length shouldBe None
  }

  /** Writes a VCF whose single variant has INFO and FORMAT fields that are not described in the header. */
  private def writeVcfWithUndescribedFields(): PathToVcf = {
    val vcf = makeTempFile(getClass.getSimpleName, VcfExtension)
    Io.writeLines(vcf, Seq(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=5000>",
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1").mkString("\t"),
      Seq("chr1", "100", ".", "A", "C", "50", "PASS", "DP=10;XF;XL=a,b;XS=0.5", "GT:XG:XH", "0/1:7:x,y").mkString("\t")
    ))
    vcf
  }

  "VcfSource" should "fail to read fields not described in the header by default" in {
    val source = VcfSource(writeVcfWithUndescribedFields())
    val ex     = intercept[RuntimeException] { source.toIndexedSeq }
    Iterator.iterate[Throwable](ex)(_.getCause).takeWhile(_ != null).map(_.getMessage).mkString(" ") should include ("not described in header")
    source.safelyClose()
  }

  it should "read an END INFO field not described in the header as an Int when allowed" in {
    val vcf = makeTempFile(getClass.getSimpleName, VcfExtension)
    Io.writeLines(vcf, Seq(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=5000>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1").mkString("\t"),
      Seq("chr1", "100", ".", "A", "<DEL>", "50", "PASS", "END=150", "GT", "0/1").mkString("\t")
    ))
    val source  = VcfSource(vcf, allowMissingFieldsInHeader=true)
    val variant = source.iterator.next()
    variant[Int]("END") shouldBe 150
    variant.end shouldBe 150
    source.safelyClose()
  }

  it should "read fields not described in the header when allowed, and write them back as they were (with INFO fields sorted)" in {
    val input   = writeVcfWithUndescribedFields()
    val source  = VcfSource(input, allowMissingFieldsInHeader=true)
    val variant = source.iterator.next()
    variant[Int]("DP") shouldBe 10
    variant[Seq[String]]("XL") should contain theSameElementsInOrderAs Seq("a", "b")
    variant[Seq[String]]("XS") should contain theSameElementsInOrderAs Seq("0.5")
    variant[String]("XF") shouldBe Variant.FlagValue
    variant.gt("s1")[Seq[String]]("XG") should contain theSameElementsInOrderAs Seq("7")
    variant.gt("s1")[Seq[String]]("XH") should contain theSameElementsInOrderAs Seq("x", "y")

    val output = makeTempFile(getClass.getSimpleName, VcfExtension)
    val writer = VcfWriter(output, source.header)
    writer.write(variant)
    writer.close()
    source.safelyClose()

    Io.readLines(output).filterNot(_.startsWith("#")).toSeq shouldBe Io.readLines(input).filterNot(_.startsWith("#")).toSeq
  }
}
