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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec, VcfBuilder}
import com.fulcrumgenomics.umi.ReviewConsensusVariants.ConsensusVariantReviewInfo
import com.fulcrumgenomics.util.{Io, Metric}
import com.fulcrumgenomics.vcf.api.{VcfContigHeader, VcfCount, VcfFieldType, VcfFormatHeader}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.IntervalList


object ReviewConsensusVariantsTest {
  val Fasta : Seq[String] =
    """
      |>chr1
      |AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      |AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      |>chr2
      |CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      |CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    """.stripMargin.linesIterator.map(_.trim).filter(_.nonEmpty).toSeq

  val Fai : Seq[String] =
    """
      |chr1    100     6       50      51
      |chr2    100     114     50      51
    """.stripMargin.linesIterator.map(_.trim.replaceAll("\\s+", "\t")).filter(_.nonEmpty).toSeq

  val Dict : Seq[String] =
    """
      |@HD    	VN:1.5 	SO:unsorted
      |@SQ    	SN:chr1	LN:100 	M5:8adc5937e635f6c9af646f0b23560fae    	AS:n/a 	UR:n/a 	SP:n/a
      |@SQ    	SN:chr2	LN:100 	M5:48ecc9d349522f836cadf615e370bc51    	AS:n/a 	UR:n/a 	SP:n/a
    """.stripMargin.linesIterator.map(_.trim.replaceAll("\\s+", "\t")).filter(_.nonEmpty).toSeq
}

class ReviewConsensusVariantsTest extends UnitSpec {
  import ReviewConsensusVariantsTest._

  def buildRefFasta(withFai: Boolean = true, withDict: Boolean = true): PathToFasta = {
    val ref  = makeTempFile("ref.", ".fa")
    Io.writeLines(ref, Fasta)

    if (withFai) {
      val fai  = ref.getParent.resolve(s"${ref.getFileName}.fai")
      fai.toFile.deleteOnExit()
      Io.writeLines(fai, Fai)
    }
    if (withDict) {
      val dict = PathUtil.replaceExtension(ref, ".dict")
      dict.toFile.deleteOnExit()
      Io.writeLines(dict, Dict)
    }
    ref
  }

  /** Sets up a trivial reference for testing against. */
  lazy val refFasta: PathToFasta = buildRefFasta()

  lazy val ref: ReferenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFasta)

  lazy val header = {
    val h = new SAMFileHeader
    h.setSequenceDictionary(ref.getSequenceDictionary)
    h
  }

  /** A VCF builder with the reference's contigs, the given samples and a FORMAT header for the allele fraction (AF). */
  private def vcfBuilder(samples: String*): VcfBuilder = {
    val contigs = VcfContigHeader.fromSequenceDictionary(header.getSequenceDictionary.fromSam)
    VcfBuilder(VcfBuilder.DefaultHeader.copy(contigs=contigs, samples=samples.toIndexedSeq))
      .withFormatHeaders(VcfFormatHeader("AF", VcfCount.OnePerAltAllele, VcfFieldType.Float, "Alt allele fraction."))
  }

  // We're going to simulate raw reads and consensuses as if there were variants at:
  //    chr1:10
  //    chr1:20
  //    chr1:30
  //    chr2:20
  lazy val (rawBam, consensusBam) : (PathToBam, PathToBam) = {
    val raw = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate))
    val con = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate), baseQuality=45)
    raw.header.setSequenceDictionary(ref.getSequenceDictionary)
    con.header.setSequenceDictionary(ref.getSequenceDictionary)

    // Add two raw read pairs and a consensus, with a mismatch at base 10 in the reference
    raw.addPair(name="A1", contig=0, start1=6, start2=50, attrs=Map("MI" -> "A"), bases1="AAAATAAAAA", bases2="AAAAAAAAAA")
    raw.addPair(name="A2", contig=0, start1=6, start2=50, attrs=Map("MI" -> "A"), bases1="AAAATAAAAG", bases2="AAAAAAAAAA")
    con.addPair(name="A" , contig=0, start1=6, start2=50, attrs=Map("MI" -> "A"), bases1="AAAATAAAAN", bases2="AAAAAAAAAA")


    // Two consensus, only one supporting the alt, with one raw read pair each, at chr1:20
    raw.addPair(name="B1", contig=0, start1=16, start2=50, attrs=Map("MI" -> "B"), bases1="AAAACAAAAA", bases2="AAAAAAAAAA")
    con.addPair(name="B" , contig=0, start1=16, start2=50, attrs=Map("MI" -> "B"), bases1="AAAACAAAAA", bases2="AAAAAAAAAA")

    raw.addPair(name="C1", contig=0, start1=17, start2=50, attrs=Map("MI" -> "C"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    con.addPair(name="C" , contig=0, start1=17, start2=50, attrs=Map("MI" -> "C"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")


    // A messy locus with consensuses with spanning deletions(D), no-calls(E), variant alleles(F) and mismatches at non-variant positions(G) at chr1:30
    raw.addPair(name="D1", contig=0, start1=25, start2=60, cigar1="4M4D6M", attrs=Map("MI" -> "D"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
    con.addPair(name="D" , contig=0, start1=25, start2=60, cigar1="4M4D6M", attrs=Map("MI" -> "D"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")

    raw.addPair(name="E1", contig=0, start1=26, start2=60, attrs=Map("MI" -> "E"), bases1="AAAANAAAAA", bases2="AAAAAAAAAA")
    con.addPair(name="E" , contig=0, start1=26, start2=60, attrs=Map("MI" -> "E"), bases1="AAAANAAAAA", bases2="AAAAAAAAAA")

    raw.addPair(name="F1", contig=0, start1=27, start2=60, attrs=Map("MI" -> "F"), bases1="AAAGAAAAAA", bases2="AAAAAAAAAA")
    con.addPair(name="F" , contig=0, start1=27, start2=60, attrs=Map("MI" -> "F"), bases1="AAAGAAAAAA", bases2="AAAAAAAAAA")

    raw.addPair(name="G1", contig=0, start1=28, start2=60, attrs=Map("MI" -> "G"), bases1="AAAAAAAGAA", bases2="AAAAAAAAAA")
    con.addPair(name="G" , contig=0, start1=28, start2=60, attrs=Map("MI" -> "G"), bases1="AAAAAAAGAA", bases2="AAAAAAAAAA")


    // A case where both ends of a pair overlap the variant position at chr2:20
    raw.addPair(name="H1", contig=1, start1=15, start2=19, attrs=Map("MI" -> "H"), bases1="CCCCCTCCCC", bases2="CTCCCCCCCC")
    con.addPair(name="H" , contig=1, start1=15, start2=19, attrs=Map("MI" -> "H"), bases1="CCCCCTCCCC", bases2="CTCCCCCCCC")


    // Some unmapped reads at the end of the file
    raw.addPair(name="X1", unmapped1 = true, unmapped2 = true, start1= -1, start2= -2, attrs=Map("MI" -> "X"))
    raw.addPair(name="X",  unmapped1 = true, unmapped2 = true, start1= -1, start2= -2, attrs=Map("MI" -> "X"))

    (raw.toTempFile(), con.toTempFile())
  }

  "ReviewConsensusVariants.toMi" should "retrieve the MI tag from reads" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamBuilder()
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "A")).foreach(toMi(_) shouldBe "A")
  }

  it should "strip the suffix from duplex read MIs" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamBuilder()
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/A")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/B")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/ReallyLongSuffix")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/456/A")).foreach(toMi(_) shouldBe "123/456")
  }

  it should "fail if a read without an MI tag is provided" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamBuilder()
    an[Exception] shouldBe thrownBy { builder.addPair(start1=1,start2=2,attrs=Map()).map(toMi) }
  }

  "ReviewConsensusVariants.toInsertString" should "return NA for anything that's not an FR pair" in {
    import ReviewConsensusVariants.toInsertString
    val builder = new SamBuilder(readLength=10)
    builder.addFrag(start=100).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addFrag(unmapped=true).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addPair(start1=100, start2=200, strand1=Plus, strand2=Plus).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addPair(start1=100, start2=100, unmapped2=true).foreach(r => toInsertString(r) shouldBe "NA")
  }

  it should "return generate sensible strings for FR pairs" in {
    import ReviewConsensusVariants.toInsertString
    val builder = new SamBuilder(readLength=10)
    builder.addPair(start1=100, start2=191).foreach(r => toInsertString(r) shouldBe "chr1:100-200 | F1R2")
    builder.addPair(start1=191, start2=100, strand1=Minus, strand2=Plus).foreach(r => toInsertString(r) shouldBe "chr1:100-200 | F2R1")
  }

  "ReviewConsensusVariants" should "create empty BAMs when given an empty interval list as input" in {
    val outBase = makeTempFile("review_consensus.", ".out")
    val conOut  = outBase.getParent.resolve(s"${outBase.getFileName}.consensus.bam")
    val rawOut  = outBase.getParent.resolve(s"${outBase.getFileName}.grouped.bam")
    val txtOut  = outBase.getParent.resolve(s"${outBase.getFileName}.txt")
    val intervals = makeTempFile("empty.", ".interval_list")
    new IntervalList(header).write(intervals)
    new ReviewConsensusVariants(input=intervals, consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=outBase).execute()

    conOut.toFile.exists() shouldBe true
    rawOut.toFile.exists() shouldBe true
    txtOut.toFile.exists() shouldBe true

    readBamRecs(rawOut) shouldBe empty
    readBamRecs(conOut) shouldBe empty
    Metric.read[ConsensusVariantReviewInfo](txtOut) shouldBe empty
  }

  it should "create empty BAMs when given an empty VCF as input" in {
    val outBase = makeTempFile("review_consensus.", ".out")
    val conOut  = outBase.getParent.resolve(s"${outBase.getFileName}.consensus.bam")
    val rawOut  = outBase.getParent.resolve(s"${outBase.getFileName}.grouped.bam")
    val txtOut  = outBase.getParent.resolve(s"${outBase.getFileName}.txt")
    val intervals = makeTempFile("empty.", ".vcf")

    vcfBuilder("s1").write(intervals)
    new ReviewConsensusVariants(input=intervals, consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=outBase).execute()

    conOut.toFile.exists() shouldBe true
    rawOut.toFile.exists() shouldBe true
    txtOut.toFile.exists() shouldBe true

    readBamRecs(rawOut) shouldBe empty
    readBamRecs(conOut) shouldBe empty
    Metric.read[ConsensusVariantReviewInfo](txtOut) shouldBe empty
  }

  it should "extract the right reads given a set of loci" in {
    val outBase = makeTempFile("review_consensus.", ".out")
    val conOut  = outBase.getParent.resolve(s"${outBase.getFileName}.consensus.bam")
    val rawOut  = outBase.getParent.resolve(s"${outBase.getFileName}.grouped.bam")
    val txtOut  = outBase.getParent.resolve(s"${outBase.getFileName}.txt")

    val vcf = vcfBuilder("tumor")
      .add(chrom="chr1", pos=10, alleles=Seq("A","T"), gts=Seq(Gt(sample="tumor", gt="A/T", attrs=Map("AF" -> Seq(0.01f)))))
      .add(chrom="chr1", pos=20, alleles=Seq("A","C"), gts=Seq(Gt(sample="tumor", gt="A/C", attrs=Map("AF" -> Seq(0.01f)))))
      .add(chrom="chr1", pos=30, alleles=Seq("A","G"), gts=Seq(Gt(sample="tumor", gt="A/G", attrs=Map("AF" -> Seq(0.01f)))))
      .add(chrom="chr2", pos=20, alleles=Seq("C","T"), gts=Seq(Gt(sample="tumor", gt="C/T", ad=Seq(100, 2))))
      .toTempFile()

    new ReviewConsensusVariants(input=vcf, consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=outBase).execute()

    conOut.toFile.exists() shouldBe true
    rawOut.toFile.exists() shouldBe true
    txtOut.toFile.exists() shouldBe true

    readBamRecs(rawOut).map(_.id) should contain theSameElementsAs Seq("A1/1", "A2/1", "B1/1", "D1/1", "E1/1", "F1/1", "H1/1", "H1/2")
    readBamRecs(conOut).map(_.id) should contain theSameElementsAs Seq("A/1", "B/1", "D/1", "E/1", "F/1", "H/1", "H/2")

    val metrics = Metric.read[ConsensusVariantReviewInfo](txtOut)
    // The metrics don't currently contain reads with spanning deletions, so D/1 is present above and absent below
    metrics.toIndexedSeq.map(_.consensus_read) should contain theSameElementsAs Seq("A/1", "B/1", "E/1", "F/1", "H/1", "H/2")
  }

  it should "fail if a FASTA sequence dictionary (.dict) is not present" in {
    val outBase   = makeTempFile("review_consensus.", ".out")
    val intervals = makeTempFile("empty.", ".vcf")
    val ref       = buildRefFasta(withDict=false)

    vcfBuilder("s1").write(intervals)
    a[ValidationException] should be thrownBy new ReviewConsensusVariants(input=intervals, consensusBam=consensusBam, groupedBam=rawBam, ref=ref, output=outBase)
  }

  it should "fail if a FASTA index (.fai) is not present" in {
    val outBase   = makeTempFile("review_consensus.", ".out")
    val intervals = makeTempFile("empty.", ".vcf")
    val ref       = buildRefFasta(withFai=false)

    vcfBuilder("s1").write(intervals)
    a[ValidationException] should be thrownBy new ReviewConsensusVariants(input=intervals, consensusBam=consensusBam, groupedBam=rawBam, ref=ref, output=outBase)
  }

  "ReviewConsensusVariants.loadPositions" should "load SNPs at or below the MAF with their genotype and filters" in {
    val vcf = vcfBuilder("s1")
      .add(chrom="chr1", pos=10, alleles=Seq("A","T"), filters=Seq("PASS"), gts=Seq(Gt(sample="s1", gt="A/T", attrs=Map("AF" -> Seq(0.05f)))))
      .add(chrom="chr1", pos=11, alleles=Seq("A","T"), gts=Seq(Gt(sample="s1", gt="A/T", attrs=Map("AF" -> Seq(0.06f)))))
      .add(chrom="chr1", pos=12, alleles=Seq("A","C","G"), filters=Seq("LowQD", "LowAB"), gts=Seq(Gt(sample="s1", gt="C|G", ad=Seq(98, 1, 1))))
      .add(chrom="chr1", pos=13, alleles=Seq("A","T"), gts=Seq(Gt(sample="s1", gt="A/T", ad=Seq(90, 10))))
      .add(chrom="chr1", pos=14, alleles=Seq("A","AT"), gts=Seq(Gt(sample="s1", gt="A/AT")))
      .add(chrom="chr1", pos=15, alleles=Seq("A","<NON_REF>"), gts=Seq(Gt(sample="s1", gt="A/A")))
      .add(chrom="chr1", pos=16, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="A/A")))
      .add(chrom="chr2", pos=20, alleles=Seq("c","t"), gts=Seq(Gt(sample="s1", gt="./.")))
      .toTempFile()

    val tool      = new ReviewConsensusVariants(input=vcf, consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=makeTempFile("review_consensus.", ".out"))
    val positions = tool.loadPositions(vcf, ref)
    positions.map(p => (p.chrom, p.start, p.refBase, p.genotype, p.filters)) should contain theSameElementsInOrderAs Seq(
      ("chr1", 10, 'A', Some("A/T"), None),
      ("chr1", 12, 'A', Some("C|G"), Some("LowAB,LowQD")),
      ("chr2", 20, 'C', Some("./."), None)
    )
  }

  it should "not filter on MAF or report genotypes when the VCF has many samples and no sample is given" in {
    val vcf = vcfBuilder("s1", "s2")
      .add(chrom="chr1", pos=10, alleles=Seq("A","T"), gts=Seq(Gt(sample="s1", gt="A/T", attrs=Map("AF" -> Seq(0.5f)))))
      .toTempFile()

    val tool = new ReviewConsensusVariants(input=vcf, consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=makeTempFile("review_consensus.", ".out"))
    tool.loadPositions(vcf, ref).map(p => (p.chrom, p.start, p.genotype)) should contain theSameElementsInOrderAs Seq(("chr1", 10, None))
  }
}
