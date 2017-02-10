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
import com.fulcrumgenomics.testing.SamRecordSetBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec, VariantContextSetBuilder}
import com.fulcrumgenomics.umi.ReviewConsensusVariants.ConsensusVariantReviewInfo
import com.fulcrumgenomics.util.{Io, Metric}
import dagr.commons.io.PathUtil
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.samtools.{SAMFileHeader, SAMRecord}

object ReviewConsensusVariantsTest {
  val Fasta : Seq[String] =
    """
      |>chr1
      |AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      |AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      |>chr2
      |CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      |CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    """.stripMargin.lines.toSeq

  val Fai : Seq[String] =
    """
      |chr1    100     6       50      51
      |chr2    100     114     50      51
    """.stripMargin.lines.map(_.trim.replaceAll("\\s+", "\t")).filter(_.nonEmpty).toSeq

  val Dict : Seq[String] =
    """
      |@HD    	VN:1.5 	SO:unsorted
      |@SQ    	SN:chr1	LN:100 	M5:8adc5937e635f6c9af646f0b23560fae    	AS:n/a 	UR:n/a 	SP:n/a
      |@SQ    	SN:chr2	LN:100 	M5:48ecc9d349522f836cadf615e370bc51    	AS:n/a 	UR:n/a 	SP:n/a
    """.stripMargin.lines.map(_.trim.replaceAll("\\s+", "\t")).filter(_.nonEmpty).toSeq
}


class ReviewConsensusVariantsTest extends UnitSpec {
  import ReviewConsensusVariantsTest._

  /** Sets up a trivial reference for testing against. */
  lazy val refFasta: PathToFasta = {
    val ref  = makeTempFile("ref.", ".fa")
    val fai  = ref.getParent.resolve(ref.getFileName + ".fai")
    val dict = PathUtil.replaceExtension(ref, ".dict")
    Seq(fai, dict).foreach(_.toFile.deleteOnExit())

    Io.writeLines(ref, Fasta)
    Io.writeLines(fai, Fai)
    Io.writeLines(dict, Dict)
    ref
  }

  lazy val ref: ReferenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFasta)

  lazy val header = {
    val h = new SAMFileHeader
    h.setSequenceDictionary(ref.getSequenceDictionary)
    h
  }

  // We're going to simulate raw reads and consensuses as if there were variants at:
  //    chr1:10
  //    chr1:20
  //    chr1:30
  //    chr2:20
  lazy val (rawBam, consensusBam) : (PathToBam, PathToBam) = {
    val raw = new SamRecordSetBuilder(readLength=10, sortOrder=SortOrder.coordinate)
    val con = new SamRecordSetBuilder(readLength=10, sortOrder=SortOrder.coordinate, baseQuality=45)
    raw.header.setSequenceDictionary(ref.getSequenceDictionary)
    con.header.setSequenceDictionary(ref.getSequenceDictionary)

    // Add two raw read pairs and a consensus, with a mismatch at base 10 in the reference
    raw.addPair(name="A1", contig=0, start1=6, start2=50, attrs=Map("MI" -> "A")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAATAAAAA" else "AAAAAAAAAA"))
    raw.addPair(name="A2", contig=0, start1=6, start2=50, attrs=Map("MI" -> "A")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAATAAAAG" else "AAAAAAAAAA"))
    con.addPair(name="A" , contig=0, start1=6, start2=50, attrs=Map("MI" -> "A")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAATAAAAN" else "AAAAAAAAAA"))


    // Two consensus, only one supporting the alt, with one raw read pair each, at chr1:20
    raw.addPair(name="B1", contig=0, start1=16, start2=50, attrs=Map("MI" -> "B")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAACAAAAA" else "AAAAAAAAAA"))
    con.addPair(name="B" , contig=0, start1=16, start2=50, attrs=Map("MI" -> "B")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAACAAAAA" else "AAAAAAAAAA"))

    raw.addPair(name="C1", contig=0, start1=17, start2=50, attrs=Map("MI" -> "C")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "AAAAAAAAAA"))
    con.addPair(name="C" , contig=0, start1=17, start2=50, attrs=Map("MI" -> "C")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "AAAAAAAAAA"))


    // A messy locus with consensuses with spanning deletions(D), no-calls(E), variant alleles(F) and mismatches at non-variant positions(G) at chr1:30
    raw.addPair(name="D1", contig=0, start1=25, start2=60, cigar1="4M4D6M", attrs=Map("MI" -> "D")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "AAAAAAAAAA"))
    con.addPair(name="D" , contig=0, start1=25, start2=60, cigar1="4M4D6M", attrs=Map("MI" -> "D")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "AAAAAAAAAA"))

    raw.addPair(name="E1", contig=0, start1=26, start2=60, attrs=Map("MI" -> "E")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAANAAAAA" else "AAAAAAAAAA"))
    con.addPair(name="E" , contig=0, start1=26, start2=60, attrs=Map("MI" -> "E")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAANAAAAA" else "AAAAAAAAAA"))

    raw.addPair(name="F1", contig=0, start1=27, start2=60, attrs=Map("MI" -> "F")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAGAAAAAA" else "AAAAAAAAAA"))
    con.addPair(name="F" , contig=0, start1=27, start2=60, attrs=Map("MI" -> "F")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAGAAAAAA" else "AAAAAAAAAA"))

    raw.addPair(name="G1", contig=0, start1=28, start2=60, attrs=Map("MI" -> "G")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAGAA" else "AAAAAAAAAA"))
    con.addPair(name="G" , contig=0, start1=28, start2=60, attrs=Map("MI" -> "G")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAGAA" else "AAAAAAAAAA"))


    // A case where both ends of a pair overlap the variant position at chr2:20
    raw.addPair(name="H1", contig=1, start1=15, start2=19, attrs=Map("MI" -> "H")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "CCCCCTCCCC" else "CTCCCCCCCC"))
    con.addPair(name="H" , contig=1, start1=15, start2=19, attrs=Map("MI" -> "H")).foreach(r => r.setReadString(if (r.getFirstOfPairFlag) "CCCCCTCCCC" else "CTCCCCCCCC"))


    // Some unmapped reads at the end of the file
    raw.addPair(name="X1", record1Unmapped = true, record2Unmapped = true, start1= -1, start2= -2, attrs=Map("MI" -> "X"))
    raw.addPair(name="X",  record1Unmapped = true, record2Unmapped = true, start1= -1, start2= -2, attrs=Map("MI" -> "X"))

    (raw.toTempFile(), con.toTempFile())
  }

  /** Returns queryname + /1 or /2 */
  private def nameAndNumber(rec: SAMRecord): String = rec.getReadName + "/" + (if (rec.getFirstOfPairFlag) "1" else "2")

  "ReviewConsensusVariants.toMi" should "retrieve the MI tag from reads" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamRecordSetBuilder()
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "A")).foreach(toMi(_) shouldBe "A")
  }

  it should "strip the suffix from duplex read MIs" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamRecordSetBuilder()
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/A")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/B")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/ReallyLongSuffix")).foreach(toMi(_) shouldBe "123")
    builder.addPair(start1=1,start2=2,attrs=Map("MI" -> "123/456/A")).foreach(toMi(_) shouldBe "123/456")
  }

  it should "fail if a read without an MI tag is provided" in {
    import ReviewConsensusVariants.toMi
    val builder = new SamRecordSetBuilder()
    an[Exception] shouldBe thrownBy { builder.addPair(start1=1,start2=2,attrs=Map()).map(toMi) }
  }

  "ReviewConsensusVariants.toInsertString" should "return NA for anything that's not an FR pair" in {
    import ReviewConsensusVariants.toInsertString
    val builder = new SamRecordSetBuilder(readLength=10)
    builder.addFrag(start=100).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addFrag(unmapped=true).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addPair(start1=100, start2=200, strand1=Plus, strand2=Plus).foreach(r => toInsertString(r) shouldBe "NA")
    builder.addPair(start1=100, start2=100, record2Unmapped=true).foreach(r => toInsertString(r) shouldBe "NA")
  }

  it should "return generate sensible strings for FR pairs" in {
    import ReviewConsensusVariants.toInsertString
    val builder = new SamRecordSetBuilder(readLength=10)
    builder.addPair(start1=100, start2=191).foreach(r => toInsertString(r) shouldBe "chr1:100-200 | F1R2")
    builder.addPair(start1=191, start2=100, strand1=Minus, strand2=Plus).foreach(r => toInsertString(r) shouldBe "chr1:100-200 | F2R1")
  }

  "ReviewConsensusVariants" should "create empty BAMs when given an empty interval list as input" in {
    val outBase = makeTempFile("review_consensus.", ".out")
    val conOut  = outBase.getParent.resolve(outBase.getFileName + ".consensus.bam")
    val rawOut  = outBase.getParent.resolve(outBase.getFileName + ".grouped.bam")
    val txtOut  = outBase.getParent.resolve(outBase.getFileName + ".txt")
    val intervals = makeTempFile("empty.", ".interval_list")
    new IntervalList(header).write(intervals.toFile)
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
    val conOut  = outBase.getParent.resolve(outBase.getFileName + ".consensus.bam")
    val rawOut  = outBase.getParent.resolve(outBase.getFileName + ".grouped.bam")
    val txtOut  = outBase.getParent.resolve(outBase.getFileName + ".txt")
    val intervals = makeTempFile("empty.", ".vcf")

    val vcfBuilder = VariantContextSetBuilder("s1")
    vcfBuilder.header.setSequenceDictionary(header.getSequenceDictionary)
    vcfBuilder.write(intervals)
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
    val conOut  = outBase.getParent.resolve(outBase.getFileName + ".consensus.bam")
    val rawOut  = outBase.getParent.resolve(outBase.getFileName + ".grouped.bam")
    val txtOut  = outBase.getParent.resolve(outBase.getFileName + ".txt")

    val vcfBuilder = new VariantContextSetBuilder(sampleNames=List("tumor"))
    vcfBuilder.header.setSequenceDictionary(this.header.getSequenceDictionary)
    vcfBuilder.addVariant(refIdx=0, start=10, variantAlleles=List("A","T"), genotypeAlleles=List("A","T"), genotypeAttributes=Map("AF" -> 0.01))
    vcfBuilder.addVariant(refIdx=0, start=20, variantAlleles=List("A","C"), genotypeAlleles=List("A","C"), genotypeAttributes=Map("AF" -> 0.01))
    vcfBuilder.addVariant(refIdx=0, start=30, variantAlleles=List("A","G"), genotypeAlleles=List("A","G"), genotypeAttributes=Map("AF" -> 0.01))
    vcfBuilder.addVariant(refIdx=1, start=20, variantAlleles=List("C","T"), genotypeAlleles=List("C","T"), genotypeAttributes=Map("AD" -> Array(100,2)))

    new ReviewConsensusVariants(input=vcfBuilder.toTempFile(), consensusBam=consensusBam, groupedBam=rawBam, ref=refFasta, output=outBase).execute()

    conOut.toFile.exists() shouldBe true
    rawOut.toFile.exists() shouldBe true
    txtOut.toFile.exists() shouldBe true

    readBamRecs(rawOut).map(nameAndNumber) should contain theSameElementsAs Seq("A1/1", "A2/1", "B1/1", "D1/1", "E1/1", "F1/1", "H1/1", "H1/2")
    readBamRecs(conOut).map(nameAndNumber) should contain theSameElementsAs Seq("A/1", "B/1", "D/1", "E/1", "F/1", "H/1", "H/2")

    val metrics = Metric.read[ConsensusVariantReviewInfo](txtOut)
    // The metrics don't currently contain reads with spanning deletions, so D/1 is present above and absent below
    metrics.toIndexedSeq.map(_.consensus_read) should contain theSameElementsAs Seq("A/1", "B/1", "E/1", "F/1", "H/1", "H/2")
  }
}
