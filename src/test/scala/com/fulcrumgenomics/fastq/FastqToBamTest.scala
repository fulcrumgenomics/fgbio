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

package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.ReadStructure
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.util.Iso8601Date

class FastqToBamTest extends UnitSpec {
  /** Writes one or more fastq records to a file. */
  def fq(recs: FastqRecord*): PathToFastq = {
    val path = makeTempFile("fastqToBamTest.", ".fq")
    val out = FastqWriter(path)
    recs.foreach(out.write)
    out.close()
    path
  }

  def rs(rs: String): ReadStructure = ReadStructure(rs)

  "FastqToBam" should "make a BAM from a single fastq file without a read structure" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="), FastqRecord("q2", "CCCCCCCCCC", "##########"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1), output=bam, sample="foo", library="bar").execute()
    val recs = readBamRecs(bam)

    recs should have size 2
    val Seq(q1, q2) = recs
    q1.getReadName shouldBe "q1"
    q1.getReadPairedFlag shouldBe false
    q1.getReadString shouldBe "AAAAAAAAAA"
    q1.getBaseQualities.forall(_ == 28) shouldBe true
    q1.getReadGroup.getSample shouldBe "foo"
    q1.getReadGroup.getLibrary shouldBe "bar"

    q2.getReadName shouldBe "q2"
    q2.getReadPairedFlag shouldBe false
    q2.getReadString shouldBe "CCCCCCCCCC"
    q2.getBaseQualities.forall(_ == 2) shouldBe true
    q2.getReadGroup.getSample shouldBe "foo"
    q2.getReadGroup.getLibrary shouldBe "bar"
  }

  it should "make a BAM from two fastq files without any read structures" in {
    val r1s = fq(FastqRecord("q1", "AAAAAAAAAA", "==========", readNumber=Some(1)))
    val r2s = fq(FastqRecord("q1", "CCCCCCCCCC", "##########", readNumber=Some(2)))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1s, r2s), output=bam, sample="pip", library="pop").execute()
    val recs = readBamRecs(bam)

    recs should have size 2
    val Seq(r1, r2) = recs
    r1.getReadName shouldBe "q1"
    r1.getReadPairedFlag shouldBe true
    r1.getFirstOfPairFlag shouldBe true
    r1.getSecondOfPairFlag shouldBe false
    r1.getReadString shouldBe "AAAAAAAAAA"
    r1.getBaseQualities.forall(_ == 28) shouldBe true
    r1.getReadGroup.getSample shouldBe "pip"
    r1.getReadGroup.getLibrary shouldBe "pop"

    r2.getReadName shouldBe "q1"
    r2.getReadPairedFlag shouldBe true
    r2.getFirstOfPairFlag shouldBe false
    r2.getSecondOfPairFlag shouldBe true
    r2.getReadString shouldBe "CCCCCCCCCC"
    r2.getBaseQualities.forall(_ == 2) shouldBe true
    r2.getReadGroup.getSample shouldBe "pip"
    r2.getReadGroup.getLibrary shouldBe "pop"
  }

  it should "make a BAM file with two fastqs with an inline UMI in R1" in {
    val r1s = fq(FastqRecord("q1", "ACGTAAAAAA", "==========", readNumber=Some(1)))
    val r2s = fq(FastqRecord("q1", "CCCCCCCCCC", "##########", readNumber=Some(2)))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1s, r2s), readStructures=Seq(rs("4M+T"), rs("+T")), output=bam, sample="s", library="l").execute()
    val recs = readBamRecs(bam)

    recs should have size 2
    val Seq(r1, r2) = recs
    r1.getReadName shouldBe "q1"
    r1.getReadString shouldBe "AAAAAA"
    r1.getBaseQualityString shouldBe "======"
    r1.getAttribute(ConsensusTags.UmiBases) shouldBe "ACGT"

    r2.getReadName shouldBe "q1"
    r2.getReadString shouldBe "CCCCCCCCCC"
    r2.getBaseQualityString shouldBe "##########"
    r2.getAttribute(ConsensusTags.UmiBases) shouldBe "ACGT"
  }

  it should "correctly handle complicated read structures with multiple UMI and sample barcode segments" in {
    val r1s = fq(FastqRecord("q1", "AAACCCTTTAAAAA", "==============", readNumber=Some(1)))
    val r2s = fq(FastqRecord("q1", "GGGTTTAAACCCCC", "##############", readNumber=Some(2)))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1s, r2s), readStructures=Seq(rs("3B3M3B5T"), rs("3B3M3B5T")), output=bam, sample="s", library="l").execute()
    val recs = readBamRecs(bam)

    recs should have size 2
    val Seq(r1, r2) = recs
    r1.getReadName shouldBe "q1"
    r1.getReadString shouldBe "AAAAA"
    r1.getBaseQualityString shouldBe "====="
    r1.getAttribute(ConsensusTags.UmiBases) shouldBe "CCC-TTT"
    r1.getAttribute("BC") shouldBe "AAA-TTT-GGG-AAA"

    r2.getReadName shouldBe "q1"
    r2.getReadString shouldBe "CCCCC"
    r2.getBaseQualityString shouldBe "#####"
    r2.getAttribute(ConsensusTags.UmiBases) shouldBe "CCC-TTT"
    r2.getAttribute("BC") shouldBe "AAA-TTT-GGG-AAA"
  }

  it should "use four fastqs to make reads" in {
    val r1s = fq(FastqRecord("q1", "AAAAAAAAAA", "==========", readNumber=Some(1)))
    val r2s = fq(FastqRecord("q1", "CCCCCCCCCC", "##########", readNumber=Some(2)))
    val r3s = fq(FastqRecord("q1", "ACGT", "????"))
    val r4s = fq(FastqRecord("q1", "GAGA", "????"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1s, r2s, r3s, r4s), readStructures=Seq(rs("+T"), rs("+T"), rs("4B"), rs("4M")), output=bam, sample="s", library="l").execute()
    val recs = readBamRecs(bam)

    recs should have size 2
    val Seq(r1, r2) = recs
    r1.getReadName shouldBe "q1"
    r1.getReadString shouldBe "AAAAAAAAAA"
    r1.getBaseQualityString shouldBe "=========="
    r1.getAttribute(ConsensusTags.UmiBases) shouldBe "GAGA"
    r1.getAttribute("BC") shouldBe "ACGT"

    r2.getReadName shouldBe "q1"
    r2.getReadString shouldBe "CCCCCCCCCC"
    r2.getBaseQualityString shouldBe "##########"
    r2.getAttribute(ConsensusTags.UmiBases) shouldBe "GAGA"
    r2.getAttribute("BC") shouldBe "ACGT"
  }

  it should "put everything in the right place in the header" in {
    val r1  = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(
      input=Seq(r1),
      readStructures=Seq(rs("10T")),
      output=bam,
      sample="foo",
      library="bar",
      readGroupId="MyRG",
      platform="Illumina",
      platformUnit=Some("pee-eww"),
      platformModel=Some("hiseq2500"),
      sequencingCenter=Some("nowhere"),
      predictedInsertSize=Some(300),
      description=Some("Some reads!"),
      runDate=Some(new Iso8601Date("")),
      comment=List("hello world", "comment two")).execute()

    val reader = SamReaderFactory.make().open(bam)
    val header = reader.getFileHeader
    header.getComments.iterator().toSeq shouldBe Seq("@CO\thello world", "@CO\tcomment two") // ewww
    header.getReadGroups.size() shouldBe 1
    header.getSortOrder shouldBe SortOrder.unsorted
    header.getGroupOrder shouldBe GroupOrder.query

    val rg = header.getReadGroups.iterator().next()
    rg.getSample shouldBe "foo"
    rg.getLibrary shouldBe "bar"
    rg.getReadGroupId shouldBe "MyRG"
    rg.getPlatform shouldBe "Illumina"
    rg.getPlatformUnit shouldBe "pee-eww"
    rg.getPlatformModel shouldBe "hiseq2500"
    rg.getSequencingCenter shouldBe "nowhere"
    rg.getPredictedMedianInsertSize shouldBe 300
    rg.getDescription shouldBe "Some reads!"

    reader.iterator().next().getReadGroup shouldBe rg
  }

  it should "sort the output bam if asked to" in {
    val r1 = fq(FastqRecord("q2", "AAAAAAAAAA", "=========="), FastqRecord("q10", "CCCCCCCCCC", "##########"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1), output=bam, sample="s", library="l", sort=true).execute()
    readBamRecs(bam).map(_.getReadName) shouldBe Seq("q10", "q2")
  }

  it should "work with fastq files with zero length reads" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="), FastqRecord("q2", "", ""))
    val r2 = fq(FastqRecord("q1", "", ""),                     FastqRecord("q2", "CCCCCCCCCC", "##########"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    new FastqToBam(input=Seq(r1, r2), output=bam, sample="s", library="l").execute()

    val recs = readBamRecs(bam)
    recs should have size 4
    recs(0).getReadString shouldBe "AAAAAAAAAA"
    recs(1).getReadString shouldBe "N"
    recs(2).getReadString shouldBe "N"
    recs(3).getReadString shouldBe "CCCCCCCCCC"
  }

  it should "fail when read names don't match up" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="))
    val r2 = fq(FastqRecord("x1", "CCCCCCCCCC", "##########"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    an[Exception] shouldBe thrownBy { new FastqToBam(input=Seq(r1, r2), output=bam, sample="s", library="l").execute() }
  }

  it should "fail when one fastq has more reads than another" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="), FastqRecord("q2", "TTTTTTTTTT", "??????????"))
    val r2 = fq(FastqRecord("q1", "CCCCCCCCCC", "##########"))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    an[Exception] shouldBe thrownBy { new FastqToBam(input=Seq(r1, r2), output=bam, sample="s", library="l").execute() }
  }

  it should "fail when the number of fastqs and read structures are mismatched" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    an[Exception] shouldBe thrownBy { new FastqToBam(input=Seq(r1), readStructures=Seq(rs("+T"), rs("+T")), output=bam, sample="s", library="l").execute() }
  }

  it should "fail when a read is not long enough for the read structure" in {
    val r1 = fq(FastqRecord("q1", "AAAAAAAAAA", "=========="))
    val bam = makeTempFile("fastqToBamTest.", ".bam")
    an[Exception] shouldBe thrownBy { new FastqToBam(input=Seq(r1), readStructures=Seq(rs("8M8S+T")), output=bam, sample="s", library="l").execute() }
  }
}
