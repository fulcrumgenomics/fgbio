/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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
package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.{Io, ReadStructure}

/**
  * Tests for AnnotateBamWithUmis
  */

class AnnotateBamWithUmisTest extends UnitSpec {
  private val dir     = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/umi")
  private val sam     = dir.resolve("annotate_umis.sam")
  private val fq      = dir.resolve("annotate_umis.fastq")
  private val umiTag  = "RX"
  private val qualTag = "QX"

  "AnnotateBamWithUmis" should "successfully add UMIs to a BAM" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag)
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
    })
  }

  it should "successfully add UMIs to a BAM in when the fastq is sorted" in {
    val out       = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, sorted=true)
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
    })
  }

  it should "successfully add UMI qualities to a BAM" in {
    val out       = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, qualAttribute=Some(qualTag))
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
      rec[String](qualTag) shouldBe rec.qualsString.substring(0,8)
    })
  }

  it should "successfully add UMI qualities to a BAM when the fastq is sorted" in {
    val out       = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, qualAttribute=Some(qualTag), sorted=true)
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
      rec[String](qualTag) shouldBe rec.qualsString.substring(0,8)
    })
  }

  it should "fail if one or more reads doesn't have a UMI" in {
    val out     = makeTempFile("with_umis.", ".bam")
    val shortFq = makeTempFile("missing_umis.", ".fq.gz")
    Io.writeLines(shortFq, Io.readLines(fq).toSeq.dropRight(8))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(shortFq), output=out, attribute=umiTag)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "fail if one or more reads doesn't have a UMI when the fastq is sorted" in {
    val out     = makeTempFile("with_umis.", ".bam")
    val shortFq = makeTempFile(s"missing_umis.", ".fq.gz")
    Io.writeLines(shortFq, Io.readLines(fq).toSeq.dropRight(8))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(shortFq), output=out, attribute=umiTag, sorted=true)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "not fail if there are extra reads in the fastq not in the bam" in {
    val out    = makeTempFile("with_umis.", ".bam")
    val longFq = makeTempFile(s"extra_umis.", ".fq.gz")
    val lines  = Io.readLines(fq) ++ Seq("@not_a_flowcell:1:1101:10060:3200/2 2:N:0:19","GATCTTGG","+","-,86,,;:")
    Io.writeLines(longFq, lines)
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(longFq), output=out, attribute=umiTag)
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
    })
  }

  it should "not fail if there are extra reads in the fastq not in the bam when the fastq is sorted" in {
    val out    = makeTempFile("with_umis.", ".bam")
    val longFq = makeTempFile(s"extra_umis.", ".fq.gz")
    val lines  = Io.readLines(fq) ++ Seq("@not_a_flowcell:1:1101:10060:3200/2 2:N:0:19","GATCTTGG","+","-,86,,;:")
    Io.writeLines(longFq, lines)
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(longFq), output=out, attribute=umiTag, sorted=true)
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
    })
  }

  it should "successfully add UMIs to a BAM with a given read structure" in {
    val out       = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, readStructure=Seq(ReadStructure("2B4M+B")))
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(2, 6)
    })
  }

  it should "successfully add UMIs to a BAM with a given read structure with multiple molecular barcodes" in {
    val out       = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, qualAttribute=Some(qualTag), readStructure=Seq(ReadStructure("2M2B+M")))
    annotator.execute()
    SamSource(out).foreach(rec => {
      val bases = rec.basesString
      val quals = rec.qualsString
      rec[String](umiTag) shouldBe bases.substring(0,2) + "-" + bases.substring(4, 8)
      rec[String](qualTag) shouldBe quals.substring(0,2) + " " + quals.substring(4, 8)
    })
  }

  it should "fail to add UMIs to a BAM wiht one unsorted FASTQ and corresponding read structures when sorted is true" in {
    val out = makeTempFile("with_umis.", ".bam")
    val reverseFq = makeTempFile(s"reverse_umis.", ".fq.gz")
    Io.writeLines(reverseFq, Io.readLines(fq).grouped(4).toSeq.reverse.flatten)
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(reverseFq), output=out, attribute=umiTag, readStructure=Seq(ReadStructure("2M4B+M")), sorted=true)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "successfully add UMIs to a BAM using all bases from multiple FASTQs" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, fq), readStructure=Seq(ReadStructure("2B+M")), output=out, attribute=umiTag)
    annotator.execute()
    SamSource(out).foreach { rec =>
      val bases = rec.basesString.substring(2,8)
      rec[String](umiTag) shouldBe s"${bases}-${bases}"
    }
  }

  it should "successfully add UMIs to a BAM with multiple FASTQs and corresponding read structures" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, fq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag)
    annotator.execute()
    SamSource(out).foreach { rec =>
      val first  = rec.basesString.substring(0,2) // read one: [2M]4B+M
      val second = rec.basesString.substring(6,8) // read one: 2M4B[+M]
      val third  = rec.basesString.substring(1,8) // read two: 1B[+M]
      rec[String](umiTag) shouldBe s"${first}-${second}-${third}"
    }
  }

  it should "successfully add UMIs to a BAM with multiple FASTQs and corresponding read structures when sorted is true" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, fq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, sorted=true)
    annotator.execute()
    SamSource(out).foreach { rec =>
      val first  = rec.basesString.substring(0,2) // read one: [2M]4B+M
      val second = rec.basesString.substring(6,8) // read one: 2M4B[+M]
      val third  = rec.basesString.substring(1,8) // read two: 1B[+M]
      rec[String](umiTag) shouldBe s"${first}-${second}-${third}"
    }
  }

  it should "fail to add UMIs to a BAM with multiple truncated FASTQs and corresponding read structures when sorted is true" in {
    val out = makeTempFile("with_umis.", ".bam")
    val shortFq = makeTempFile(s"missing_umis.", ".fq.gz")
    Io.writeLines(shortFq, Io.readLines(fq).toSeq.dropRight(8))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(shortFq, shortFq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, sorted=true)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "fail to add UMIs to a BAM with multiple unsorted FASTQs and corresponding read structures when sorted is true" in {
    val out = makeTempFile("with_umis.", ".bam")
    val reverseFq = makeTempFile(s"reverse_umis.", ".fq.gz")
    Io.writeLines(reverseFq, Io.readLines(fq).grouped(4).toSeq.reverse.flatten)
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(reverseFq, reverseFq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, sorted=true)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "fail to add UMIs to a BAM with multiple FASTQs with extra lines and corresponding read structures when sorted is true" in {
    val out = makeTempFile("with_umis.", ".bam")
    val longFq = makeTempFile(s"extra_umis.", ".fq.gz")
    Io.writeLines(longFq, Io.readLines(fq))
    Io.writeLines(longFq, Seq("@not_a_flowcell:1:1101:10060:3200/2 2:N:0:19","GATCTTGG","+","-,86,,;:"))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(longFq, longFq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, sorted=true)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "fail if FASTQs are not sorted the same" in {
    val out = makeTempFile("with_umis.", ".bam")
    val reverseFq = makeTempFile(s"reverse_umis.", ".fq.gz")
    Io.writeLines(reverseFq, Io.readLines(fq).grouped(4).toSeq.reverse.flatten)
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, reverseFq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, sorted=true)
    val result = intercept[Exception] { annotator.execute() }
    println(result)
    result.getMessage should include ("out of sync")
  }
}
