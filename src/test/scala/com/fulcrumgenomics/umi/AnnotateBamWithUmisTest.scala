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
  private val dir = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/umi")
  private val sam = dir.resolve("annotate_umis.sam")
  private val fq  = dir.resolve("annotate_umis.fastq")
  private val umiTag    = "RX"

  "AnnotateBamWithUmis" should "successfully add UMIs to a BAM" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag)
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(0,8)
    })
  }

  it should "fail if one or more reads doesn't have a UMI" in {
    val out     = makeTempFile("with_umis.", ".bam")
    val shortFq = makeTempFile("missing_umis.", ".fq.gz")
    Io.writeLines(shortFq, Io.readLines(fq).toSeq.dropRight(8))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(shortFq), output=out, attribute=umiTag)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }

  it should "successfully add UMIs to a BAM with a given read structure in" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq), output=out, attribute=umiTag, readStructure=Seq(ReadStructure("2B4M+B")))
    annotator.execute()
    SamSource(out).foreach(rec => {
      rec[String](umiTag) shouldBe rec.basesString.substring(2, 6)
    })
  }

  it should "successfully add UMIs to a BAM using all bases from multiple FASTQs" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, fq), readStructure=Seq(ReadStructure("2B+M")), output=out, attribute=umiTag, delimiter="x")
    annotator.execute()
    SamSource(out).foreach { rec =>
      val bases = rec.basesString.substring(2,8)
      rec[String](umiTag) shouldBe s"${bases}x${bases}"
    }
  }

  it should "successfully add UMIs to a BAM with multiple FASTQs and corresponding read structures" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=Seq(fq, fq), readStructure=Seq(ReadStructure("2M4B+M"), ReadStructure("1B+M")), output=out, attribute=umiTag, delimiter="x")
    annotator.execute()
    SamSource(out).foreach { rec =>
      val first  = rec.basesString.substring(0,2) // read one: [2M]4B+M
      val second = rec.basesString.substring(6,8) // read one: 2M4B[+M]
      val third  = rec.basesString.substring(1,8) // read two: 1B[+M]
      rec[String](umiTag) shouldBe s"${first}x${second}x${third}"
    }
  }
}
