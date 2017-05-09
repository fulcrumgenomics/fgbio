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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.commons.io.PathUtil
import htsjdk.samtools.SamReaderFactory

/**
  * Tests for AnnotateBamWithUmis
  */
class AnnotateBamWithUmisTest extends UnitSpec {
  val dir = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/umi")
  val sam = dir.resolve("annotate_umis.sam")
  val fq  = dir.resolve("annotate_umis.fastq")
  val umiTag    = "RX"
  val umiLength = 8

  "AnnotateBamWithUmis" should "successfully add UMIs to a BAM in" in {
    val out = makeTempFile("with_umis.", ".bam")
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=fq, output=out, attribute=umiTag)
    annotator.execute()
    SamReaderFactory.make().open(out.toFile).iterator().foreach(rec => {
      rec.getStringAttribute(umiTag) shouldBe rec.getReadString.take(8)
    })
  }

  "AnnotateBamWithUmis" should "fail if one or more reads doesn't have a UMI" in {
    val out     = makeTempFile("with_umis.", ".bam")
    val shortFq = makeTempFile("missing_umis.", ".fq.gz")
    Io.writeLines(shortFq, Io.readLines(fq).toSeq.dropRight(8))
    val annotator = new AnnotateBamWithUmis(input=sam, fastq=shortFq, output=out, attribute=umiTag)
    an[FailureException] shouldBe thrownBy { annotator.execute() }
  }
}
