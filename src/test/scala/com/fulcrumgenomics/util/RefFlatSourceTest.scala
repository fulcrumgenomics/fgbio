/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import java.nio.file.Paths
import java.util.Collections

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.GeneAnnotations.Exon
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}

class RefFlatSourceTest extends UnitSpec {

  // Adapted from http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refFlat.txt.gz
  private val RefFlatFile = Paths.get("src/test/resources/com/fulcrumgenomics/util/refFlat.txt.gz")

  private val Dictionary = {
    new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("chr1", 249250621)))
  }

  "RefFlatSource" should "read valid refFlat from various kinds of input resources" in {
    val rf = RefFlatFile
    val lines = Io.readLines(rf)

    Seq(
      RefFlatSource(rf, Some(Dictionary)),
      RefFlatSource(rf.toFile, Some(Dictionary)),
      RefFlatSource(Io.toInputStream(rf), Some(Dictionary)),
      RefFlatSource(lines, Some(Dictionary))
    ).foreach(source => {
      val genes = source.toSeq

      // # of genes
      genes.length shouldBe 142

      // # of transcripts across all genes
      genes.map(_.size).sum shouldBe 360

      // # of exons across all transcripts across all genes
      genes.map { gene => gene.map(_.exons.size).sum }.sum shouldBe 4613

      // verify the first gene
      val gene = genes.head
      gene.contig shouldBe "chr1"
      gene.start shouldBe 141638944
      gene.end shouldBe 141655128
      gene.negativeStrand shouldBe true
      gene.name shouldBe "ANKRD20A12P"
      gene.size shouldBe 1

      // verify the first transcript in the first gene
      val transcript = gene.head
      transcript.name shouldBe "NR_046228"
      transcript.start shouldBe 141638944
      transcript.end shouldBe 141655128
      transcript.cdsStart shouldBe 141655129
      transcript.cdsEnd shouldBe 141655128
      transcript.exons.size shouldBe 5

      //verify the first exon in the first transcript
      val exon = transcript.exons.head
      exon shouldBe Exon(141638944, 141639396)
    })
  }

  it should "filter genes where chromosomes differ across transcripts" in {
    val lines = Iterator(
      Seq("ACKR4", "NM_178445-1", "chr1", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t"),
      Seq("ACKR4", "NM_178445-2", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t")
    )
    val source = RefFlatSource(lines, dict=None).toSeq
    source should have size 1
    source.head should have size 1
  }

  it should "filter genes where strands differ across transcripts" in {
    val lines = Iterator(
      Seq("ACKR4", "NM_178445-1", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t"),
      Seq("ACKR4", "NM_178445-2", "chr3", "-", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t")
    )
    val source = RefFlatSource(lines, dict=None).toSeq
    source should have size 1
    source.head should have size 1
  }

  it should "fail if the # of exon starts or ends do not equal the exon count for a transcript" in {
    val startsMismatch = Iterator(
      Seq("ACKR4", "NM_178445", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670,133801680", "133804175").mkString("\t")
    )
    val endsMismatch = Iterator(
      Seq("ACKR4", "NM_178445", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175,133801685").mkString("\t")
    )
    val bothMismatch = Iterator(
      Seq("ACKR4", "NM_178445", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670,133801680", "133804175,133804185").mkString("\t")
    )
    an[Exception] should be thrownBy RefFlatSource(startsMismatch, dict=None)
    an[Exception] should be thrownBy RefFlatSource(endsMismatch, dict=None)
    an[Exception] should be thrownBy RefFlatSource(bothMismatch, dict=None)
  }

  it should "filter out genes with chromosomes not in the sequence dictionary" in {
    val lines = Iterator(
      Seq("ACKR4", "NM_178445-1", "chr1", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t"),
      Seq("ACKR4", "NM_178445-1", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175").mkString("\t")
    )
    RefFlatSource(lines, dict=Some(Dictionary)) should have size 1
  }
}
