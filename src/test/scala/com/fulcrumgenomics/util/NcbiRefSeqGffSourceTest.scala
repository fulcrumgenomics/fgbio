/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import java.nio.file.Paths
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.GeneAnnotations.Exon
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}
import org.scalatest.OptionValues

class NcbiRefSeqGffSourceTest extends UnitSpec with OptionValues {
  // Excerpted from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
  private val GffFile = Paths.get("src/test/resources/com/fulcrumgenomics/util/human.gff.gz")

  private val Chr1 = new SAMSequenceRecord("chr1", 249250621)
  Chr1.setAttribute("AN", "1,NC_000001.11")

  private val AltChr1 = new SAMSequenceRecord("1", 249250621)
  AltChr1.setAttribute("AN", "chr1,NC_000001.11")

  private val DictEmpty = new SAMSequenceDictionary(Collections.emptyList())
  private val DictChr1  = new SAMSequenceDictionary(Collections.singletonList(Chr1))
  private val DictAlt1  = new SAMSequenceDictionary(Collections.singletonList(AltChr1))


  "NcbiRefSeqSource" should "auto-map the accession to chr1 when given an empty sequence dictionary" in {
    val source = NcbiRefSeqGffSource(GffFile, includeXs=true, dict=DictEmpty)
    source should have size 7

    // Pseudo-gene should not have been included
    source.get("DDX11L1") shouldBe None

    // Check a micro-RNA for details
    val mir = source.get("MIR6859-1").value
    mir.loci should have size 1
    mir.loci.head.chrom shouldBe "chr1"
    mir.loci.head.start shouldBe 17369
    mir.loci.head.end   shouldBe 17436
    mir.loci.head.transcripts should have size 1
    mir.loci.head.transcripts.head.chrom shouldBe "chr1"
    mir.loci.head.transcripts.head.start shouldBe 17369
    mir.loci.head.transcripts.head.end   shouldBe 17436
    mir.loci.head.transcripts.head.cdsStart shouldBe None
    mir.loci.head.transcripts.head.cdsEnd shouldBe None
    mir.loci.head.transcripts.head.negativeStrand shouldBe true
    mir.loci.head.transcripts.head.exons shouldBe Seq(Exon(17369, 17436))

    // Check a lncRNA for specifics
    val lnc = source.get("MIR1302-2HG").value
    lnc.loci should have size 1
    lnc.loci.head.transcripts should have size 1
    lnc.loci.head.transcripts.head.chrom shouldBe "chr1"
    lnc.loci.head.transcripts.head.start shouldBe 29926
    lnc.loci.head.transcripts.head.end   shouldBe 31295
    lnc.loci.head.transcripts.head.cdsStart shouldBe None
    lnc.loci.head.transcripts.head.cdsEnd shouldBe None
    lnc.loci.head.transcripts.head.negativeStrand shouldBe false
    lnc.loci.head.transcripts.head.exons shouldBe Seq(Exon(29926, 30039), Exon(30564, 30667), Exon(30976, 31295))

    // Check a coding gene for
    val gene = source.get("OR4F5").value
    gene.loci should have size 1
    gene.loci.head.transcripts should have size 1
    gene.loci.head.transcripts.head.name shouldBe "NM_001005484.1"
    gene.loci.head.transcripts.head.chrom shouldBe "chr1"
    gene.loci.head.transcripts.head.start shouldBe 69091
    gene.loci.head.transcripts.head.end   shouldBe 70008
    gene.loci.head.transcripts.head.cdsStart.value shouldBe 69091
    gene.loci.head.transcripts.head.cdsEnd.value   shouldBe 70008
    gene.loci.head.transcripts.head.negativeStrand shouldBe false

    // Check a gene that has multiple transcripts
    val g2 = source.get("LOC100996442").value
    g2.loci should have size 1
    g2.loci.head.transcripts should have size 14
  }

  it should "still load all genes when given a dictionary that has all the used chroms in it" in {
    val source = NcbiRefSeqGffSource(GffFile, includeXs=true, dict=DictChr1)
    source should have size 7
    for (gene <- source; locus <- gene.loci; tx <- locus) {
      locus.chrom shouldBe "chr1"
      tx.chrom shouldBe "chr1"
    }
  }

  it should "map the chromosome name using the dictionary" in {
    val source = NcbiRefSeqGffSource(GffFile, includeXs=true, dict=DictAlt1)
    source should have size 7
    for (gene <- source; locus <- gene.loci; tx <- locus) {
      locus.chrom shouldBe "1"
      tx.chrom shouldBe "1"
    }
  }

  it should "exclude experimental transcripts (and genes with only exp transcripts)" in {
    val source = NcbiRefSeqGffSource(GffFile, includeXs=false, dict=DictChr1)
    source should have size 5

    for (gene <- source; locus <- gene.loci; tx <- locus) {
      tx.name.charAt(0) should not be 'X'
    }
  }
}
