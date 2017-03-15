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

package com.fulcrumgenomics.vcf.filtration

import com.fulcrumgenomics.bam.{BaseEntry, Pileup, PileupBuilder, PileupEntry}
import com.fulcrumgenomics.testing.SamRecordSetBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec, VariantContextSetBuilder}
import htsjdk.variant.variantcontext.{Genotype, VariantContext}

class EndRepairArtifactLikelihoodFilterTest extends UnitSpec {
  private val A = 'A'.toByte
  private val C = 'C'.toByte
  private val G = 'G'.toByte
  private val T = 'T'.toByte

  /** Function to make a single well-formed genotype. */
  def singleGenotype(alleles: String*): Genotype = {
    val as = alleles.toList
    val builder = new VariantContextSetBuilder().addVariant(start=100, variantAlleles=as, genotypeAlleles=as)
    builder.head.getGenotype(0)
  }

  "appliesTo" should "return false for any event that is not a SNP" in {
    val filter  = new EndRepairArtifactLikelihoodFilter()
    filter.appliesTo(singleGenotype("A",  "AT")) shouldBe false // Insertion
    filter.appliesTo(singleGenotype("TA", "A"))  shouldBe false // Deletion
    filter.appliesTo(singleGenotype("AT", "GC")) shouldBe false // MNP
  }

  it should "return true for *>A and *>T and false for all other SNVs " in {
    val filter  = new EndRepairArtifactLikelihoodFilter()
    filter.appliesTo(singleGenotype("A",  "C")) shouldBe false
    filter.appliesTo(singleGenotype("A",  "G")) shouldBe false
    filter.appliesTo(singleGenotype("A",  "T")) shouldBe true
    filter.appliesTo(singleGenotype("C",  "A")) shouldBe true
    filter.appliesTo(singleGenotype("C",  "G")) shouldBe false
    filter.appliesTo(singleGenotype("C",  "T")) shouldBe true
    filter.appliesTo(singleGenotype("G",  "A")) shouldBe true
    filter.appliesTo(singleGenotype("G",  "C")) shouldBe false
    filter.appliesTo(singleGenotype("G",  "T")) shouldBe true
    filter.appliesTo(singleGenotype("T",  "A")) shouldBe true
    filter.appliesTo(singleGenotype("T",  "C")) shouldBe false
    filter.appliesTo(singleGenotype("T",  "G")) shouldBe false
  }

  "isArtifactCongruent" should "return false for any base in a read that is not near the ends" in {
    val filter  = new EndRepairArtifactLikelihoodFilter(distance=5)
    val recs    = new SamRecordSetBuilder(readLength=50).addPair(start1=101, start2=101)
    recs.foreach(_.setReadString("ACACA"*10))

    for (r <- recs; pos <- Range.inclusive(106, 145)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe false
    }
  }

  it should "return false for any base at the ends of the insert that is not allele-matched" in {
    val filter = new EndRepairArtifactLikelihoodFilter(distance=5)
    val recs   = new SamRecordSetBuilder(readLength=50).addPair(start1=101, start2=101)
    recs.foreach(_.setReadString("CACAC" + ("N"*40) + "GTGTG"))

    for (r <- recs; pos <- Range.inclusive(101, 105)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe false
    }

    for (r <- recs; pos <- Range.inclusive(146, 150)) {
      filter.isArtifactCongruent(G, T, BaseEntry(r, pos-101), pos) shouldBe false
    }
  }

  it should "return true for bases at the ends of the insert that are allele-matched" in {
    val filter = new EndRepairArtifactLikelihoodFilter(distance=5)
    val recs   = new SamRecordSetBuilder(readLength=50).addPair(start1=101, start2=101)
    recs.foreach(_.setReadString("GTGTG" + ("N" * 40) + "CACAC"))

    for (r <- recs; pos <- Range.inclusive(101, 105)) {
      filter.isArtifactCongruent(G, T, BaseEntry(r, pos-101), pos) shouldBe true
    }
    for (r <- recs; pos <- Range.inclusive(146, 150)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe true
    }
  }

  "priors" should "return values that prefer artifacts at low MAFs and vice versa" in {
    val filter = new EndRepairArtifactLikelihoodFilter(distance=5)
    val pileup = Pileup("chr1", 0, 100, Seq.empty[PileupEntry])
    val gt     = singleGenotype("A", "T")

    filter.priors(pileup, maf=0.01) shouldBe (0.0004, 0.9996)

    Range.inclusive(1, 100).sliding(2).foreach { case Seq(maf1, maf2) =>
      val (pMut1, pArt1) = filter.priors(pileup, maf1/100.0)
      val (pMut2, pArt2) = filter.priors(pileup, maf2/100.0)

      pMut1 should be <= pMut2
      pArt1 should be >= pArt2

      Seq(pMut1, pMut2, pArt1, pArt2).foreach { p =>
        p should be  > 0.0
        p should be <= 1.0
      }
    }
  }

  "annotations" should "compute a non-significant p-value when data is distributed throughout the reads" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new EndRepairArtifactLikelihoodFilter(distance=5)
    val builder = new SamRecordSetBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("G"*50))
    }

    // And a bit of alt allele
    for (start <- Range.inclusive(1, 25); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("T"*50))
    }

    val refName     = builder.dict.getSequence(0).getSequenceName
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.InfoKey) shouldBe true
    annotations(filter.InfoKey).asInstanceOf[Double] should be > 0.5
  }

  "it" should "compute a significant p-value when data is heavily biased" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new EndRepairArtifactLikelihoodFilter(distance=5)
    val builder = new SamRecordSetBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele at various positions
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 5); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("G"*50))
    }

    // And a bit of alt allele, all at the start of the read
    for (start <- Range.inclusive(21, 25); strand <- Seq(Plus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("T"*50))
    }

    val refName     = builder.dict.getSequence(0).getSequenceName
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.InfoKey) shouldBe true
    annotations(filter.InfoKey).asInstanceOf[Double] should be < 1e-6
  }

  "it" should "compute an intermediate p-value when data is heavily biased for both ref and alt" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new EndRepairArtifactLikelihoodFilter(distance=5)
    val builder = new SamRecordSetBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele, 5/6ths in the first 5bp of the read
    for (start <- Range.inclusive(20, 25); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("G"*50))
    }

    // And a bit of alt allele, all at the start of the read
    for (start <- Range.inclusive(21, 25); strand <- Seq(Plus)) {
      builder.addFrag(start=start, strand=strand).foreach(r => r.setReadString("T"*50))
    }

    val refName     = builder.dict.getSequence(0).getSequenceName
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.InfoKey) shouldBe true
    annotations(filter.InfoKey).asInstanceOf[Double] should be < 1e-4
  }

  "filters" should "only set the filter if a threshold was provide and the pvalue is <= threshold" in {
    val filterNoThreshold   = new EndRepairArtifactLikelihoodFilter(distance=5, pValueThreshold=None)
    val filterWithThreshold = new EndRepairArtifactLikelihoodFilter(distance=5, pValueThreshold=Some(0.001))
    def ann(pvalue: Double) = Map(filterNoThreshold.InfoKey -> pvalue)

    filterNoThreshold.filters(ann(1   )) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-3)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-4)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-5)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(0   )) should contain theSameElementsAs Seq()

    filterWithThreshold.filters(ann(1      )) should contain theSameElementsAs Seq()
    filterWithThreshold.filters(ann(0.01   )) should contain theSameElementsAs Seq()
    filterWithThreshold.filters(ann(0.001  )) should contain theSameElementsAs Seq(filterWithThreshold.FilterKey)
    filterWithThreshold.filters(ann(0.00099)) should contain theSameElementsAs Seq(filterWithThreshold.FilterKey)
    filterWithThreshold.filters(ann(1e-20  )) should contain theSameElementsAs Seq(filterWithThreshold.FilterKey)
    filterWithThreshold.filters(ann(0      )) should contain theSameElementsAs Seq(filterWithThreshold.FilterKey)
  }
}
