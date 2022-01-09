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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.{BaseEntry, PileupBuilder}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.vcf.api.{AlleleSet, Genotype}

class ReadEndSomaticVariantFilterTest extends UnitSpec {
  private val A = 'A'.toByte
  private val C = 'C'.toByte
  private val G = 'G'.toByte
  private val T = 'T'.toByte

  /** Function to make a single well-formed genotype. */
  def singleGenotype(alleles: String*): Genotype = {
    val as = AlleleSet(alleles:_*)
    Genotype(as, "s1", as.toIndexedSeq)
  }

  "filters" should "only set the filter if a threshold was provide and the pvalue is <= threshold" in {
    val filterNoThreshold   = new ATailingArtifactLikelihoodFilter(distance=5, pValueThreshold=None)
    val filterWithThreshold = new ATailingArtifactLikelihoodFilter(distance=5, pValueThreshold=Some(0.001))
    def ann(pvalue: Double) = Map(filterNoThreshold.readEndInfoLine.id -> pvalue)

    filterNoThreshold.filters(ann(1   )) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-3)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-4)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(1e-5)) should contain theSameElementsAs Seq()
    filterNoThreshold.filters(ann(0   )) should contain theSameElementsAs Seq()

    filterWithThreshold.filters(ann(1      )) should contain theSameElementsAs Seq()
    filterWithThreshold.filters(ann(0.01   )) should contain theSameElementsAs Seq()
    filterWithThreshold.filters(ann(0.001  )) should contain theSameElementsAs Seq(filterWithThreshold.readEndFilterLine.id)
    filterWithThreshold.filters(ann(0.00099)) should contain theSameElementsAs Seq(filterWithThreshold.readEndFilterLine.id)
    filterWithThreshold.filters(ann(1e-20  )) should contain theSameElementsAs Seq(filterWithThreshold.readEndFilterLine.id)
    filterWithThreshold.filters(ann(0      )) should contain theSameElementsAs Seq(filterWithThreshold.readEndFilterLine.id)
  }

  "EndRepairFillInArtifactLikelihoodFilter.appliesTo" should "return false for any event that is not a SNP" in {
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance = 15)
    filter.appliesTo(singleGenotype("G",  "GT")) shouldBe false // Insertion
    filter.appliesTo(singleGenotype("GT", "G"))  shouldBe false // Deletion
    filter.appliesTo(singleGenotype("CT", "GC")) shouldBe false // MNP
  }

  it should "return false for homozygous genotypes" in {
    val mg = Genotype(AlleleSet("A"), "s1", AlleleSet("A", "A").toIndexedSeq)
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance = 15)
    filter.appliesTo(mg) shouldBe false
  }

  it should "return false for events with multiple alt alleles if any alts are not SNVs" in {
    val mg = Genotype(AlleleSet("A", "G"), "s1", AlleleSet("A", "G", "TT").toIndexedSeq)
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance = 15)
    filter.appliesTo(mg) shouldBe false
  }

  it should "return true for any event that is a SNP" in {
    val bases = Seq("A", "C", "G", "T")
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance = 15)
    for (ref <- bases; alt <- bases; if ref != alt) yield filter.appliesTo(singleGenotype(ref, alt)) shouldBe true
  }

  it should "return true for events with multiple alt alleles if they are all SNVs" in {
    val mg = Genotype(AlleleSet("A", "G"), "s1", AlleleSet("A", "G", "T").toIndexedSeq)
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance = 15)
    filter.appliesTo(mg) shouldBe true
  }

  "EndRepairFillInArtifactLikelihoodFilter.isArtifactCongruent" should "return false for any base that is not within the defined read end" in {
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance=15)
    val recs    = new SamBuilder(readLength=50).addPair(start1=101, start2=101, bases1="ACACA"*10, bases2="ACACA"*10)

    for (r <- recs; pos <- Range.inclusive(116, 135)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe false
    }
  }

  it should "return true for any base within the defined read end" in {
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance=15)
    val recs    = new SamBuilder(readLength=50).addPair(start1=101, start2=101, bases1="ACACA"*10, bases2="ACACA"*10)

    for (r <- recs; pos <- Range.inclusive(101, 115)) {
      filter.isArtifactCongruent(G, T, BaseEntry(r, pos-101), pos) shouldBe true
    }
    for (r <- recs; pos <- Range.inclusive(136, 150)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe true
    }
  }

  "EndRepairFillInArtifactLikelihoodFilter.annotations" should "compute a non-significant p-value when data is distributed throughout the reads" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance=15)
    val builder = new SamBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="G"*50)
    }

    // And a bit of alt allele
    for (start <- Range.inclusive(1, 25); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="T"*50)
    }

    val refName     = builder.dict(0).name
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.readEndInfoLine.id) shouldBe true
    annotations(filter.readEndInfoLine.id).asInstanceOf[Double] should be > 0.5
  }

  it should "compute a significant p-value when data is heavily biased" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new EndRepairFillInArtifactLikelihoodFilter(distance=15)
    val builder = new SamBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele at various positions
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 5); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="G"*50)
    }

    // And a bit of alt allele, all at the start of the read
    for (start <- Range.inclusive(11, 25); strand <- Seq(Plus)) {
      builder.addFrag(start=start, strand=strand, bases="T"*50)
    }

    val refName     = builder.dict(0).name
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.readEndInfoLine.id) shouldBe true
    annotations(filter.readEndInfoLine.id).asInstanceOf[Double] should be < 1e-6
  }


  "ATailingArtifactLikelihoodFilter.appliesTo" should "return false for any event that is not a SNP" in {
    val filter  = new ATailingArtifactLikelihoodFilter()
    filter.appliesTo(singleGenotype("A",  "AT")) shouldBe false // Insertion
    filter.appliesTo(singleGenotype("TA", "A"))  shouldBe false // Deletion
    filter.appliesTo(singleGenotype("AT", "GC")) shouldBe false // MNP
  }

  it should "return true for *>A and *>T and false for all other SNVs " in {
    val filter  = new ATailingArtifactLikelihoodFilter()
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

  "ATailingArtifactLikelihoodFilter.isArtifactCongruent" should "return false for any base in a read that is not near the ends" in {
    val filter  = new ATailingArtifactLikelihoodFilter(distance=5)
    val recs    = new SamBuilder(readLength=50).addPair(start1=101, start2=101, bases1="ACACA"*10, bases2="ACACA"*10)

    for (r <- recs; pos <- Range.inclusive(106, 145)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe false
    }
  }

  it should "return false for any base at the ends of the insert that is not allele-matched" in {
    val filter = new ATailingArtifactLikelihoodFilter(distance=5)
    val recs   = new SamBuilder(readLength=50).addPair(start1=101, start2=101)
    recs.foreach(_.bases = "CACAC" + ("N" * 40) + "GTGTG")

    for (r <- recs; pos <- Range.inclusive(101, 105)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe false
    }

    for (r <- recs; pos <- Range.inclusive(146, 150)) {
      filter.isArtifactCongruent(G, T, BaseEntry(r, pos-101), pos) shouldBe false
    }
  }

  it should "return true for bases at the ends of the insert that are allele-matched" in {
    val filter = new ATailingArtifactLikelihoodFilter(distance=5)
    val recs   = new SamBuilder(readLength=50).addPair(start1=101, start2=101)
    recs.foreach(_.bases = "GTGTG" + ("N" * 40) + "CACAC")

    for (r <- recs; pos <- Range.inclusive(101, 105)) {
      filter.isArtifactCongruent(G, T, BaseEntry(r, pos-101), pos) shouldBe true
    }
    for (r <- recs; pos <- Range.inclusive(146, 150)) {
      filter.isArtifactCongruent(C, A, BaseEntry(r, pos-101), pos) shouldBe true
    }
  }

  "ATailingArtifactLikelihoodFilter.annotations" should "compute a non-significant p-value when data is distributed throughout the reads" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new ATailingArtifactLikelihoodFilter(distance=5)
    val builder = new SamBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="G"*50)
    }

    // And a bit of alt allele
    for (start <- Range.inclusive(1, 25); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="T"*50)
    }

    val refName     = builder.dict(0).name
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.readEndInfoLine.id) shouldBe true
    annotations(filter.readEndInfoLine.id).asInstanceOf[Double] should be > 0.5
  }

  it should "compute a significant p-value when data is heavily biased" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new ATailingArtifactLikelihoodFilter(distance=5)
    val builder = new SamBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele at various positions
    for (start <- Range.inclusive(1, 25); i <- Range.inclusive(1, 5); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="G"*50)
    }

    // And a bit of alt allele, all at the start of the read
    for (start <- Range.inclusive(21, 25); strand <- Seq(Plus)) {
      builder.addFrag(start=start, strand=strand, bases="T"*50)
    }

    val refName     = builder.dict(0).name
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.readEndInfoLine.id) shouldBe true
    annotations(filter.readEndInfoLine.id).asInstanceOf[Double] should be < 1e-6
  }

  it should "compute an intermediate p-value when data is heavily biased for both ref and alt" in {
    // Simulate a G>T non-artifact at bp 25
    val filter  = new ATailingArtifactLikelihoodFilter(distance=5)
    val builder = new SamBuilder(readLength=50, baseQuality=40)

    // Add a ton of reference allele, 5/6ths in the first 5bp of the read
    for (start <- Range.inclusive(20, 25); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
      builder.addFrag(start=start, strand=strand, bases="G"*50)
    }

    // And a bit of alt allele, all at the start of the read
    for (start <- Range.inclusive(21, 25); strand <- Seq(Plus)) {
      builder.addFrag(start=start, strand=strand, bases="T"*50)
    }

    val refName     = builder.dict(0).name
    val pile        = new PileupBuilder(dict=builder.dict, mappedPairsOnly=false).build(builder.iterator, refName, 25)
    val annotations = filter.annotations(pile, singleGenotype("G", "T"))
    annotations.contains(filter.readEndInfoLine.id) shouldBe true
    annotations(filter.readEndInfoLine.id).asInstanceOf[Double] should be < 1e-4
  }

  it should "not throw an exception if there is no alt allele coverage or the pileup is empty" in {
    val filter  = new ATailingArtifactLikelihoodFilter(distance=10)
    val builder = new SamBuilder(readLength=50, baseQuality=40)
    val refName = builder.dict(0).name

    { // Test with just an empty pileup
      val pile        = new PileupBuilder(dict = builder.dict, mappedPairsOnly = false).build(builder.iterator, refName, 25)
      val annotations = filter.annotations(pile, singleGenotype("G", "T"))
      annotations.contains(filter.readEndInfoLine.id) shouldBe true
    }

    { // And again with just a bunch of ref allele
      for (start <- Range.inclusive(1, 20); i <- Range.inclusive(1, 10); strand <- Seq(Plus, Minus)) {
        builder.addFrag(start = start, strand = strand, bases = "G" * 50)
      }
      val pile        = new PileupBuilder(dict = builder.dict, mappedPairsOnly = false).build(builder.iterator, refName, 25)
      val annotations = filter.annotations(pile, singleGenotype("G", "T"))
      annotations.contains(filter.readEndInfoLine.id) shouldBe true
    }
  }
}
