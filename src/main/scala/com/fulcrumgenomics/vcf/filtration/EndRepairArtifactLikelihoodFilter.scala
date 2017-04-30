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

import java.lang.Math.{min, pow}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.{Bams, BaseEntry, Pileup, PileupEntry}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.util.NumericTypes.{LogProbability => LnProb}
import htsjdk.samtools.util.SequenceUtil
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.vcf.{VCFFilterHeaderLine, VCFInfoHeaderLine}
/**
  * Filter that examines SNVs with the alt allele being A or T to see if they are likely the result
  * of aberrant end-repair and A-base addition.
  *
  * Occurs when end-repair is performed in a way that chews back single-stranded 3' overhangs
  * to create a blunt end, but over-digests and ends up creating a recessed 3' end which subsequently
  * (and incorrectly) gets filled in with one or more As during the A-base addition step. The result
  * is one or more base substitutions to A at the 3' end of molecules, which after PCR also shows up
  * as substitutions to T at the 5' end.
  *
  * @param distance The distance from the end of the template that defines the region where this
  *                 artifact may occur
  * @param pValueThreshold the pvalue at and below which mutations are filtered as being likely
  *                        caused by the artifact
  */
class EndRepairArtifactLikelihoodFilter(val distance: Int = 2,
                                        pValueThreshold: Option[Double] = None) extends SomaticVariantFilter with LazyLogging {
  val InfoKey           = "ERAP"
  val InfoDescription   = "P-Value for the End-Repair Artifact test."
  val FilterKey         = "EndRepairArtifact"
  val FilterDescription = "Variant is likely an artifact caused by end-repair and A-base addition."

  override val VcfInfoLines: Traversable[VCFInfoHeaderLine]     = Seq(vcfInfoLine(InfoKey, InfoDescription))
  override val VcfFilterLines: Traversable[VCFFilterHeaderLine] = Seq(vcfFilterLine(FilterKey, FilterDescription))

  /** Only applies to het SNVs where the alt allele is A or T. */
  override def appliesTo(gt: Genotype): Boolean = {
    gt.isHet &&
      !gt.isHetNonRef &&
      gt.getAlleles.forall(_.length() == 1) &&
      gt.getAlleles.filter(_.isNonReference).exists(a => a.getBaseString.equalsIgnoreCase("A") || a.getBaseString.equalsIgnoreCase("T"))
  }

  /** Calculates the p-value and returns it in the map of annotations. */
  override def annotations(pileup: Pileup[PileupEntry], gt: Genotype): Map[String,Any] = {
    if (!appliesTo(gt)) {
      Map.empty
    }
    else {
      val refAllele = SequenceUtil.upperCase(gt.getAlleles.filter(_.isReference).map(_.getBases()(0)).next())
      val altAllele = SequenceUtil.upperCase(gt.getAlleles.filter(_.isNonReference).map(_.getBases()(0)).next())
      annotations(pileup, refAllele, altAllele)
    }
  }

  /**
    * Workhorse method that does the majority of the calculation of the p-value.
    */
  private[filtration] def annotations(pileup: Pileup[PileupEntry], refAllele: Byte, altAllele: Byte) = {
    var (refGood, refBad, altGood, altBad) = (0, 0, 0, 0)
    pileup.baseIterator.filter(_.base == refAllele).foreach { rec =>
      if (isArtifactCongruent(refAllele, altAllele, rec, pileup.pos)) refBad += 1
      else refGood += 1
    }

    val lnFractionCongruent   = LnProb.toLogProbability(refBad / (refGood + refBad).toDouble)
    val lnFractionIncongruent = LnProb.toLogProbability(refGood / (refGood + refBad).toDouble)
    var llMutation = 1.0
    var llArtifact = 1.0

    // Now we can iterate over the alt bases and determine their likelihood under the two models
    pileup.baseIterator.filter(_.base == altAllele).foreach { rec =>
      val congruent    = isArtifactCongruent(refAllele, altAllele, rec, pileup.pos)
      val logProbError = LnProb.fromPhredScore(rec.qual)
      val logProbRight = LnProb.not(logProbError)

      if (congruent) {
        altBad += 1
        llArtifact += logProbRight
        llMutation += LnProb.or(LnProb.and(logProbRight, lnFractionCongruent), LnProb.and(logProbError, lnFractionIncongruent))
      }
      else {
        altGood += 1
        llArtifact += logProbError
        llMutation += LnProb.or(LnProb.and(logProbRight, lnFractionIncongruent), LnProb.and(logProbError, lnFractionCongruent))
      }
    }

    // Apply a prior based on the MAF
    val totalObs = altGood + altBad + refGood + refBad
    val maf      = if (totalObs > 0) (altGood + altBad) / totalObs.toDouble else 0
    val (priorMutation, priorArtifact) = priors(pileup, maf)
    llMutation += LnProb.toLogProbability(priorMutation)
    llArtifact += LnProb.toLogProbability(priorArtifact)

    val total = LnProb.or(llMutation, llArtifact)
    val pMutation = LnProb.expProb(LnProb.normalizeByLogProbability(llMutation, total))
    val pArtifact = LnProb.expProb(LnProb.normalizeByLogProbability(llArtifact, total))

    Map(InfoKey -> pMutation)
  }

  /** Calculates a pair of priors:
    *   The first is the prior that the genotype is the result of a true somatic mutations
    *   The second is the prior that the genotype is the result of the end-repair artifact process
    *
    * Default implementation gives a mutation prior of {{{(maf * 2)^2}}} capped at 0.9999 and an
    * artifact prior of 1 - mutation prior.  This has the effect of creating priors that are heavily
    * skewed towards the artifact at very low MAFs, and towards mutations as the MAF increases towards 0.5.
    */
  private[filtration] def priors(pileup: Pileup[PileupEntry], maf: Double) : (Double, Double) = {
    require(maf >= 0 && maf <= 1, s"Maf must be between 0 and 1. Observed maf ${maf} at ${pileup.refName}:${pileup.pos}.")

    val m = if (maf != 0) maf else {
      logger.warning(s"No alt allele observations found with selected cutoffs at: ${pileup.refName}:${pileup.pos}")
      1.0 / pileup.depth
    }

    val priorMutation = min(pow(m * 2, 2), 0.9999)
    val priorArtifact = 1 - priorMutation
    (priorMutation, priorArtifact)
  }

  /**
    * Returns true if a base in the pileup (either ref or alt) is congruent with the artifact.
    * I.e. the base is within the first or last few bases of the template _and_ the base as read
    * makes sense given the position.
    *
    * Given a ref allele and alt allele (alt will always be A or T when invoking this method) we infer
    * how the alleles are paired.  E.g. if we are give ref=C alt=A then we pair (C,A) and (G,T)
    * together and think about C and A in similar ways, and G and T in similar ways.
    *
    * Next, it's important to know if we are looking at an observation at the beginning or end of
    * the read (e.g. if the read goes all the way through the template), as this informs how the
    * artifact will present.  See the two examples below for a C>A artifact:
    *
    * Example 1: congruent with the artifact
    *     Ref(+)     CCCCCCCCCCCCCCCCCC
    *     R1(+)   5'-CCCCCCCCCCCCCCCCCa-3'
    *     R2(-)   3'-GGGGGGGGGGGGGGGGGt-5' (in BAM would also be CCCCCCCCCCCCCCCCCa)
    *
    * Example 2: incongruent with the artifact
    *     Ref(+)     CCCCCCCCCCCCCCCCCC
    *     R1(+)   5'-aCCCCCCCCCCCCCCCCC-3'
    *     R2(-)   3'-tGGGGGGGGGGGGGGGGG-5' (in BAM would also be aCCCCCCCCCCCCCCCCC)
    *
    * In Example 1 we see a C>A substitution at the end of the template, which R1 observes at a C>A substitution
    * _at the end of the read_ and R2 observes as a G>T substitution _at the beginning of the read_.  This is
    * congruent with the model of over-digestion of the 3' end and addition of one or more A nucleotides.
    *
    * In Example 2 we see a C>A substitution at the beginning of the template, which R1 observes at a C>A substitution
    * _at the beginning of the read_ and R2 observes as a G>T substitution _at the end of the read_.  This is
    * NOT congruent with the model of over-digestion of the 3' end and addition of one or more A nucleotides.
    *
    * Therefore for alt observations to be congruent we require that we observe a T (as read by the sequencer)
    * at the start of the read or an A (as read by the sequencer) at the end of the read (where the read spans
    * the whole template).
    *
    * Then for ref observations to be congruent we require an observation of the T-paired allele at the start of
    * the read or the A-paired allele at the end of the read - i.e. they are observations that if converted to
    * alt observations would be congruent.
    */
  private[filtration] def isArtifactCongruent(refAllele: Byte, altAllele: Byte, entry : BaseEntry, pileupPosition: Int): Boolean = {
    require(entry.base == refAllele || entry.base == altAllele, "Pileup base is neither ref or alt.")
    val isRef                = entry.base == refAllele
    val positionInRead       = entry.positionInReadInReadOrder
    val positionFromOtherEnd = Bams.positionFromOtherEndOfTemplate(entry.rec, pileupPosition)

    // Pick the position that's closest to an end, and define the congruent base by that end
    val (pos, congruentAlt) = {
      if (positionFromOtherEnd.exists(_ < positionInRead)) (positionFromOtherEnd.get, 'A')
      else                                                 (positionInRead, 'T')
    }

    // If we're not in the danger zone it can't be the artifact
    if (pos > this.distance) {
      false
    }
    else {
      val congruentRef = if (congruentAlt == altAllele) refAllele else SequenceUtil.complement(refAllele)
      entry.baseInReadOrientation == (if (isRef) congruentRef else congruentAlt)
    }
  }

  /** Applies a single filter if a) a threshold was provided at construction, b) the computed
    * pvalue is present in the annotations, and c) the pvalue <= threshold.
    */
  override def filters(annotations: Map[String, Any]): Traversable[String] = {
    (this.pValueThreshold, annotations.get(InfoKey).map(_.asInstanceOf[Double])) match {
      case (Some(threshold), Some(pvalue)) if pvalue <= threshold => List(FilterKey)
      case _ => Nil
    }
  }
}
