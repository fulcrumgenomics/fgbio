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

import java.util

import com.fulcrumgenomics.util.MathUtil
import com.fulcrumgenomics.util.NumericTypes._
import htsjdk.samtools.util.SequenceUtil

import scala.math.{max, min}

object ConsensusCaller {
  type Base = Byte

  private val NoCall = 'N'.toByte
  private val DnaBasesUpperCase: Array[Base] = Array('A', 'C', 'G', 'T').map(_.toByte)
  private val DnaBaseCount  = DnaBasesUpperCase.length

  /** Returns a copy of the DNA bases in the order they are used in all internal arrays. */
  def DnaBases: Array[Base] = DnaBasesUpperCase.clone()

}

/**
  * Generic consensus caller class that can be used to produce consensus base calls and qualities from
  * a pileup of raw base calls and qualities.
  *
  * The consensus caller sees the process of going from a DNA source molecule in its original, pristine, state
  * to a sequenced base as having three phases each with their own distinct error profiles:
  *   1) The phase whereby the source molecule is harvested (e.g. cells extracted and lysed) to the point where
  *   some kind of molecular identifier has been attached that will allow for identification of replicates that
  *   are generated from the same original source molecule.  Errors in this phase will be present in all copies
  *   of the molecule that are prepared for sequencing (except where a second error reverts the change).
  *
  *   2) Everything between phases 1 & 3! Generally including any sample preparation activities after a molecular
  *   identifier has been attached but prior to sequencing.  Errors introduced in this phase will be present in
  *   some fraction of the molecules available at sequencing.
  *
  *   3) Sequencing of the molecule (or clonal cluster of molecules) on a sequencer.  E.g. the process of base-by-base
  *   resynthesis and sequencing on an Illumina sequencer _after_ cluster amplification.  Errors in this phase are
  *   captured by the raw base quality scores from the sequencer.
  *
  * @param errorRatePreLabeling an estimate of the error rate (i.e. rate of base substitutions) caused by DNA damage
  *                             prior to any labeling of source molecules. Estimates the rate at which errors from
  *                             Phase 1 described above would be observed if the rest of the process were error-free.
  * @param errorRatePostLabeling an estimate of the error rate (i.e. rate of base substitutions) caused by DNA damage
  *                              post-labeling. Estimates the errors from Phase 2 which would be non-uniform across
  *                              replicates of the same source molecule.
  * @param rawBaseQualityShift shift the incoming raw base quality by this amount before use. Useful if the
  *                            raw base qualities are believed to be uniformly over-estimated. E.g. a value of 5
  *                            would cause an incoming base quality of 30 to be shifted to 25.
  * @param maxRawBaseQuality the maximum raw base quality to allow after any shift has been applied. Base qualities
  *                          higher than this value are capped at this value.
  */
class ConsensusCaller(errorRatePreLabeling:  PhredScore,
                      errorRatePostLabeling: PhredScore,
                      rawBaseQualityShift:   PhredScore = 0.toByte,
                      maxRawBaseQuality:     PhredScore = PhredScore.MaxValue
                      ) {
  import ConsensusCaller._

  assert(maxRawBaseQuality >= PhredScore.MinValue, s"maxBaseQuality must be >= ${PhredScore.MinValue}")
  assert(maxRawBaseQuality <= PhredScore.MaxValue, s"maxBaseQuality must be <= ${PhredScore.MaxValue}")

  private val LnErrorRatePreLabeling  = LogProbability.fromPhredScore(errorRatePreLabeling)
  private val LnErrorRatePostLabeling = LogProbability.fromPhredScore(errorRatePostLabeling)

  /**
    * An inner class for tracking the likelihoods for the consensus for a single base.
    */
  class ConsensusBaseBuilder {
    private val observations = new Array[Int](DnaBaseCount)
    private val likelihoods = new Array[LogProbability](DnaBaseCount)

    /** Resets the likelihoods to p=1 so that the builder can be re-used. */
    def reset(): Unit = {
      util.Arrays.fill(observations, 0)
      util.Arrays.fill(this.likelihoods, LnOne)
    }

    /** Adds a base and un-adjusted base quality to the consensus likelihoods. */
    def add(base: Base, qual: PhredScore): Unit = add(base, pError=phredToAdjustedLogProbError(qual), pTruth=phredToAdjustedLogProbCorrect(qual))

    /** Adds a base with adjusted error and truth probabilities to the consensus likelihoods. */
    def add(base: Base, pError: LogProbability, pTruth: LogProbability) = {
      val b = SequenceUtil.upperCase(base)
      if (b != 'N') {
        var i = 0
        while (i < DnaBaseCount) {
          val candidateBase = DnaBasesUpperCase(i)

          if (base == candidateBase) {
            likelihoods(i) += pTruth
            observations(i) += 1
          }
          else {
            likelihoods(i) += LogProbability.normalizeByScalar(pError, 3)
          }

          i += 1
        }
      }
    }

    /**
      * Returns the number of reads that contributed evidence to the consensus. The value is equal
      * to the number of times add() was called with non-ambiguous bases.
      */
    def contributions: Int = this.observations.sum

    /** Gets the number of observations of the base in question. */
    def observations(base: Base): Int = base match {
      case 'A' => this.observations(0)
      case 'C' => this.observations(1)
      case 'G' => this.observations(2)
      case 'T' => this.observations(3)
      case x   => throw new IllegalArgumentException("Unsupported base: " + x.toChar)
    }

    /** Call the consensus base and quality score given the current set of likelihoods. */
    def call() : (Base, PhredScore) = {
      // get the sum of the likelihoods
      // pick the base with the maximum posterior
      val lls  = likelihoods
      val likelihoodSum   = LogProbability.or(lls)
      val (maxLikelihood, maxLlIndex) = MathUtil.maxWithIndex(lls, requireUniqueMaximum=true)

      maxLlIndex match {
        case -1 =>
          (NoCall, PhredScore.MinValue) // Multiple equally likely bases
        case i  =>
          val maxPosterior    = LogProbability.normalizeByLogProbability(maxLikelihood, likelihoodSum)
          val pConsensusError = LogProbability.not(maxPosterior) // convert to probability of the called consensus being wrong

          // Factor in the pre-UMI error rate.
          // Pr(error) = Pr(any pre-UMI error AND correct consensus) + Pr(no pre-UMI error AND any error in consensus)
          //               + Pr(pre-UMI error AND error in consensus, that do not give us the correct bases)
          // The last term tries to capture the case where a pre-UMI error modifies the base (ex. A->C) but a sequencing
          // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
          val p = LogProbability.probabilityOfErrorTwoTrials(LnErrorRatePreLabeling, pConsensusError)
          val q = PhredScore.cap(PhredScore.fromLogProbability(p))
          val base = DnaBasesUpperCase(i)
          (base, q)
      }
    }

    /**
      * Gives the calculated likelihoods of the four bases given the data and the model. The likelihoods
      * returned factor in both the base observations and the probability of error prior to applying the
      * labels.
      */
    def logLikelihoods: Array[LogProbability] = {
      val lls             = likelihoods
      val likelihoodSum   = LogProbability.or(lls)
      val posteriors      = lls.map(l => LogProbability.normalizeByLogProbability(l, likelihoodSum))
      val errors          = posteriors.map(LogProbability.not)
      val errorsTwoTrials = errors.map(e => LogProbability.probabilityOfErrorTwoTrials(LnErrorRatePreLabeling, e))
      errorsTwoTrials.map(LogProbability.not)
    }
  }


  /** Pre-computes the the log-scale probabilities of an error for each a phred-scaled base quality from 0-127. */
  private val phredToAdjustedLogProbError: Array[LogProbability] = Range(0, Byte.MaxValue).toArray.map(q => {
    val aq = min(this.maxRawBaseQuality, max(PhredScore.MinValue, q - this.rawBaseQualityShift)).toByte
    val e1 = LogProbability.fromPhredScore(this.errorRatePostLabeling)
    val e2 = LogProbability.fromPhredScore(aq)
    LogProbability.probabilityOfErrorTwoTrials(e1, e2)
  })

  /** Pre-computes the the log-scale probabilities of an not an error for each a phred-scaled base quality from 0-127. */
  private val phredToAdjustedLogProbCorrect: Array[Double] = phredToAdjustedLogProbError.map(LogProbability.not)

  /**
    * Returns the adjusted probability of error given a base quality. This is the probability that the base was
    * sequenced incorrectly or that the template had the wrong base due to an error during sample-preparation
    * post-labeling (phase 2 described above).
    *
    * @param qual the raw base quality as a Phred scaled number
    * @return the LogProbability that the base is incorrect
    */
  def adjustedErrorProbability(qual: PhredScore): LogProbability = {
    if (qual < PhredScore.MinValue) throw new IllegalArgumentException(s"Cannot adjust score lower than ${PhredScore.MinValue}")
    else this.phredToAdjustedLogProbError(qual)
  }

  /**
    * Returns 1 - #adjustedErrorProbability
    *
    * @param qual the raw base quality as a Phred scaled number
    * @return the LogProbability that the base is incorrect
    */
  def adjustedTruthProbability(qual: PhredScore): LogProbability = {
    if (qual < PhredScore.MinValue) throw new IllegalArgumentException(s"Cannot adjust score lower than ${PhredScore.MinValue}")
    else this.phredToAdjustedLogProbCorrect(qual)
  }

  /** Returns a new builder that can be used to call the consensus at one or more sites serially. */
  def builder(): ConsensusBaseBuilder = new ConsensusBaseBuilder
}
