/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import enumeratum.EnumEntry
import htsjdk.samtools.util.SequenceUtil


/** Consensus calls the overlapping mapped bases in read pairs.
  *
  * This will iterate through the mapped bases that overlap between the read and mate in a read pair.  If the read and
  * mate agree at a given reference position, then read and mate base will not change and the base quality returned
  * is controlled by `agreementStrategy` (see [[AgreementStrategy]]). If they disagree at a given reference position,
  * then the base and quality returned is controlled by `disagreementStrategy` (see [[DisagreementStrategy]]).
  *
  * If either read base is a no-call (e.g. N), then the read bases for each read respectively will not be changed.
  *
  * @param agreementStrategy the strategy to use when the two bases agree.  This only affects the base qualities.
  * @param disagreementStrategy the disagreement strategy to use when the two bases disagree.  This will affect both the
  *                             bases and the qualities.
  */
class OverlappingBasesConsensusCaller(agreementStrategy: AgreementStrategy = AgreementStrategy.Consensus,
                                      disagreementStrategy: DisagreementStrategy = DisagreementStrategy.Consensus) {
  import OverlappingBasesConsensusCaller.NoCall

  /** Consensus calls the overlapping bases if and only if the template is a paired end where both ends map with at
    * least one base overlapping.
    *
    * @param template the template to potentially correct.
    * @return summary statistics about how many bases were examined and modified
    */
  def call(template: Template): CorrectionStats = (template.r1, template.r2) match {
    case (Some(r1), Some(r2)) if r1.mapped && r2.mapped && r1.matesOverlap.contains(true) => call(r1, r2)
    case _ => CorrectionStats()
  }

  /** Consensus calls the overlapping bases if and only if both ends map with at least one base overlapping.
    *
    * @param r1 the first read in the pair
    * @param r2 the second read in the pair
    * @return summary statistics about how many bases were examined and modified
    */
  def call(r1: SamRecord, r2: SamRecord): CorrectionStats = {
    require(r1.mapped && r2.mapped && r1.paired && r2.paired && r1.name == r2.name && r1.matesOverlap.contains(true))
    require(r1.firstOfPair && r2.secondOfPair)

    val iter             = new ReadMateAndRefPosIterator(rec=r1, mate=r2).buffered
    var overlappingBases = 0
    var r1Corrected      = 0
    var r2Corrected      = 0

    // Walk through the iterators by reference position
    iter.foreach { item =>
      val base1 = r1.bases(item.readPos - 1)
      val base2 = r2.bases(item.matePos - 1)
      if (!SequenceUtil.isNoCall(base1) && !SequenceUtil.isNoCall(base2)) {
        // consensus call
        val (newBase1, newQual1, newBase2, newQual2) = consensusCall(
          base1 = base1,
          qual1 = r1.quals(item.readPos - 1),
          base2 = base2,
          qual2 = r2.quals(item.matePos - 1)
        )

        if (newBase1 != NoCall && newBase1 != base1) r1Corrected += 1
        if (newBase2 != NoCall && newBase2 != base2) r2Corrected += 1

        r1.bases(item.readPos - 1) = newBase1
        r1.quals(item.readPos - 1) = newQual1
        r2.bases(item.matePos - 1) = newBase2
        r2.quals(item.matePos - 1) = newQual2

        overlappingBases += 1 // only one for the template
      }
    }

    CorrectionStats(overlappingBases, r1Corrected, r2Corrected)
  }

  private def consensusCall(base1: Byte, qual1: PhredScore,
                            base2: Byte, qual2: PhredScore): (Byte, PhredScore, Byte, PhredScore) = {
    if (base1 == base2) {
      val (q1, q2) = agreementStrategy.qual(qual1, qual2)
      (base1, q1, base2, q2)
    }
    else {
      disagreementStrategy.baseAndQual(base1, qual1, base2, qual2)
    }
  }
}

object OverlappingBasesConsensusCaller {
  val NoCall: Byte = 'N'.toByte
  val NoCallQual: PhredScore = PhredScore.MinValue

  /** Returns an iterator over records in the SAM source that are in query group order, where overlapping read pairs
    * are consensus called.  The input will be re-sorted if not already query sorted or grouped.  Statistics for how
    * may bases and templates had overlaps and were modified are logged to the given logger.
    *
    * @param in the input [[SamSource]] from which to read
    * @param logger the logger, to which output statistics are written
    * @param agreementStrategy the strategy to use when the two bases agree.  This only affects the base qualities.
    * @param disagreementStrategy the disagreement strategy to use when the two bases disagree.  This will affect both the
    *                             bases and the qualities.
    */
  def iterator(in: SamSource,
               logger: Logger,
               agreementStrategy: AgreementStrategy = AgreementStrategy.Consensus,
               disagreementStrategy: DisagreementStrategy = DisagreementStrategy.Consensus): Iterator[SamRecord] = {
    val templateMetric = CallOverlappingConsensusBasesMetric(kind=CountKind.Templates)
    val basesMetric    = CallOverlappingConsensusBasesMetric(kind=CountKind.Bases)
    val caller         = new OverlappingBasesConsensusCaller(
      agreementStrategy    = agreementStrategy,
      disagreementStrategy = disagreementStrategy
    )
    val templateIterator = Bams.templateIterator(in=in).flatMap { template =>
      caller.call(template)
      // update metrics
      templateMetric.total += 1
      basesMetric.total += template.primaryReads.map(_.length).sum
      template.allReads
    }
    new SelfClosingIterator(templateIterator, closer = () => {
      val pctTemplates = if (templateMetric.overlapping == 0) 0 else 100 * templateMetric.corrected / templateMetric.overlapping.toDouble
      val pctBases     = if (basesMetric.overlapping == 0) 0 else 100 * basesMetric.corrected / basesMetric.overlapping.toDouble
      logger.info("Consensus calling overlapping read pairs statistics:")
      logger.info(f"    ${templateMetric.overlapping}%,d overlapping templates")
      logger.info(f"    ${templateMetric.corrected}%,d corrected templates (${pctTemplates}%.2f%%)")
      logger.info(f"    ${basesMetric.overlapping}%,d overlapping bases")
      logger.info(f"    ${basesMetric.corrected}%,d corrected bases (${pctBases}%.2f%%)")
    })
  }
}

/** Statistics for consensus calling overlapping bases in a read pair
  *
  * @param overlappingBases the number of bases that overlap between R1 and R2
  * @param r1CorrectedBases the number of bases modified in R1
  * @param r2CorrectedBases the number of bases modified in R2
  */
case class CorrectionStats(overlappingBases: Int = 0, r1CorrectedBases: Int = 0, r2CorrectedBases: Int = 0)


/** Trait that entries in [[AgreementStrategy]] will extend. */
sealed trait AgreementStrategy extends EnumEntry {
  def qual(qual1: PhredScore, qual2: PhredScore): (PhredScore, PhredScore)
}


/** Enum to represent the strategy to consensus call when both bases agree. */
object AgreementStrategy extends FgBioEnum[AgreementStrategy] {
  /** Call the consensus base and return a new base quality that is the sum of the two base qualities. */
  case object Consensus extends AgreementStrategy {
    def qual(qual1: PhredScore, qual2: PhredScore): (PhredScore, PhredScore) = {
      val qual = PhredScore.cap(qual1 + qual2)
      (qual, qual)
    }
  }

  /** Call the consensus base and return a new base quality that is the maximum of the two base qualities. */
  case object MaxQual extends AgreementStrategy {
    def qual(qual1: PhredScore, qual2: PhredScore): (PhredScore, PhredScore) = {
      val qual = PhredScore.cap(Math.max(qual1, qual2))
      (qual, qual)
    }
  }

  /** Leave the bases and base qualities unchanged. */
  case object PassThrough extends AgreementStrategy {
    def qual(qual1: PhredScore, qual2: PhredScore): (PhredScore, PhredScore) = (qual1, qual2)
  }

  override def values: IndexedSeq[AgreementStrategy] = super.findValues
}

/** Trait that entries in [[DisagreementStrategy]] will extend. */
sealed trait DisagreementStrategy extends EnumEntry {
  def baseAndQual(base1: Byte, qual1: PhredScore, base2: Byte, qual2: PhredScore): (Byte, PhredScore, Byte, PhredScore)
}

/** Enum to represent the strategy to consensus call when both bases disagree.  Masking a base will make the base
  * an "N" with base quality "2". */
object DisagreementStrategy extends FgBioEnum[DisagreementStrategy] {
  import OverlappingBasesConsensusCaller.{NoCall, NoCallQual}

  /** Mask both bases. */
  case object MaskBoth extends DisagreementStrategy {
    def baseAndQual(base1: Byte, qual1: PhredScore, base2: Byte, qual2: PhredScore): (Byte, PhredScore, Byte, PhredScore) = {
      (NoCall, NoCallQual, NoCall, NoCallQual)
    }
  }

  /** Mask the base with the lowest base quality, with the other base unchanged.  If the base qualities are the same,
    * mask both bases. */
  case object MaskLowerQual extends DisagreementStrategy {
    def baseAndQual(base1: Byte, qual1: PhredScore, base2: Byte, qual2: PhredScore): (Byte, PhredScore, Byte, PhredScore) = {
      //                       r1-base  r1-qual     2-base  r2-qual
      if (qual1 > qual2)      (base1,   qual1,      NoCall, NoCallQual)
      else if (qual2 > qual1) (NoCall,  NoCallQual, base2,  qual2)
      else                    (NoCall,  NoCallQual, NoCall, NoCallQual)
    }
  }

  /** Consensus call the base.  If the base qualities are the same, mask both bases.  Otherwise, call the base
    * with the highest base quality and return a new base quality that is the difference between the highest and lowest
    * base quality. */
  case object Consensus extends DisagreementStrategy {
    def baseAndQual(base1: Byte, qual1: PhredScore, base2: Byte, qual2: PhredScore): (Byte, PhredScore, Byte, PhredScore) = {
      val (base, qual) = {
        if (qual1 > qual2)      (base1, PhredScore.cap(qual1 - qual2))
        else if (qual2 > qual1) (base2, PhredScore.cap(qual2 - qual1))
        else                    (NoCall, NoCallQual)
      }
      (base, qual, base, qual)
    }
  }

  override def values: IndexedSeq[DisagreementStrategy] = super.findValues
}