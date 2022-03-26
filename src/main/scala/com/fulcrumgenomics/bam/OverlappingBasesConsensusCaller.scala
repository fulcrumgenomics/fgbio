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

import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.util.NumericTypes.PhredScore


/** Consensus calls the overlapping mapped bases in read pairs.
  *
  * This will iterate through the mapped bases that overlap between the read and mate in a read pair.  If the read and
  * mate agree at a given reference position, then read and mate base will not change and the base quality returned
  * is controlled by `maxQualOnAgreement`. If they disagree at a given reference position, then the base and quality
  * returned is controlled by `maskDisagreements`.
  *
  * @param maskDisagreements if the read and mate bases disagree at a given reference position, true to mask (make
  *                          'N') the read and mate bases, otherwise pick the base with the highest base quality and
  *                          return a base quality that's the difference between the higher and lower base quality.
  * @param maxQualOnAgreement    if the read and mate bases agree at a given reference position, true to for the
  *                              resulting base quality to be the maximum base quality, otherwise the sum of the base
  *                              qualities.
  */
class OverlappingBasesConsensusCaller(maskDisagreements: Boolean = false,
                                      maxQualOnAgreement: Boolean = false) {
  private val NoCall: Byte = 'N'.toByte
  private val NoCallQual: PhredScore = PhredScore.MinValue

  private val r1BasesBuilder = Array.newBuilder[Byte]
  private val r1QualsBuilder = Array.newBuilder[PhredScore]
  private val r2BasesBuilder = Array.newBuilder[Byte]
  private val r2QualsBuilder = Array.newBuilder[PhredScore]
  private val r1Builders     = Seq(r1BasesBuilder, r1QualsBuilder)
  private val r2Builders     = Seq(r2BasesBuilder, r2QualsBuilder)

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

    // Clear and resize the builders
    r1Builders.foreach { builder =>
      builder.clear()
      builder.sizeHint(r1.length)
    }
    r2Builders.foreach { builder =>
      builder.clear()
      builder.sizeHint(r2.length)
    }

    val r1Bases = r1.bases
    val r1Quals = r1.quals
    val r2Bases = r2.bases
    val r2Quals = r2.quals

    // Initialize the iters
    val r1Iter        = new MateOverlappingReadAndRefPosIterator(r1, r2).buffered
    val r2Iter        = new MateOverlappingReadAndRefPosIterator(r2, r1).buffered
    var r1LastReadPos = r1Iter.head.readPos - 1
    var r2LastReadPos = r2Iter.head.readPos - 1
    var overlappingBases = 0

    // add all the bases prior to the overlapping bases
    r1BasesBuilder.addAll(r1Bases.take(r1LastReadPos))
    r1QualsBuilder.addAll(r1Quals.take(r1LastReadPos))
    r2BasesBuilder.addAll(r2Bases.take(r2LastReadPos))
    r2QualsBuilder.addAll(r2Quals.take(r2LastReadPos))

    // Walk through the iterators by reference position
    while (r1Iter.hasNext && r2Iter.hasNext) {
      val r1Head = r1Iter.head
      val r2Head = r2Iter.head
      if (r1Head.refPos < r2Head.refPos) {
        r1BasesBuilder.addAll(r1Bases.slice(from=r1LastReadPos + 1, until=r1Head.readPos + 1))
        r1QualsBuilder.addAll(r1Quals.slice(from=r1LastReadPos + 1, until=r1Head.readPos + 1))
        r1LastReadPos = r1Head.readPos
        r1Iter.next()
      }
      else if (r2Head.refPos < r1Head.refPos) {
        r2BasesBuilder.addAll(r2Bases.slice(from=r2LastReadPos + 1, until=r2Head.readPos + 1))
        r2QualsBuilder.addAll(r2Quals.slice(from=r2LastReadPos + 1, until=r2Head.readPos + 1))
        r2LastReadPos = r2Head.readPos
        r2Iter.next()
      }
      else { // both reads are mapped to the same reference bases, so consensus call
        // add read bases from insertions
        r1BasesBuilder.addAll(r1Bases.slice(from = r1LastReadPos + 1, until = r1Head.readPos))
        r1QualsBuilder.addAll(r1Quals.slice(from = r1LastReadPos + 1, until = r1Head.readPos))
        r2BasesBuilder.addAll(r2Bases.slice(from = r2LastReadPos + 1, until = r2Head.readPos))
        r2QualsBuilder.addAll(r2Quals.slice(from = r2LastReadPos + 1, until = r2Head.readPos))

        overlappingBases += 1

        // consensus call current
        // TODO: skip consensus calling if either of the bases is an Nocall
        val (base, qual) = consensusCall(
          base1 = r1Bases(r1Head.readPos - 1),
          qual1 = r1Quals(r1Head.readPos - 1),
          base2 = r2Bases(r2Head.readPos - 1),
          qual2 = r2Quals(r2Head.readPos - 1)
        )

        // Only add the consensus call if we aren't masking disagreements or its a no-call
        if (maskDisagreements && base != NoCall) {
          r1BasesBuilder.addOne(r1Bases(r1LastReadPos))
          r1QualsBuilder.addOne(r1Quals(r1LastReadPos))
          r2BasesBuilder.addOne(r2Bases(r2LastReadPos))
          r2QualsBuilder.addOne(r2Quals(r2LastReadPos))
        } else {
          r1BasesBuilder.addOne(base)
          r1QualsBuilder.addOne(qual)
          r2BasesBuilder.addOne(base)
          r2QualsBuilder.addOne(qual)
        }

        // consume the current
        r1LastReadPos = r1Head.readPos
        r2LastReadPos = r2Head.readPos
        r1Iter.next()
        r2Iter.next()
      }
    }

    // add any remaining read and reference bases
    r1BasesBuilder.addAll(r1Bases.drop(r1LastReadPos))
    r1QualsBuilder.addAll(r1Quals.drop(r1LastReadPos))
    r2BasesBuilder.addAll(r2Bases.drop(r2LastReadPos))
    r2QualsBuilder.addAll(r2Quals.drop(r2LastReadPos))

    require(r1BasesBuilder.length == r1Bases.length)
    require(r2BasesBuilder.length == r2Bases.length)
    require(r1QualsBuilder.length == r1Quals.length)
    require(r2QualsBuilder.length == r2Quals.length)

    r1.bases = r1BasesBuilder.result()
    r1.quals = r1QualsBuilder.result()
    r2.bases = r2BasesBuilder.result()
    r2.quals = r2QualsBuilder.result()

    val r1Corrected = r1Bases.zip(r1.bases).count { case (left, right) => left != right }
    val r2Corrected = r2Bases.zip(r2.bases).count { case (left, right) => left != right }

    CorrectionStats(overlappingBases, r1Corrected, r2Corrected)
  }

  private def consensusCall(base1: Byte, qual1: PhredScore, base2: Byte, qual2: PhredScore): (Byte, PhredScore) = {
    // Capture the raw consensus base prior to masking it to N, so that we can compute
    // errors vs. the actually called base.
    val (rawBase: Byte, rawQual: PhredScore) = {
      if      (base1 == base2 && maxQualOnAgreement)   (base1,  Math.max(qual1, qual2).toByte) // use the maximum base quality
      else if (base1 == base2 && !maxQualOnAgreement)  (base1,  PhredScore.cap(qual1 + qual2)) // use the sum of base qualities
      else if (maskDisagreements)                      (NoCall, NoCallQual)                    // disagreements are no-calls
      else if (qual1 > qual2)                          (base1,  PhredScore.cap(qual1 - qual2))
      else if (qual2 > qual1)                          (base2,  PhredScore.cap(qual2 - qual1))
      else                                             (base1,  PhredScore.MinValue)
    }

    // Then mask it if appropriate
    if (base1 == NoCall || base2 == NoCall || rawQual == PhredScore.MinValue) (NoCall, NoCallQual) else (rawBase, rawQual)
  }
}

object OverlappingBasesConsensusCaller {
  /** Returns an iterator over records in the SAM source that are in query group order, where overlapping read pairs
    * are consensus called.  The input will be re-sorted if not already query sorted or grouped.  Statistics for how
    * may bases and templates had overlaps and were modified are logged to the given logger.
    *
    * @param in the input [[SamSource]] from which to read
    * @param logger the logger, to which output statistics are written
    * @param maskDisagreements see [[OverlappingBasesConsensusCaller]]
    * @param maxQualOnAgreement see [[OverlappingBasesConsensusCaller]]
    */
  def iterator(in: SamSource,
               logger: Logger,
               maskDisagreements: Boolean = false,
               maxQualOnAgreement: Boolean = false): Iterator[SamRecord] = {
    val templateMetric = CallOverlappingConsensusBasesMetric(kind=CountKind.Templates)
    val basesMetric    = CallOverlappingConsensusBasesMetric(kind=CountKind.Bases)
    val caller         = new OverlappingBasesConsensusCaller(
      maskDisagreements  = maskDisagreements,
      maxQualOnAgreement = maxQualOnAgreement
    )
    val templateIterator = Bams.templateIterator(in=in).flatMap { template =>
      caller.call(template);
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

