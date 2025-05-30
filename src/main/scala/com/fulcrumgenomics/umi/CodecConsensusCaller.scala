/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.bam.{ClippingMode, SamRecordClipper}
import com.fulcrumgenomics.umi.DuplexConsensusCaller.{FilterFragments, FilterMinReads}
import com.fulcrumgenomics.umi.UmiConsensusCaller.{ReadType, SourceRead}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.Sequences

/**
  * TODO: write docs
  *
  * @param readNamePrefix the prefix to apply to all consensus read names
  * @param readGroupId    the read group ID to apply to all created consensus reads
  * @param minInputBaseQuality the minimum input base quality score to use a raw read's base
  * @param errorRatePreUmi the estimated rate of errors in the DNA prior to attaching UMIs
  * @param errorRatePostUmi the estimated rate of errors in the DNA post attaching UMIs
  * @param minReadsPerStrand
  * @param maxReads
  * @param minDuplexLength
  */
class CodecConsensusCaller(readNamePrefix: String,
                           readGroupId: String = "A",
                           minInputBaseQuality: PhredScore = DuplexConsensusCaller.MinInputBaseQuality,
                           errorRatePreUmi: PhredScore = DuplexConsensusCaller.ErrorRatePreUmi,
                           errorRatePostUmi: PhredScore = DuplexConsensusCaller.ErrorRatePostUmi,
                           val minReadsPerStrand: Int = 1,
                           maxReadsPerStrand: Int = VanillaUmiConsensusCallerOptions.DefaultMaxReads,
                           minDuplexLength: Int = 1
                          ) extends DuplexConsensusCaller(
  readNamePrefix      = readNamePrefix,
  readGroupId         = readGroupId,
  minInputBaseQuality = minInputBaseQuality,
  qualityTrim         = false,
  errorRatePreUmi     = errorRatePreUmi,
  errorRatePostUmi    = errorRatePostUmi,
  minReads            = Seq(1, 1, 1),
  maxReadsPerStrand   = maxReadsPerStrand
) {
  require(this.minReadsPerStrand >= 1, "minReadsPerStrand must be at least 1.")

  private val FilterDuplexTooShort = s"Top/Bottom Strand Overlap < ${minDuplexLength} bases."
  private val FilterIndelError = "Indel error between strands in duplex overlap."
  private val clipper = new SamRecordClipper(ClippingMode.Hard, autoClipAttributes = false)

  /** Returns the MI tag as is because in CODEC there is no /A or /B unlike classic duplex-seq */
  override protected[umi] def sourceMoleculeId(rec: SamRecord): String = rec[String]("MI")

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SamRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  override protected def consensusSamRecordsFromSamRecords(recs: Seq[SamRecord]): Seq[SamRecord] = {
    val (pairs, frags) = recs.partition(_.paired)
    rejectRecords(frags, FilterFragments)

    if (pairs.isEmpty) {
      Nil
    }
    else {
      // Extract just the primary alignments and then clip where they extend past the ends of their mates
      // TODO: Remove the isFrPair here and handle chimeric reads or reads with one unmapped etc.
      val primaries = pairs.filterNot(r => r.secondary || r.supplementary).filter(_.isFrPair)
      primaries.groupBy(_.name).foreach { case (_, Seq(rec, mate)) => clipper.clipExtendingPastMateEnds(rec, mate) }

      val r1s = filterToMostCommonAlignment(primaries.filter(_.firstOfPair).map(toSourceReadForCodec))
      val r2s = filterToMostCommonAlignment(primaries.filter(_.secondOfPair).map(toSourceReadForCodec))
      if (r1s.length < this.minReadsPerStrand || r2s.length < this.minReadsPerStrand) {
        rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), FilterMinReads)
        Nil
      }
      else {
        // Calculate the max overlap between R1 and R2, and if it's too short then filter those reads
        val longestR1Alignment = r1s.flatMap(_.sam).maxBy(_.cigar.lengthOnTarget)
        val longestR2Alignment = r2s.flatMap(_.sam).maxBy(_.cigar.lengthOnTarget)

        val (longestPosAln, longestNegAln) =
          if (longestR1Alignment.positiveStrand) (longestR1Alignment, longestR2Alignment)
          else (longestR2Alignment, longestR1Alignment)

        val (overlapStart, overlapEnd) = (longestNegAln.start, longestPosAln.end)

        // If the overlap isn't long enough then reject the records and return no consensus reads
        if (overlapEnd - overlapStart + 1 < minDuplexLength) {
          rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), FilterDuplexTooShort)
          Nil
        }
        else if (
        // Check that the start and end of the overlap are in phase with one another in the longest alignments
          longestR1Alignment.readPosAtRefPos(overlapStart, returnLastBaseIfDeleted = true)
            - longestR2Alignment.readPosAtRefPos(overlapStart, returnLastBaseIfDeleted = true)
            !=
            longestR1Alignment.readPosAtRefPos(overlapEnd, returnLastBaseIfDeleted = true)
              - longestR2Alignment.readPosAtRefPos(overlapEnd, returnLastBaseIfDeleted = true)) {
          rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), FilterIndelError)
          Nil
        }
        else {
          // Call the single-stranded consensus reads
          val r1Consensus = this.ssCaller.consensusCall(r1s).getOrElse(unreachable("Failed to make SS consensus."))
          val r2Consensus = this.ssCaller.consensusCall(r2s).getOrElse(unreachable("Failed to make SS consensus."))

          // Re-Reverse the negative strand consensus so it's easier to line up with the other read
          if (longestR1Alignment.negativeStrand) r1Consensus.revcomp() else r2Consensus.revcomp()

          // Calculate the length of the consensus
          computeConsensusLength(longestPosAln, longestNegAln) match {
            case -1 =>
              rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), FilterIndelError)
              Nil
            case n if n < r1Consensus.length || n < r2Consensus.length =>
              // TODO remove this branch once SamRecordClipper has been updated
              rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), "Temporary Issue: Fix Clipping Issue")
              Nil
            case consensusLength =>
              // Pad out the vanilla consensus reads to the full length
              val paddedR1 = r1Consensus.padded(newLength = consensusLength, left = longestR1Alignment.negativeStrand, qual = 0)
              val paddedR2 = r2Consensus.padded(newLength = consensusLength, left = longestR2Alignment.negativeStrand, qual = 0)

              val consensus = duplexConsensus(Some(paddedR1), Some(paddedR2), sourceReads = None)
              if (longestR1Alignment.negativeStrand) consensus.foreach(_.revcomp())

              val umis = recs.iterator.map(r => r[String](ConsensusTags.UmiBases)).toSeq
              consensus
                .iterator
                .map(c => createSamRecord(c, ReadType.Fragment, umis))
                .toSeq
          }
        }
      }
    }
  }

  /**
    * Computes the length of the consensus read that should be generated. The provided positive and
    * negative strand alignments _must_ overlap.
    *
    * Returns -1 if the length couldn't be computed because the overlap of the reads ends on an indel.
    */
  private def computeConsensusLength(pos: SamRecord, neg: SamRecord): Int = {
    val refOverlapEnd = pos.end
    val posReadPos = pos.readPosAtRefPos(refOverlapEnd, returnLastBaseIfDeleted=false)
    val negReadPos = neg.readPosAtRefPos(refOverlapEnd, returnLastBaseIfDeleted=false)
    if (posReadPos == 0 || negReadPos == 0) -1 else {
      posReadPos + neg.length - negReadPos
    }
  }

  private def toSourceReadForCodec(rec: SamRecord): SourceRead = {
    // TODO: handle reads that map beyond the end of their mate
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    val cigar = if (rec.positiveStrand) rec.cigar else rec.cigar.reverse

    if (rec.negativeStrand) {
      Sequences.revcomp(bases)
      Sequences.reverse(quals)
    }

    SourceRead(
      id    = sourceMoleculeId(rec),
      bases = bases,
      quals = quals,
      cigar = cigar,
      sam   = Some(rec)
    )
  }
}
