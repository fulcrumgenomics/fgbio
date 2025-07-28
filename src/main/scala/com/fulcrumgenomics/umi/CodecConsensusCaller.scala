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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamWriter}
import com.fulcrumgenomics.bam.{ClippingMode, SamRecordClipper}
import com.fulcrumgenomics.umi.DuplexConsensusCaller.{DuplexConsensusRead, FilterFragments, FilterMinReads}
import com.fulcrumgenomics.umi.UmiConsensusCaller.{ReadType, SourceRead}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.Sequences

/**
  * Consensus caller for CODEC sequencing[1].  In CODEC each read-pair has an R1 which is generated from one
  * strand of the original duplex and an R2 which is generated from the opposite strand of the original duplex.
  * Therefore, even a single read-pair can generate a duplex consensus.  However, consensus yield is significantly
  * affected by how much the reads overlap as the overlapping region is the only place where information from
  * both strands is available.
  *
  * The result of consensus calling CODEC data is a single fragment read per duplex!
  *
  * The caller works approximately as follows (with checks at many points to see if sufficient data remains to
  * generate a duplex within the given parameters):
  *   - Read pairs are trimmed to remove any bases that extend past the beginning of their mates
  *   - Read 1s and Read 2s are _separately_ filtered for compatible CIGARs
  *   - Single strand consensus reads are formed from R1 and R2 separately
  *   - One of the SS reads is reverse complemented to bring both reads into the same orientation
  *   - The single strand reads are padded out with `n`s at Q0 if they do not fully overlap
  *   - A _single_ duplex consensus read is generated from the pair of single strand consensus reads
  *   - If the incoming R1s were on the negative strand of the genome, the duplex consensus is reverse complemented
  *
  * The CODEC consensus caller relies heavily on information about how the reads align to the genome in order to
  * synchronize the reads for duplex calling.  At this time read pairs whose alignments do not overlap on the genome,
  * including unmapped pairs, pairs with only one read mapped, chimeric pairs and read pairs whose insert size prevents
  * the reads from overlapping, are rejected by the CODEC caller.
  *
  * [1] https://doi.org/10.1038/s41588-023-01376-0
  *
  * @param readNamePrefix the prefix to apply to all consensus read names
  * @param readGroupId    the read group ID to apply to all created consensus reads
  * @param minInputBaseQuality the minimum base quality for the consensus caller to use a base in a raw read
  * @param errorRatePreUmi the estimated rate of errors in the DNA prior to attaching UMIs
  * @param errorRatePostUmi the estimated rate of errors in the DNA post attaching UMIs
  * @param minReadsPerStrand the minimum number of reads to form a single strand consensus from R1 or R2
  * @param maxReadsPerStrand the maximum number of reads to use when building a single strand consensus read before
  *                          triggering random downsampling
  * @param minDuplexLength the minimum length of the duplex region of a consensus read for the consensus read to be
  *                        built and emitted
  * @param singleStrandQual Reduce quality scores in single stranded regions of the consensus read to the given quality
  * @param outerBasesQual Reduce the first and last `outer-base-length` bases to the given quality
  * @param outerBasesLength The number of bases at the start and end of the read to reduce quality over *if*
  *                        outerBasesQual is specified
  * @param maxDuplexDisagreements filter out reads with more than this number of disagreements between the two single
  *                               strand consensus reads (only counting where both consensus reads make a call)
  * @param maxDuplexDisagreementRate filter out reads with more than this rate of disagreements between the two single
  *                                  strand consensus reads (only counting where both consensus reads make a call)
  */
class CodecConsensusCaller(readNamePrefix: String,
                           readGroupId: String = "A",
                           minInputBaseQuality: PhredScore = DuplexConsensusCaller.MinInputBaseQuality,
                           errorRatePreUmi: PhredScore = DuplexConsensusCaller.ErrorRatePreUmi,
                           errorRatePostUmi: PhredScore = DuplexConsensusCaller.ErrorRatePostUmi,
                           val minReadsPerStrand: Int = 1,
                           maxReadsPerStrand: Int = VanillaUmiConsensusCallerOptions.DefaultMaxReads,
                           minDuplexLength: Int = 1,
                           val singleStrandQual: Option[PhredScore] = None,
                           val outerBasesQual: Option[PhredScore] = None,
                           val outerBasesLength: Int = 5,
                           val maxDuplexDisagreements: Int = Int.MaxValue,
                           val maxDuplexDisagreementRate: Double = 1.0,
                           rejectsWriter: Option[SamWriter] = None
                          ) extends DuplexConsensusCaller(
  readNamePrefix      = readNamePrefix,
  readGroupId         = readGroupId,
  minInputBaseQuality = minInputBaseQuality,
  qualityTrim         = false,
  errorRatePreUmi     = errorRatePreUmi,
  errorRatePostUmi    = errorRatePostUmi,
  // We set minReads to all 1s here because the CODEC caller will do all the necessary filtering/limiting
  // ahead of calling the ss consensus calling process, and we absolutely don't want the SS caller to
  // trim off any regions that drop below the `minReadsPerStrand` threshold as that would significantly
  // complicate the logic of stitching R1 and R2 together.
  minReads            = Seq(1, 1, 1),
  maxReadsPerStrand   = maxReadsPerStrand,
  rejectsWriter       = rejectsWriter
) {
  require(this.minReadsPerStrand >= 1, "minReadsPerStrand must be at least 1.")

  private val FilterDuplexTooShort = s"Top/Bottom Strand Overlap < ${minDuplexLength} bases."
  private val FilterIndelError = "Indel error between strands in duplex overlap."
  private val FilterHighDuplexDisagreement = "Duplex had too many disagreements between strands."
  private val clipper = new SamRecordClipper(ClippingMode.Hard, autoClipAttributes = false)

  /** Returns the MI tag as is because in CODEC there is no /A or /B unlike classic duplex-seq */
  override protected[umi] def sourceMoleculeId(rec: SamRecord): String = rec[String](ConsensusTags.MolecularId)

  /** Returns a clone of this consensus caller in a state where no previous reads were processed.  I.e. all counters
    * are set to zero.  */
  override def emptyClone(): CodecConsensusCaller = {
    new CodecConsensusCaller(
      readNamePrefix      = this.readNamePrefix,
      readGroupId         = this.readGroupId,
      minInputBaseQuality = this.minInputBaseQuality,
      errorRatePreUmi     = this.errorRatePreUmi,
      errorRatePostUmi    = this.errorRatePostUmi,
      minReadsPerStrand   = this.minReadsPerStrand,
      maxReadsPerStrand   = this.maxReadsPerStrand,
      minDuplexLength     = this.minDuplexLength,
      singleStrandQual    = this.singleStrandQual,
      outerBasesQual      = this.outerBasesQual,
      outerBasesLength    = this.outerBasesLength,
      maxDuplexDisagreements    = this.maxDuplexDisagreements,
      maxDuplexDisagreementRate = this.maxDuplexDisagreementRate,
      rejectsWriter             = this.rejectsWriter
    )
  }

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
        val longestR1Alignment = r1s.iterator.flatMap(_.sam).maxBy(_.cigar.lengthOnTarget)
        val longestR2Alignment = r2s.iterator.flatMap(_.sam).maxBy(_.cigar.lengthOnTarget)

        val (longestPosAln, longestNegAln) =
          if (longestR1Alignment.positiveStrand) (longestR1Alignment, longestR2Alignment)
          else (longestR2Alignment, longestR1Alignment)

        // Calculate the overlapping region in reference space; this calculation only works because we
        // can assert from checks above that we have an FR pair that sequences towards each other.
        val (overlapStart, overlapEnd) = (longestNegAln.start, longestPosAln.end)
        val overlapLength = overlapEnd - overlapStart + 1

        // If the overlap isn't long enough then reject the records and return no consensus reads
        if (overlapLength < minDuplexLength) {
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

              // Calculate number and rate of duplex errors ... but don't really if the thresholds won't need them
              val (duplexErrors, duplexErrorRate) = {
                if (this.maxDuplexDisagreementRate >= 1 && this.maxDuplexDisagreements >= overlapLength) (0, 0.0)
                else {
                  val (duplexBases, duplexErrors) = computeDuplexCounts(paddedR1, paddedR2)
                  val duplexErrorRate = if (duplexBases > 0) duplexErrors.toDouble / duplexBases.toDouble else 0
                  (duplexErrors, duplexErrorRate)
                }
              }

              if (duplexErrors > this.maxDuplexDisagreements || duplexErrorRate > this.maxDuplexDisagreementRate) {
                rejectRecords((r1s.view ++ r2s.view).flatMap(_.sam), FilterHighDuplexDisagreement)
                Nil
              }
              else {
                val consensus = duplexConsensus(Some(paddedR1), Some(paddedR2), sourceReads = None)
                if (longestR1Alignment.negativeStrand) consensus.foreach(_.revcomp())
                consensus.foreach(maskCodecConsensusQuals)

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

  /** Applies quality masking to the consensus read based on the singleStrandQual and outerBasesQual parameters. */
  private[umi] def maskCodecConsensusQuals(rec: DuplexConsensusRead): Unit = {
    val quals   = rec.quals

    // Mask the ends of the reads if desired
    if (this.outerBasesLength > 0) {
      val lastIdx = quals.length - 1
      this.outerBasesQual.foreach { q =>
        forloop (from=0, until=Math.min(this.outerBasesLength, quals.length)) { i =>
          quals(i) = q
          quals(lastIdx - i) = q
        }
      }
    }

    // Mask single-stranded regions of the consensus
    this.singleStrandQual.foreach { q =>
      val a = rec.abConsensus.bases
      val b = rec.baConsensus.getOrElse(unreachable("Codec consensus found without a second strand consensus.")).bases

      forloop (from=0, until=quals.length) { i =>
        if (a(i) == 'n' || a(i) == 'N' || b(i) == 'N' || b(i) == 'n') {
          quals(i) = q
        }
      }
    }
  }

  /**
    * Calculates and returns two counts.  The first is the number of bases that have non-N coverage in both the
    * a and b consensus reads.  The second is the subset of those bases that do not agree between the a and b reads.
    *
    * @param a one of the two single stranded consensus reads
    * @param b the other single strand consensus reads
    * @return a tuple of (duplexBases, duplexDisagreements)
    */
  private[umi] def computeDuplexCounts(a: VanillaConsensusRead, b: VanillaConsensusRead): (Int, Int) = {
    var duplexBases  = 0
    var duplexErrors = 0

    forloop (from=0, until=a.length) { i =>
      val aBase = a.bases(i)
      val bBase = b.bases(i)

      if (aBase != 'n' && aBase != 'N' && bBase != 'n' && bBase != 'N') {
        duplexBases += 1

        if (aBase != bBase) duplexErrors += 1
      }
    }

    (duplexBases, duplexErrors)
  }
}
