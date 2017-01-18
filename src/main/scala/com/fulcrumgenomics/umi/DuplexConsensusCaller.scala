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

package com.fulcrumgenomics.umi

import java.lang.Math.min

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.umi.DuplexConsensusCaller._
import com.fulcrumgenomics.umi.UmiConsensusCaller.{SimpleRead, SourceRead}
import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType.{ReadType, _}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.SAMRecord

/**
  * Container for constant values and types used by the [[DuplexConsensusCaller]]
  */
object DuplexConsensusCaller {
  /** Various default values for the consensus caller. */
  val ErrorRatePreUmi: PhredScore         = 45.toByte
  val ErrorRatePostUmi: PhredScore        = 40.toByte
  val MinInputBaseQuality: PhredScore     = 15.toByte
  val NoCall: Byte = 'N'.toByte

  /**
    * Stores information about a consensus read.  Bases, arrays, and the two single
    * strand consensus reads must all be the same length.
    *
    * @param bases the base calls of the consensus read
    * @param quals the calculated phred-scaled quality scores of the bases
    * @param errors the total number of raw read errors vs. the final consensus sequence at each base
    * @param abConsensus the AB consensus read from which the duplex was constructed
    * @param baConsensus the AB consensus read from which the duplex was constructed
    */
  case class DuplexConsensusRead(id: String,
                                 bases: Array[Byte],
                                 quals: Array[Byte],
                                 errors: Array[Short],
                                 abConsensus: VanillaConsensusRead,
                                 baConsensus: VanillaConsensusRead) extends SimpleRead {
    require(bases.length == quals.length,  "Bases and qualities are not the same length.")
    require(bases.length == errors.length, "Bases and errors are not the same length.")
    require(bases.length == abConsensus.length, "Bases and AB consensus are not the same length.")
    require(bases.length == baConsensus.length, "Bases and BA consensus are not the same length.")
  }
}


/**
  * Creates duplex consensus reads from SAMRecords that have been grouped by their source molecule
  * but not yet by source strand.
  *
  * Filters incoming bases by quality before building the duplex.
  *
  * Output reads and bases are constructed only if there is at least one read from each source
  * molecule strand.  Otherwise no filtering is performed.
  *
  * Note that a consequence of the above is that the output reads can be shorter than _some_ of
  * the input reads if the input reads are of varying length; they will be the length at which
  * there is coverage from both source strands.
  *
  * @param readNamePrefix the prefix to apply to all consensus read names
  * @param readGroupId    the read group ID to apply to all created consensus reads
  * @param minInputBaseQuality the minimum input base quality score to use a raw read's base
  * @param errorRatePreUmi the estimated rate of errors in the DNA prior to attaching UMIs
  * @param errorRatePostUmi the estimated rate of errors in the DNA post attaching UMIs
  */
class DuplexConsensusCaller(override val readNamePrefix: String,
                            override val readGroupId: String    = "A",
                            val minInputBaseQuality: PhredScore = DuplexConsensusCaller.MinInputBaseQuality,
                            val errorRatePreUmi: PhredScore     = DuplexConsensusCaller.ErrorRatePreUmi,
                            val errorRatePostUmi: PhredScore    = DuplexConsensusCaller.ErrorRatePostUmi
                           ) extends UmiConsensusCaller[DuplexConsensusRead] {

  private val ssCaller = new VanillaUmiConsensusCaller(readNamePrefix="x", options=new VanillaUmiConsensusCallerOptions(
      errorRatePreUmi         = this.errorRatePreUmi,
      errorRatePostUmi        = this.errorRatePostUmi,
      minReads                = 1,
      minInputBaseQuality     = this.minInputBaseQuality,
      minConsensusBaseQuality = PhredScore.MinValue,
      producePerBaseTags      = true
    ))

  /**
    * Returns the MI tag minus the trailing suffix that identifies /A vs /B
    */
  override protected[umi] def sourceMoleculeId(rec: SAMRecord): String = {
    val mi = rec.getStringAttribute(ConsensusTags.MolecularId)
    require(mi != null, s"Read ${rec.getReadName} is missing it's ${ConsensusTags.MolecularId} tag.")
    val index = mi.lastIndexOf('/')
    require(index > 0, s"Read ${rec.getReadName}'s ${ConsensusTags} tag doesn't look like a duplex id: ${mi}")
    mi.substring(0, index)
  }

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SAMRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  override protected def consensusSamRecordsFromSamRecords(recs: Seq[SAMRecord]): Seq[SAMRecord] = {
    val (pairs, frags) = recs.partition(_.getReadPairedFlag)
    rejectRecords(frags)

    // Group the reads by /A vs. /B and ensure that /A is the first group and /B the second
    val groups = pairs.groupBy(_.getStringAttribute(ConsensusTags.MolecularId)).toSeq.sortBy { case (mi, rs) => mi }.map(_._2)

    groups match {
      case Seq() =>
        Nil
      case Seq(ab) =>
        rejectRecords(ab)
        Nil
      case Seq(ab, ba) =>
        // Fragments have no place in duplex land (and are filtered out previously anyway)!
        val (_, abR1s, abR2s) = subGroupRecords(ab)
        val (_, baR1s, baR2s) = subGroupRecords(ba)

        // Filter by common indel pattern with AB and BA together
        val filteredXs = filterToMostCommonAlignment(abR1s ++ baR2s)
        val filteredYs = filterToMostCommonAlignment(abR2s ++ baR1s)

        // Then split then back apart for SS calling and turn them into SourceReads
        val filteredAbR1s = filteredXs.filter(_.getFirstOfPairFlag).map(toSourceRead(_, this.minInputBaseQuality))
        val filteredAbR2s = filteredYs.filter(_.getSecondOfPairFlag).map(toSourceRead(_, this.minInputBaseQuality))
        val filteredBaR1s = filteredYs.filter(_.getFirstOfPairFlag).map(toSourceRead(_, this.minInputBaseQuality))
        val filteredBaR2s = filteredXs.filter(_.getSecondOfPairFlag).map(toSourceRead(_, this.minInputBaseQuality))

        // Call the single-stranded consensus reads
        val abR1Consensus = ssCaller.consensusCall(filteredAbR1s)
        val abR2Consensus = ssCaller.consensusCall(filteredAbR2s)
        val baR1Consensus = ssCaller.consensusCall(filteredBaR1s)
        val baR2Consensus = ssCaller.consensusCall(filteredBaR2s)

        // Call the duplex reads
        val duplexR1 = duplexConsensus(abR1Consensus, baR2Consensus, filteredAbR1s ++ filteredBaR2s)
        val duplexR2 = duplexConsensus(abR2Consensus, baR1Consensus, filteredAbR2s ++ filteredBaR1s)

        // Convert to SAMRecords and return
        (duplexR1, duplexR2) match {
          case (Some(r1), Some(r2)) =>
            Seq(createSamRecord(r1, FirstOfPair), createSamRecord(r2, SecondOfPair))
          case _                    =>
            rejectRecords(recs)
            Nil
        }
      case gs =>
        unreachable("SAMRecords supplied with more than two distinct MI values.")
    }
  }

  /**
    * Constructs a duplex consensus read from a pair of single strand consensus reads.
    * If either of the incoming reads are undefined, the duplex read will be undefined.
    *
    * @param ab the single-strand consensus from one strand
    * @param ba the single-strand consensus from the other strand
    * @return a duplex consensus if one can be built
    */
  private[umi] def duplexConsensus(ab: Option[VanillaConsensusRead],
                                   ba: Option[VanillaConsensusRead],
                                   sourceReads: Seq[SourceRead]): Option[DuplexConsensusRead] = {
    (ab, ba) match {
      case (Some(a), Some(b)) =>
        val len = min(a.length, b.length)
        val id  = a.id
        val bases  = new Array[Byte](len)
        val quals  = new Array[Byte](len)
        val errors = new Array[Short](len)

        forloop(from=0, until=len) { i =>
          val aBase = a.bases(i)
          val bBase = b.bases(i)
          val aQual = a.quals(i)
          val bQual = b.quals(i)

          val (base, qual) = {
            if      (aBase == NoCall || bBase == NoCall) (NoCall, PhredScore.MinValue)
            else if (aBase == bBase) (aBase, (aQual + bQual).toByte)
            else if (aQual > bQual)  (aBase, (aQual - bQual).toByte)
            else if (bQual > aQual)  (bBase, (bQual - aQual).toByte)
            else (NoCall, PhredScore.MinValue)
          }

          bases(i)  = base
          quals(i)  = qual
          errors(i) = min(sourceReads.count(s => s.bases(i) != NoCall && s.bases(i) != base), Short.MaxValue).toShort
        }

        Some(DuplexConsensusRead(id=id, bases, quals, errors, a.truncate(bases.length), b.truncate(bases.length)))
      case _ =>
        None
    }
  }

  /**
    * Creates a SAMRecord with a ton of additional tags annotating the duplex read.
    */
  override protected def createSamRecord(read: DuplexConsensusRead, readType: ReadType): SAMRecord = {
    val rec = super.createSamRecord(read, readType)
    val ab = read.abConsensus
    val ba = read.baConsensus

    // Calculate the total depths across both SS consensus reads
    val totalDepths = read.abConsensus.depths.zip(read.baConsensus.depths).map(x => x._1 + x._2)

    { import ConsensusTags.PerRead._
      rec.setAttribute(RawReadCount,       totalDepths.max)
      rec.setAttribute(MinRawReadCount,    totalDepths.min)
      rec.setAttribute(RawReadErrorRate,   sum(read.errors) / totalDepths.sum.toFloat)
      rec.setAttribute(AbRawReadCount,     read.abConsensus.depths.max.toInt)
      rec.setAttribute(BaRawReadCount,     read.baConsensus.depths.max.toInt)
      rec.setAttribute(AbMinRawReadCount,  read.abConsensus.depths.min.toInt)
      rec.setAttribute(BaMinRawReadCount,  read.baConsensus.depths.min.toInt)
      rec.setAttribute(AbRawReadErrorRate, sum(read.abConsensus.errors) / sum(read.abConsensus.depths).toFloat)
      rec.setAttribute(BaRawReadErrorRate, sum(read.baConsensus.errors) / sum(read.baConsensus.depths).toFloat)
    }

    { import ConsensusTags.PerBase._
      rec.setAttribute(AbRawReadCount,  read.abConsensus.depths)
      rec.setAttribute(BaRawReadCount,  read.baConsensus.depths)
      rec.setAttribute(AbRawReadErrors, read.abConsensus.errors)
      rec.setAttribute(BaRawReadErrors, read.baConsensus.errors)
    }

    rec
  }
}
