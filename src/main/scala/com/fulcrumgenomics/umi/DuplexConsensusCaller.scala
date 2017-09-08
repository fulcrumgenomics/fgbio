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
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.umi.DuplexConsensusCaller._
import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType.{ReadType, _}
import com.fulcrumgenomics.umi.UmiConsensusCaller.{SimpleRead, SourceRead}
import com.fulcrumgenomics.util.NumericTypes.PhredScore

/**
  * Container for constant values and types used by the [[DuplexConsensusCaller]]
  */
object DuplexConsensusCaller {
  /** Various default values for the consensus caller. */
  val ErrorRatePreUmi: PhredScore         = 45.toByte
  val ErrorRatePostUmi: PhredScore        = 40.toByte
  val MinInputBaseQuality: PhredScore     = 15.toByte
  val NoCall: Byte = 'N'.toByte
  val NoCallQual   = PhredScore.MinValue

  /** Additional filter strings used when rejecting reads. */
  val FilterAorBOnly    = "Only Seeing A or B Reads"
  val FilterFragments   = "Being Fragment/Non-Paired Reads"
  val FilterSsConsensus = "Only Generating One Strand Consensus"
  val FilterCollision   = "Potential collision between independent duplex molecules"

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
  * Creates duplex consensus reads from SamRecords that have been grouped by their source molecule
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
  * @param trim if true, quality trim reads in addition to masking. If false just mask.
  * @param errorRatePreUmi the estimated rate of errors in the DNA prior to attaching UMIs
  * @param errorRatePostUmi the estimated rate of errors in the DNA post attaching UMIs
  */
class DuplexConsensusCaller(override val readNamePrefix: String,
                            override val readGroupId: String    = "A",
                            val minInputBaseQuality: PhredScore = DuplexConsensusCaller.MinInputBaseQuality,
                            val trim: Boolean = false,
                            val errorRatePreUmi: PhredScore     = DuplexConsensusCaller.ErrorRatePreUmi,
                            val errorRatePostUmi: PhredScore    = DuplexConsensusCaller.ErrorRatePostUmi
                           ) extends UmiConsensusCaller[DuplexConsensusRead] with LazyLogging {

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
  override protected[umi] def sourceMoleculeId(rec: SamRecord): String = {
    rec.get[String](ConsensusTags.MolecularId) match {
      case None =>
        throw new IllegalStateException(s"Read ${rec.name} is missing it's ${ConsensusTags.MolecularId} tag.")
      case Some(mi) =>
        val index = mi.lastIndexOf('/')
        require(index > 0, s"Read ${rec.name}'s ${ConsensusTags} tag doesn't look like a duplex id: ${mi}")
        mi.substring(0, index)
    }
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

    // Group the reads by /A vs. /B and ensure that /A is the first group and /B the second
    val groups = pairs.groupBy(r => r[String](ConsensusTags.MolecularId)).toSeq.sortBy { case (mi, _) => mi }.map(_._2)

    groups match {
      case Seq() =>
        Nil
      case Seq(ab) =>
        rejectRecords(ab, FilterAorBOnly)
        Nil
      case Seq(ab, ba) =>
        // Fragments have no place in duplex land (and are filtered out previously anyway)!
        val (_, abR1s, abR2s) = subGroupRecords(ab)
        val (_, baR1s, baR2s) = subGroupRecords(ba)

        // Get all the alignments to one end of the source molecule
        val singleStrand1 = abR1s ++ baR2s
        val singleStrand2 = abR2s ++ baR1s

        // The orientation of AB and BA reads should be:
        // AB R1: +  AB R2: -
        // BA R1: -  BA R2: +
        // or vice versa (AB-R1:-, AB-R2:+, AB-R1:-, AB-R2: +
        // Therefore, AB-R1s and BA-R2s should be on the same strand, and the same for AB-R2s and BA-R1s
        // Check for this explicitly here.
        if (singleStrand1.nonEmpty) {
          val ss1Flag = singleStrand1.head.negativeStrand
          if (singleStrand1.exists(_.negativeStrand != ss1Flag)) {
            val ss1MI   = singleStrand1.head.apply[String](ConsensusTags.MolecularId)
            rejectRecords(ab ++ ba, FilterCollision)
            logger.debug(s"Not all AB-R1s and BA-R2s were on the same strand for molecule with id: $ss1MI")
            return Nil
          }
        }
        if (singleStrand2.nonEmpty) {
          val ss2Flag = singleStrand2.head.negativeStrand
          if (singleStrand2.exists(_.negativeStrand != ss2Flag)) {
            val ss2MI   = singleStrand2.head.apply[String](ConsensusTags.MolecularId)
            rejectRecords(ab ++ ba, FilterCollision)
            logger.debug(s"Not all AB-R2s and BA-R1s were on the same strand for molecule with id: $ss2MI")
            return Nil
          }
        }

        // Filter by common indel pattern with AB and BA together
        val filteredXs = filterToMostCommonAlignment((abR1s ++ baR2s).flatMap(toSourceRead(_, this.minInputBaseQuality, this.trim)))
        val filteredYs = filterToMostCommonAlignment((abR2s ++ baR1s).flatMap(toSourceRead(_, this.minInputBaseQuality, this.trim)))

        // Then split then back apart for SS calling
        val filteredAbR1s = filteredXs.filter(_.sam.exists(_.firstOfPair))
        val filteredAbR2s = filteredYs.filter(_.sam.exists(_.secondOfPair))
        val filteredBaR1s = filteredYs.filter(_.sam.exists(_.firstOfPair))
        val filteredBaR2s = filteredXs.filter(_.sam.exists(_.secondOfPair))

        // Call the single-stranded consensus reads
        val abR1Consensus = ssCaller.consensusCall(filteredAbR1s)
        val abR2Consensus = ssCaller.consensusCall(filteredAbR2s)
        val baR1Consensus = ssCaller.consensusCall(filteredBaR1s)
        val baR2Consensus = ssCaller.consensusCall(filteredBaR2s)

        // Call the duplex reads
        val duplexR1Sources = filteredAbR1s ++ filteredBaR2s
        val duplexR2Sources = filteredAbR2s ++ filteredBaR1s
        val duplexR1 = duplexConsensus(abR1Consensus, baR2Consensus, duplexR1Sources)
        val duplexR2 = duplexConsensus(abR2Consensus, baR1Consensus, duplexR2Sources)

        // Convert to SamRecords and return
        (duplexR1, duplexR2) match {
          case (Some(r1), Some(r2)) =>
            Seq(
              createSamRecord(r1, FirstOfPair, toUmiBases(duplexR1Sources)),
              createSamRecord(r2, SecondOfPair, toUmiBases(duplexR2Sources))
            )
          case _                    =>
            // NB: some reads may have been rejected already in filterToMostCommonAlignment, so just
            //     reject those records that survived the initial filtering.
            val remainingRecs = filteredXs ++ filteredYs
            rejectRecords(remainingRecs.flatMap(_.sam), FilterSsConsensus)
            Nil
        }
      case _ =>
        unreachable("SamRecords supplied with more than two distinct MI values.")
    }
  }

  /** Extracts the UMI bases for each source read, ignoring reads without UMI bases */
  private def toUmiBases(reads: Seq[SourceRead]): Seq[String] = {
    // The source read may not have the RX tag, and if so, we ignore them through
    // the use of .flatMap and .get, with the latter returning an `Some(bases)` if the
    // UMI bases are present, `None` otherwise.
    reads.flatMap(_.sam).flatMap(_.get[String](ConsensusTags.UmiBases))
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

          // Capture the raw consensus base prior to masking it to N, so that we can compute
          // errors vs. the actually called base.
          val (rawBase, rawQual) = {
            if      (aBase == bBase) (aBase, (aQual + bQual).toByte)
            else if (aQual > bQual)  (aBase, (aQual - bQual).toByte)
            else if (bQual > aQual)  (bBase, (bQual - aQual).toByte)
            else                     (aBase, PhredScore.MinValue)
          }

          // Then mask it if appropriate
          val (base, qual) = if (aBase == NoCall || bBase == NoCall || rawQual == PhredScore.MinValue) (NoCall, NoCallQual) else (rawBase, rawQual)

          bases(i)  = base
          quals(i)  = qual

          errors(i) = min(sourceReads.count(s => s.length > i && isError(s.bases(i), rawBase)), Short.MaxValue).toShort
        }

        Some(DuplexConsensusRead(id=id, bases, quals, errors, a.truncate(bases.length), b.truncate(bases.length)))
      case _ =>
        None
    }
  }

  /** Function that returns true if the pair of bases are both valid/called bases and do not match each other. */
  @inline private def isError(lhs: Byte, rhs: Byte): Boolean = lhs != NoCall && rhs != NoCall && lhs != rhs

  /**
    * Creates a SamRecord with a ton of additional tags annotating the duplex read.
    */
  override protected def createSamRecord(read: DuplexConsensusRead, readType: ReadType, umis: Seq[String] = Seq.empty): SamRecord = {
    val rec = super.createSamRecord(read, readType, umis)
    val ab = read.abConsensus
    val ba = read.baConsensus

    // Calculate the total depths across both SS consensus reads
    val totalDepths = read.abConsensus.depths.zip(read.baConsensus.depths).map(x => x._1 + x._2)

    { import ConsensusTags.PerRead._
      rec(RawReadCount)       = totalDepths.max
      rec(MinRawReadCount)    = totalDepths.min
      rec(RawReadErrorRate)   = sum(read.errors) / totalDepths.sum.toFloat
      rec(AbRawReadCount)     = read.abConsensus.depths.max.toInt
      rec(BaRawReadCount)     = read.baConsensus.depths.max.toInt
      rec(AbMinRawReadCount)  = read.abConsensus.depths.min.toInt
      rec(BaMinRawReadCount)  = read.baConsensus.depths.min.toInt
      rec(AbRawReadErrorRate) = sum(read.abConsensus.errors) / sum(read.abConsensus.depths).toFloat
      rec(BaRawReadErrorRate) = sum(read.baConsensus.errors) / sum(read.baConsensus.depths).toFloat
    }

    { import ConsensusTags.PerBase._
      rec(AbRawReadCount)  = read.abConsensus.depths
      rec(BaRawReadCount)  = read.baConsensus.depths
      rec(AbRawReadErrors) = read.abConsensus.errors
      rec(BaRawReadErrors) = read.baConsensus.errors
    }

    rec
  }
}
