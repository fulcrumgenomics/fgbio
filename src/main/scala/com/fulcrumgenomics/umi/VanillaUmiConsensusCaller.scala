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
 *
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef.forloop
import com.fulcrumgenomics.bam.api.{SamRecord, SamWriter}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType._
import com.fulcrumgenomics.umi.UmiConsensusCaller._
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.SAMTag

import java.util
import scala.util.Random

/**
  * Holds the defaults for consensus caller options.
  */
object VanillaUmiConsensusCallerOptions {
  /** Various default values for the consensus caller. */
  val DefaultTag: String                         = ConsensusTags.MolecularId
  val DefaultErrorRatePreUmi: PhredScore         = 45.toByte
  val DefaultErrorRatePostUmi: PhredScore        = 40.toByte
  val DefaultMinInputBaseQuality: PhredScore     = 10.toByte
  val DefaultMinConsensusBaseQuality: PhredScore = 40.toByte
  val DefaultMinReads: Int                       = 2
  val DefaultMaxReads: Int                       = Int.MaxValue
  val DefaultProducePerBaseTags: Boolean         = true
  val DefaultQualityTrim: Boolean                = false
}

/**
  * Holds the parameters/options for consensus calling.
  */
case class VanillaUmiConsensusCallerOptions
(
  tag: String                         = DefaultTag,
  errorRatePreUmi: PhredScore         = DefaultErrorRatePreUmi,
  errorRatePostUmi: PhredScore        = DefaultErrorRatePostUmi,
  minInputBaseQuality: PhredScore     = DefaultMinInputBaseQuality,
  qualityTrim: Boolean                = DefaultQualityTrim,
  minConsensusBaseQuality: PhredScore = DefaultMinConsensusBaseQuality,
  minReads: Int                       = DefaultMinReads,
  maxReads: Int                       = DefaultMaxReads,
  producePerBaseTags: Boolean         = DefaultProducePerBaseTags
)


/**
  * Stores information about a consensus read.  All four arrays are of equal length.
  *
  * Depths and errors that have values exceeding Short.MaxValue (32767) will be called
  * to Short.MaxValue.
  *
  * @param bases the base calls of the consensus read
  * @param quals the calculated phred-scaled quality scores of the bases
  * @param depths the number of raw reads that contributed to the consensus call at each position
  * @param errors the number of contributing raw reads that disagree with the final consensus base at each position
  * @param sourceReads optionally the source reads that went into calling this consensus
  */
case class VanillaConsensusRead(id: String, bases: Array[Byte], quals: Array[Byte], depths: Array[Short], errors: Array[Short], sourceReads: Option[Seq[SourceRead]] = None) extends SimpleRead {
  require(bases.length == quals.length,  "Bases and qualities are not the same length.")
  require(bases.length == depths.length, "Bases and depths are not the same length.")
  require(bases.length == errors.length, "Bases and errors are not the same length.")

  /** Truncates the read to the given length. If len > current length, the read is returned at current length. */
  def truncate(len: Int): VanillaConsensusRead = {
    if (len >= this.length) this
    else this.copy(bases=bases.take(len), quals=quals.take(len), depths=depths.take(len), errors=errors.take(len))
  }

  /**
    * Modifies all the bases, quals, depths and errors arrays *in place* to reverse complement the sequence
    * and related information in the consensus read instance.
    *
    * WARNING: modifies the record in place!
    */
  def revcomp(): Unit = {
    Sequences.revcomp(this.bases)
    Sequences.reverse(quals)
    Sequences.reverse(depths)
    Sequences.reverse(errors)
  }

  /**
    * Pads a consensus reads by adding bases to the either the left or right of the existing sequence.
    *
    * @param newLength the new total length of the consensus read - must be at least the existing length
    * @param left if true pad to the left of the existing sequence, else pad to the right of the existing sequence
    * @param base the base to use to fill in the padded sequence
    * @param qual the qual to use to fill in the qualities
    */
  def padded(newLength: Int, left: Boolean, base: Byte = 'n', qual: Byte = 2): VanillaConsensusRead = {
    require(newLength >= this.bases.length, "Cannot pad to a length shorter than current length.")
    if (newLength == this.bases.length) this else {
      val newBases  = new Array[Byte](newLength)
      val newQuals  = new Array[Byte](newLength)
      val newDepths = new Array[Short](newLength)
      val newErrors = new Array[Short](newLength)

      val oldLength = this.bases.length
      val addedLength = newLength - oldLength

      val destPos = if (left) addedLength else 0
      System.arraycopy(this.bases, 0,  newBases,  destPos, oldLength)
      System.arraycopy(this.quals, 0,  newQuals,  destPos, oldLength)
      System.arraycopy(this.depths, 0, newDepths, destPos, oldLength)
      System.arraycopy(this.errors, 0, newErrors, destPos, oldLength)

      val startIndex = if (left) 0 else oldLength
      util.Arrays.fill(newBases, startIndex, startIndex + addedLength, base)
      util.Arrays.fill(newQuals, startIndex, startIndex + addedLength, qual)

      VanillaConsensusRead(
        id          = this.id,
        bases       = newBases,
        quals       = newQuals,
        depths      = newDepths,
        errors      = newErrors,
        sourceReads = sourceReads,
      )
    }
  }
}

/** Calls consensus reads by grouping consecutive reads with the same SAM tag.
  *
  * Consecutive reads with the SAM tag are partitioned into fragments, first of pair, and
  * second of pair reads, and a consensus read is created for each partition.  A consensus read
  * for a given partition may not be returned if any of the conditions are not met (ex. minimum
  * number of reads, minimum mean consensus base quality, ...).
  * */
class VanillaUmiConsensusCaller(override val readNamePrefix: String,
                                override val readGroupId: String = "A",
                                val options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions(),
                                override val cellTag: Option[String] = Some(SAMTag.CB.name),
                                override val rejectsWriter: Option[SamWriter] = None,
                               ) extends UmiConsensusCaller[VanillaConsensusRead] with LazyLogging {

  initializeRejectCounts(_.usedByVanilla)

  private val NotEnoughReadsQual: PhredScore = 0.toByte // Score output when masking to N due to insufficient input reads
  private val TooLowQualityQual: PhredScore = 2.toByte  // Score output when masking to N due to too low consensus quality

  private val caller = new ConsensusCaller(errorRatePreLabeling  = options.errorRatePreUmi,
                                           errorRatePostLabeling = options.errorRatePostUmi)

  /** Map from input qual score to output qual score in the case where there is only one read going into the consensus. */
  private val SingleInputConsensusQuals: Array[Byte] = Range.inclusive(0, PhredScore.MaxValue.toInt).map { q =>
    val lnProbOne = LogProbability.fromPhredScore(q)
    val lnProbTwo = LogProbability.fromPhredScore(Math.min(this.options.errorRatePreUmi.toInt, this.options.errorRatePostUmi.toInt))
    PhredScore.fromLogProbability(LogProbability.probabilityOfErrorTwoTrials(lnProbOne, lnProbTwo))
  }.toArray

  private val random = new Random(42)

  /** Returns a clone of this consensus caller in a state where no previous reads were processed.  I.e. all counters
    * are set to zero.*/
  def emptyClone(): VanillaUmiConsensusCaller = {
    new VanillaUmiConsensusCaller(
      readNamePrefix = readNamePrefix,
      readGroupId    = readGroupId,
      options        = options,
      cellTag        = this.cellTag,
      rejectsWriter  = this.rejectsWriter
    )
  }

  /** Returns the value of the SAM tag directly. */
  override def sourceMoleculeId(rec: SamRecord): String = rec(this.options.tag)

  /** Takes in all the SamRecords for a single source molecule and produces consensus records. */
  override protected def consensusSamRecordsFromSamRecords(recs: Seq[SamRecord]): Seq[SamRecord] = {
    val cellBarcode: Option[String] = this.cellTag.flatMap { tag =>
      val barcodes = recs.flatMap(_.get[String](tag)).distinct
      require(barcodes.length <= 1, s"Multiple different cell barcodes found for tag $tag: $barcodes")
      barcodes.headOption
    }

    // partition the records to which end of a pair it belongs, or if it is a fragment read.
    val (fragments, firstOfPair, secondOfPair) = subGroupRecords(recs)
    val builder = IndexedSeq.newBuilder[SamRecord]

    // fragment
    consensusFromSamRecords(records=fragments).map { case frag =>
      builder += createSamRecord(
        read        = frag,
        readType    = Fragment,
        umis        = frag.sourceReads.getOrElse(Seq.empty).flatMap(rec => rec.sam.flatMap(_.get[String](ConsensusTags.UmiBases))),
        cellBarcode = cellBarcode,
      )
    }

    // pairs
    (consensusFromSamRecords(firstOfPair), consensusFromSamRecords(secondOfPair)) match {
      case (None, Some(_))      => rejectRecords(secondOfPair, RejectionReason.OrphanConsensus)
      case (Some(_), None)      => rejectRecords(firstOfPair,  RejectionReason.OrphanConsensus)
      case (None, None)         => rejectRecords(firstOfPair ++ secondOfPair, RejectionReason.OrphanConsensus)
      case (Some(read1), Some(read2)) =>
        builder += createSamRecord(
          read        = read1,
          readType    = FirstOfPair,
          umis        = read1.sourceReads.getOrElse(Seq.empty).flatMap(rec => rec.sam.flatMap(_.get[String](ConsensusTags.UmiBases))),
          cellBarcode = cellBarcode,
        )
        builder += createSamRecord(
          read        = read2,
          readType    = SecondOfPair,
          umis        = read2.sourceReads.getOrElse(Seq.empty).flatMap(rec => rec.sam.flatMap(_.get[String](ConsensusTags.UmiBases))),
          cellBarcode = cellBarcode,
        )
    }

    builder.result()
  }

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned.
   *
   * The source reads returned as part of the consensus read are the _filtered_ set of source reads that are used to
   * create the consensus call.  The read may be filtered for reasons including but not limited to: if the read is too
   * short after quality trimming, if the read does not share the most common alignment of the read sequence to the
   * reference, or if the number of reads are capped.
   * */
  protected[umi] def consensusFromSamRecords(records: Seq[SamRecord]): Option[VanillaConsensusRead] = {
    if (records.size < this.options.minReads) {
      rejectRecords(records, RejectionReason.InsufficientSupport)
      None
    }
    else {
      val sourceRecords   = records.flatMap(toSourceRead(_, this.options.minInputBaseQuality, this.options.qualityTrim))
      val filteredRecords = filterToMostCommonAlignment(sourceRecords)

      if (filteredRecords.size < records.size) {
        val r = records.head
        val n = if (r.paired && r.secondOfPair) "/2" else "/1"
        val m = r[String](this.options.tag)
        val discards = records.size - filteredRecords.size
        logger.debug("Discarded ", discards, "/", records.size, " records due to mismatched alignments for ", m, n)
      }

      if (filteredRecords.size >= this.options.minReads) {
        consensusCall(filteredRecords)
      } else {
        rejectRecords(filteredRecords.flatMap(_.sam), RejectionReason.InsufficientSupport)
        None
      }
    }
  }

  /** Creates a consensus read from the given read and qualities sequences.
    * If no consensus read was created, None is returned.
    *
    * The same number of base sequences and quality sequences should be given.
    * */
  private[umi] def consensusCall(reads: Seq[SourceRead]): Option[VanillaConsensusRead] = {
    // check to see if we have enough reads.
    if (reads.size < this.options.minReads) {
      None
    }
    else {
      // First limit to max reads if necessary
      val capped  = if (reads.size <= this.options.maxReads) reads else this.random.shuffle(reads).take(this.options.maxReads)
      // get the most likely consensus bases and qualities
      val consensusLength = consensusReadLength(capped, this.options.minReads)
      val consensusBases  = new Array[Base](consensusLength)
      val consensusQuals  = new Array[PhredScore](consensusLength)
      val consensusDepths = new Array[Short](consensusLength)
      val consensusErrors = new Array[Short](consensusLength)

      if (capped.length == 1) {
        val inBases = capped.head.bases
        val inQuals = capped.head.quals

        forloop (from=0, until=consensusLength) { i =>
          val rawBase      = inBases(i)
          val rawQual      = SingleInputConsensusQuals(inQuals(i).toInt)
          val (base, qual) = if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual) else (rawBase, rawQual)

          consensusBases(i)  = base
          consensusQuals(i)  = qual
          consensusDepths(i) = if (rawBase == NoCall) 0 else 1
          consensusErrors(i) = 0
        }
      }
      else {
        var positionInRead = 0
        val builder = this.caller.builder()
        while (positionInRead < consensusLength) {
          // Add the evidence from all reads that are long enough to cover this base
          capped.foreach { read =>
            if (read.length > positionInRead) {
              val base = read.bases(positionInRead)
              val qual = read.quals(positionInRead)
              if (base != NoCall) builder.add(base=base, qual=qual)
            }
          }

          val depth = builder.contributions // NB: cache this value, as it is re-computed each time

          // Call the consensus and do any additional filtering
          val (rawBase, rawQual) = builder.call()
          val (base, qual) = {
            if (depth < this.options.minReads)       (NoCall, NotEnoughReadsQual)
            else if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual)
            else (rawBase, rawQual)
          }

          consensusBases(positionInRead) = base
          consensusQuals(positionInRead) = qual

          // Generate the values for depth and count of errors
          val errors = if (rawBase == NoCall) depth else depth - builder.observations(rawBase)
          consensusDepths(positionInRead) = if (depth  > Short.MaxValue) Short.MaxValue else depth.toShort
          consensusErrors(positionInRead) = if (errors > Short.MaxValue) Short.MaxValue else errors.toShort

          // Get ready for the next pass
          builder.reset()
          positionInRead += 1
        }
      }

      Some(VanillaConsensusRead(
        id          = capped.head.id,
        bases       = consensusBases,
        quals       = consensusQuals,
        depths      = consensusDepths,
        errors      = consensusErrors,
        sourceReads = Some(capped)))
    }
  }

  /**
    * Calculates the length of the consensus read that should be produced. The length is calculated
    * as the maximum length at which minReads reads still have bases.
    *
    * @param reads the set of reads being fed into the consensus
    * @param minReads the minimum number of reads required
    * @return the length of consensus read that should be created
    */
  protected def consensusReadLength(reads: Seq[SourceRead], minReads: Int): Int = {
    require(reads.size >= minReads, "Too few reads to create a consensus.")
    reads.map(_.length).sortBy(len => -len).drop(minReads-1).head
  }

  /** Creates a `SamRecord` from the called consensus base and qualities. */
  override protected def createSamRecord(
    read: VanillaConsensusRead,
    readType: ReadType,
    umis: Seq[String]           = Seq.empty,
    cellBarcode: Option[String] = None,
  ): SamRecord = {
    val rec = super.createSamRecord(read, readType, umis, cellBarcode)
    // Set some additional information tags on the read
    rec(ConsensusTags.PerRead.RawReadCount)     = read.depths.max.toInt
    rec(ConsensusTags.PerRead.MinRawReadCount)  = read.depths.min.toInt
    rec(ConsensusTags.PerRead.RawReadErrorRate) = sum(read.errors) / sum(read.depths).toFloat
    if (this.options.producePerBaseTags) {
      rec(ConsensusTags.PerBase.RawReadCount)  = read.depths
      rec(ConsensusTags.PerBase.RawReadErrors) = read.errors
    }

    rec
  }
}
