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

import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType._
import com.fulcrumgenomics.umi.UmiConsensusCaller._
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import com.fulcrumgenomics.commons.util.LazyLogging
import htsjdk.samtools._

import scala.collection.mutable.ListBuffer
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
  */
case class VanillaConsensusRead(id: String, bases: Array[Byte], quals: Array[Byte], depths: Array[Short], errors: Array[Short]) extends SimpleRead {
  require(bases.length == quals.length,  "Bases and qualities are not the same length.")
  require(bases.length == depths.length, "Bases and depths are not the same length.")
  require(bases.length == errors.length, "Bases and errors are not the same length.")

  /** Truncates the read to the given length. If len > current length, the read is returned at current length. */
  def truncate(len: Int): VanillaConsensusRead = {
    if (len >= this.length) this
    else this.copy(bases=bases.take(len), quals=quals.take(len), depths=depths.take(len), errors=errors.take(len))
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
                                val rejects: Option[SAMFileWriter] = None
                               ) extends UmiConsensusCaller[VanillaConsensusRead] with LazyLogging {

  private val NotEnoughReadsQual: PhredScore = 0.toByte // Score output when masking to N due to insufficient input reads
  private val TooLowQualityQual: PhredScore = 2.toByte  // Score output when masking to N due to too low consensus quality
  private val DnaBasesUpperCase: Array[Byte] = Array('A', 'C', 'G', 'T').map(_.toByte)
  private val LogThree = LogProbability.toLogProbability(3.0)

  private val caller = new ConsensusCaller(errorRatePreLabeling  = options.errorRatePreUmi,
                                           errorRatePostLabeling = options.errorRatePostUmi)

  private val random = new Random(42)

  /** Returns the value of the SAM tag directly. */
  override def sourceMoleculeId(rec: SAMRecord): String = rec.getStringAttribute(this.options.tag)

  /** Takes in all the SAMRecords for a single source molecule and produces consensus records. */
  override protected def consensusSamRecordsFromSamRecords(recs: Seq[SAMRecord]): Seq[SAMRecord] = {
    // partition the records to which end of a pair it belongs, or if it is a fragment read.
    val (fragments, firstOfPair, secondOfPair) = subGroupRecords(recs)
    val buffer = new ListBuffer[SAMRecord]()

    // fragment
    consensusFromSamRecords(records=fragments).map { frag => buffer += createSamRecord(read=frag, readType=Fragment) }

    // pairs
    (consensusFromSamRecords(firstOfPair), consensusFromSamRecords(secondOfPair)) match {
      case (None, Some(r2))     => rejectRecords(secondOfPair, UmiConsensusCaller.FilterOrphan)
      case (Some(r1), None)     => rejectRecords(firstOfPair,  UmiConsensusCaller.FilterOrphan)
      case (None, None)         => rejectRecords(firstOfPair ++ secondOfPair, UmiConsensusCaller.FilterOrphan)
      case (Some(r1), Some(r2)) =>
        buffer += createSamRecord(r1, FirstOfPair)
        buffer += createSamRecord(r2, SecondOfPair)
    }

    buffer
  }

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  protected[umi] def consensusFromSamRecords(records: Seq[SAMRecord]): Option[VanillaConsensusRead] = {
    if (records.size < this.options.minReads) {
      rejectRecords(records, UmiConsensusCaller.FilterInsufficientSupport)
      None
    }
    else {
      val sourceRecords   = records.flatMap(toSourceRead(_, this.options.minInputBaseQuality, this.options.qualityTrim))
      val filteredRecords = filterToMostCommonAlignment(sourceRecords)

      if (filteredRecords.size < records.size) {
        val r = records.head
        val n = if (r.getReadPairedFlag && r.getSecondOfPairFlag) "/2" else "/1"
        val m = r.getAttribute(this.options.tag)
        val discards = records.size - filteredRecords.size
        logger.debug("Discarded ", discards, "/", records.size, " records due to mismatched alignments for ", m, n)
      }

      if (filteredRecords.size >= this.options.minReads) consensusCall(filteredRecords) else {
        rejectRecords(filteredRecords.flatMap(_.sam), UmiConsensusCaller.FilterInsufficientSupport)
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

        // Call the consensus and do any additional filtering
        val (rawBase, rawQual) = builder.call()
        val (base, qual) = {
          if (builder.contributions < this.options.minReads)       (NoCall, NotEnoughReadsQual)
          else if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual)
          else (rawBase, rawQual)
        }

        consensusBases(positionInRead) = base
        consensusQuals(positionInRead) = qual

        // Generate the values for depth and count of errors
        val depth  = builder.contributions
        val errors = if (rawBase == NoCall) depth else depth - builder.observations(rawBase)
        consensusDepths(positionInRead) = if (depth  > Short.MaxValue) Short.MaxValue else depth.toShort
        consensusErrors(positionInRead) = if (errors > Short.MaxValue) Short.MaxValue else errors.toShort

        // Get ready for the next pass
        builder.reset()
        positionInRead += 1
      }

      Some(VanillaConsensusRead(id=capped.head.id, bases=consensusBases, quals=consensusQuals, depths=consensusDepths, errors=consensusErrors))
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

  /** If a reject writer was provided, emit the reads to that writer. */
  override protected def rejectRecords(recs: Traversable[SAMRecord], reason: String): Unit = {
    super.rejectRecords(recs, reason)
    this.rejects.foreach(rej => recs.foreach(rej.addAlignment))
  }

  /** Creates a `SAMRecord` from the called consensus base and qualities. */
  override protected def createSamRecord(read: VanillaConsensusRead, readType: ReadType): SAMRecord = {
    val rec = super.createSamRecord(read, readType)
    // Set some additional information tags on the read
    rec.setAttribute(ConsensusTags.PerRead.RawReadCount,     read.depths.max.toInt)
    rec.setAttribute(ConsensusTags.PerRead.MinRawReadCount,  read.depths.min.toInt)
    rec.setAttribute(ConsensusTags.PerRead.RawReadErrorRate, sum(read.errors) / sum(read.depths).toFloat)
    if (this.options.producePerBaseTags) {
      rec.setAttribute(ConsensusTags.PerBase.RawReadCount,  read.depths)
      rec.setAttribute(ConsensusTags.PerBase.RawReadErrors, read.errors)
    }

    rec
  }
}
