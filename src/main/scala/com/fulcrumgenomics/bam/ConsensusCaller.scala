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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.bam.ConsensusCallerOptions._
import com.fulcrumgenomics.util.{PhredScore, LogDouble, ProgressLogger}
import com.fulcrumgenomics.util.PhredScore._
import com.fulcrumgenomics.util.LogDouble._
import com.fulcrumgenomics.util.LogProbability._

import dagr.commons.CommonsDef._
import htsjdk.samtools._
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object ConsensusCallerOptions {
  /** Various default values for the consensus caller. */
  val DefaultTag: String                             = "RX"
  val DefaultErrorRatePreUmi: LogDouble              = 45.fromPhredScore
  val DefaultErrorRatePostUmi: LogDouble             = 40.fromPhredScore
  val DefaultMaxBaseQuality: LogDouble               = 40.fromPhredScore
  val DefaultBaseQualityShift: Double                = 10
  val DefaultMinConsensusBaseQuality: LogDouble      = 13.fromPhredScore
  val DefaultMinReads: Int                           = 1
  val DefaultMinMeanConsensusBaseQuality: LogDouble  = 13.fromPhredScore
  val DefaultRequireConsensusForBothPairs: Boolean   = true
}

/** Holds the parameters/options for consensus calling. */
case class ConsensusCallerOptions(tag: String                            = DefaultTag,
                                  errorRatePreUmi: LogDouble             = DefaultErrorRatePreUmi,
                                  errorRatePostUmi: LogDouble            = DefaultErrorRatePostUmi,
                                  maxBaseQuality: LogDouble              = DefaultMaxBaseQuality,
                                  baseQualityShift: Double               = DefaultBaseQualityShift,
                                  minConsensusBaseQuality: LogDouble     = DefaultMinConsensusBaseQuality,
                                  minReads: Int                          = DefaultMinReads,
                                  minMeanConsensusBaseQuality: LogDouble = DefaultMinMeanConsensusBaseQuality,
                                  requireConsensusForBothPairs: Boolean  = DefaultRequireConsensusForBothPairs
                                 )

/** Stores information about a consensus read. */
case class ConsensusRead(bases: String, quals: String)

object ConsensusCaller {
  private val DnaBasesUpperCase = Array('A', 'C', 'G', 'T')


  private val ThreeLogDouble = 3.0.toLogDouble
  private val TwoThirdsLogDouble = 2.0.toLogDouble / ThreeLogDouble

  /** Creates a consensus read from the given read string and qualities tuples.  If no consensus read was created, None is returned. */
  def consensusFromBasesAndQualities(basesAndQualities: Seq[(String, Seq[LogDouble])],
                                     options: ConsensusCallerOptions = new ConsensusCallerOptions()
                                     ): Option[ConsensusRead] = {
    if (basesAndQualities.exists { case (b, q) => b.length != q.length }) {
      throw new IllegalArgumentException("Found a read with the bases and qualities having different length")
    }
    // check to see if we have enough reads.
    if (basesAndQualities.length < options.minReads) return None

    // extract the bases and qualities, and adjust the qualities based on the given options.s
    val baseStrings: Seq[String] = basesAndQualities.map { _._1.toUpperCase }
    val qualSeqs: Seq[Seq[LogDouble]] = basesAndQualities.map(bq => bq._2)
      .map { quals =>
        adjustBaseQualities(
          quals            = quals,
          maxBaseQuality   = options.maxBaseQuality,
          baseQualityShift = options.baseQualityShift,
          errorRatePostUmi = options.errorRatePostUmi
        )
      } // NB: binary Phred, not ASCII

    // get the most likely consensus bases and qualities
    val (consensusBases, consensusQualities) = consensusCalls(
      baseStrings             = baseStrings,
      qualSeqs                = qualSeqs,
      errorRatePreUmi         = options.errorRatePreUmi,
      minConsensusBaseQuality = options.minConsensusBaseQuality,
      minReads                = options.minReads
    )

    // check that the mean base quality is high enough
    if (LogDouble.mean(consensusQualities:_*).toPhredScore < options.minMeanConsensusBaseQuality.toPhredScore) None
    else Some(ConsensusRead(bases=consensusBases, quals=consensusQualities.map(_.toPhredScoreChar).mkString))
  }

  /** Creates a consensus read from the given read and quality string tuples.  If no consensus read was created, None is returned.
    * Currently used for testing. */
  private[bam] def consensusFromStringBasesAndQualities(basesAndQualities: Seq[(String, String)],
                                     options: ConsensusCallerOptions = new ConsensusCallerOptions()
                                    ): Option[ConsensusRead] = {
    consensusFromBasesAndQualities(
      basesAndQualities = basesAndQualities.map { case (b, q) => (b, q.map(SAMUtils.fastqToPhred).map(fromPhredScore(_))) },
      options           = options
    )
  }

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  def consensusFromSamRecords(records: Seq[SAMRecord],
                              options: ConsensusCallerOptions = new ConsensusCallerOptions()
                              ): Option[ConsensusRead] = {
    val baseStrings = records.map { rec =>
      if (rec.getReadNegativeStrandFlag) SequenceUtil.reverseComplement(rec.getReadString.toUpperCase)
      else rec.getReadString.toUpperCase
    }
    val qualSeqs = records.map { rec =>
      //rec.getBaseQualities
      if (rec.getReadNegativeStrandFlag) rec.getBaseQualities.reverse.map(fromPhredScore(_)).toSeq
      else rec.getBaseQualities.map(fromPhredScore(_)).toSeq
    }
    consensusFromBasesAndQualities(basesAndQualities=baseStrings.zip(qualSeqs), options = options)
  }

  /** Get the most likely consensus bases and qualities. */
  private[bam] def consensusCalls(baseStrings: Seq[String],
                                  qualSeqs: Seq[Seq[LogDouble]], // probability of an error
                                  errorRatePreUmi: LogDouble         = DefaultErrorRatePreUmi,
                                  minConsensusBaseQuality: LogDouble = DefaultMinConsensusBaseQuality,
                                  minReads: Int                      = DefaultMinReads): (String, Seq[LogDouble]) = {
    val maxReadLength = baseStrings.map(_.length).max

    // go through each position in the read(s).
    val (consensusBases, consensusErrorProbabilities) = Seq.range(0, maxReadLength, 1).map { baseIdx =>
      val basesAndQualsAtIdx = baseStrings.zip(qualSeqs).filter { case (baseString, qualSeq) =>
          baseIdx < baseString.length && baseString(baseIdx) != 'N'
      }.map { case (baseString, qualSeq) =>
        (baseString(baseIdx), qualSeq(baseIdx))
      }

      // check to see if we heave a enough reads to produce a consensus base.
      if (basesAndQualsAtIdx.length < minReads) (SequenceUtil.N.toChar, OneProbability)
      else {
        // Get the likelihood of the data given each candidate consensus base
        val likelihoods = DnaBasesUpperCase.map { candidateBase =>
          // Get the likelihood of this single observation given the specific candidate consensus base
          basesAndQualsAtIdx.map { case (base, pError) =>
            if (base == candidateBase) pError.oneMinus() // 1.0 - Pr(Error)
            else pError / ThreeLogDouble //  Pr(Error) for this specific base, assuming the error distributes uniformly across the other three bases
          }.foldLeft(OneProbability)((likelihoodSum, likelihood) => likelihoodSum * likelihood)
        }
        // normalize, assumes a uniform prior, so omits from the above calculation
        val likelihoodSum = likelihoods.foldLeft(Zero)((a, b) => a + b)
        val posteriors = likelihoods.map { likelihood => likelihood / likelihoodSum }

        // Find the consensus base with the maximum posterior.  Since the probabilities are in log-space, still find the maximum.
        val maxPosterior    = posteriors.max
        val maxPosteriorIdx = posteriors.indexOf(maxPosterior)
        val pConsensusError = maxPosterior.oneMinus() // convert to probability of the called consensus being wrong
        // Masks a base if the phred score would be too low
        if (pConsensusError.toPhredScoreInt < minConsensusBaseQuality.toPhredScoreInt) (SequenceUtil.N.toChar, pConsensusError)
        else (DnaBasesUpperCase(maxPosteriorIdx), pConsensusError)
      }
    }.unzip

    // Factor in the pre-UMI error rate.
    val consensusErrorProbabilitiesScaled = consensusErrorProbabilities.map { pConsensusError =>
      // Pr(error) = Pr(any pre-UMI error AND correct consensus) + Pr(no pre-UMI error AND any error in consensus)
      //               + Pr(pre-UMI error AND error in consensus, that do not give us the correct bases)
      // The last term tries to capture the case where a pre-UMI error modifies the base (ex. A->C) but a sequencing
      // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
      val p = probabilityOfErrorTwoTrials(errorRatePreUmi, pConsensusError)
      // Cap the quality
      Math.min(PhredScore.MaxValue, p.toPhredScoreInt).fromPhredScore
    }
    (consensusBases.mkString, consensusErrorProbabilitiesScaled)
  }

  /** Adjusts the given base qualities.  The base qualities are first shifted by `baseQualityShift`, then capped using
    * `maxBaseQuality`, and finally the `errorRatePostUmi` is incorporated.
    */
  private[bam] def adjustBaseQualities(quals: Seq[LogDouble],
                                       maxBaseQuality: LogDouble   = DefaultMaxBaseQuality,
                                       baseQualityShift: Double = DefaultBaseQualityShift,
                                       errorRatePostUmi: LogDouble = DefaultErrorRatePostUmi
                                      ): Seq[LogDouble] = {
    quals.map { qual =>
      // shift the base qualities, then cap it.
      if (qual.toPhredScore < baseQualityShift) Zero
      else {
        val shiftedQual = (qual.toPhredScore - baseQualityShift).fromPhredScore
        if (maxBaseQuality.toPhredScore < shiftedQual.toPhredScore) maxBaseQuality
        else shiftedQual
      }
    }.map { qual =>
      // Pr(err) = Pr(no post-UMI error AND any sequencing error) + Pr(any post-UMI error and no sequencing error) +
      //             Pr(post-UMI error AND sequencing error)
      // The last term tries to capture the case where a post-UMI error modifies the base (ex. A->C) but a sequencing
      // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
      probabilityOfErrorTwoTrials(errorRatePostUmi, qual)
    }
  }

  /** Computes the probability of seeing an error in the base sequence if there are two independent error processes.
    * We sum three terms:
    * 1. the probability of an error in trial one and no error in trial two: Pr(A=Error, B=NoError).
    * 2. the probability of no error in trial one and an error in trial two: Pr(A=NoError, B=Error).
    * 3. the probability of an error in both trials, but when the second trial does not reverse the error in first one, which
    *    for DNA (4 bases) would only occur 2/3 times: Pr(A=x->y, B=y->z) * Pr(x!=z | x!=y, y!=z, x,y,z \in {A,C,G,T})
    */
  private[bam] def probabilityOfErrorTwoTrials(prErrorTrialOne: LogDouble, prErrorTrialTwo: LogDouble): LogDouble = {
    val pr1 = prErrorTrialOne       * prErrorTrialTwo.oneMinus()
    val pr2 = prErrorTrialOne.oneMinus() * prErrorTrialTwo
    val pr3 = prErrorTrialOne       * prErrorTrialTwo       * TwoThirdsLogDouble
    pr1 + pr2 + pr3
  }

  /** The type of consensus read to output. */
  private object ReadType extends Enumeration {
    val Fragment, FirstOfPair, SecondOfPair = Value
  }

  /** Gets the longest common prefix of the given strings, None if there is only one string or if there is an empty string. */
  @annotation.tailrec
  private def longestCommonPrefix(strs: Iterable[String], accu: Option[String] = None): Option[String] = {
    if (strs.exists(_.isEmpty) || strs.size <= 1) accu
    else {
      val first = strs.head.head
      if (strs.tail.exists(_.head != first)) accu
      else longestCommonPrefix(strs.map(_.tail), Some(accu.getOrElse("") + first))
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
class ConsensusCaller
( input: Iterator[SAMRecord],
  val header: SAMFileHeader,
  val readNamePrefix: Option[String]   = None,
  val readGroupId: String              = "A",
  val options: ConsensusCallerOptions  = new ConsensusCallerOptions(),
  val rejects: Option[SAMFileWriter]   = None,
  val progress: Option[ProgressLogger] = None
) extends Iterator[SAMRecord] {
  import ConsensusCaller.ReadType
  import ConsensusCaller.ReadType._
  import ConsensusCaller._

  private val iter = input.buffered
  private val nextConsensusRecords: mutable.Queue[SAMRecord] = mutable.Queue[SAMRecord]() // one per UMI group
  private var readIdx = 1

  /** True if there are more consensus reads, false otherwise. */
  def hasNext(): Boolean = this.nextConsensusRecords.nonEmpty || (this.iter.nonEmpty && advance())

  /** Returns the next consensus read. */
  def next(): SAMRecord = {
    if (!this.hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false.")
    this.nextConsensusRecords.dequeue()
  }

  /** Consumes records until a consensus read can be created, or no more input records are available. Returns
    * true if a consensus read was created, false otherwise. */
  @annotation.tailrec
  private def advance(): Boolean = {
    // get the records to create the consensus read
    val buffer = nextGroupOfRecords()

    // partition the records to which end of a pair it belongs, or if it is a fragment read.
    val (fragments, firstOfPair, secondOfPair) = subGroupRecords(records = buffer)

    // track if we are successful creating any consensus reads
    var success = false

    // fragment
    ConsensusCaller.consensusFromSamRecords(records=fragments, options=options) match {
      case None       => // reject
        rejectRecords(records=fragments);
      case Some(frag) => // output
        this.createAndEnqueueSamRecord(records=fragments, read=frag, readName=nextReadName(fragments), readType=Fragment)
        success = true
    }

    // pairs
    val needBothPairs = options.requireConsensusForBothPairs // for readability later
    val firstOfPairConsensus  = ConsensusCaller.consensusFromSamRecords(records=firstOfPair, options=options)
    val secondOfPairConsensus = ConsensusCaller.consensusFromSamRecords(records=secondOfPair, options=options)
    (firstOfPairConsensus, secondOfPairConsensus) match {
      case (None, None)                                       => // reject
        rejectRecords(records=firstOfPair ++ secondOfPair)
      case (Some(_), None) | (None, Some(_)) if needBothPairs => // reject
        rejectRecords(records=firstOfPair ++ secondOfPair)
      case (firstRead, secondRead)                            => // output
        this.createAndEnqueueSamRecordPair(firstRecords=firstOfPair, firstRead=firstRead, secondRecords=secondOfPair, secondRead=secondRead)
        success = true
    }

    if (success) true // consensus created
    else if (this.iter.isEmpty) false // no more records, don't try again
    else this.advance() // no consensus, but more records, try again
  }

  private def rejectRecords(records: Seq[SAMRecord]): Unit = this.rejects.foreach(rej => records.foreach(rej.addAlignment))

  /** Adds a SAM record from the underlying iterator to the buffer if either the buffer is empty or the SAM tag is
    * the same for the records in the buffer as the next record in the input iterator.  Returns true if a record was
    * added, false otherwise.
    */
  private def nextGroupOfRecords(): List[SAMRecord] = {
    if (this.iter.isEmpty) Nil
    else {
      val tagToMatch = this.iter.head.getStringAttribute(options.tag)
      val buffer = ListBuffer[SAMRecord]()
      while (this.iter.hasNext && this.iter.head.getStringAttribute(options.tag) == tagToMatch) buffer += this.iter.next
      buffer.toList
    }
  }

  /** Split records into those that should make a single-end consensus read, first of pair consensus read,
    * and second of pair consensus read, respectively.  The default method is to use the SAM flag to find
    * unpaired reads, first of pair reads, and second of pair reads.
    */
  protected def subGroupRecords(records: Seq[SAMRecord]): (Seq[SAMRecord], Seq[SAMRecord],Seq[SAMRecord]) = {
    val fragments    = records.filter { rec => !rec.getReadPairedFlag }
    val firstOfPair  = records.filter { rec => rec.getReadPairedFlag && rec.getFirstOfPairFlag }
    val secondOfPair = records.filter { rec => rec.getReadPairedFlag && rec.getSecondOfPairFlag }
    (fragments, firstOfPair, secondOfPair)
  }

  /** Returns the next read name with format "<prefix>:<idx>", where "<prefix>" is either the supplied prefix or the
    * longest common prefix of all read names, and "<idx>" is the 1-based consensus read index.  If no prefix was found,
    * "CONSENSUS" is used.  If no records are given, the empty string is returned.
    */
  private def nextReadName(records: Seq[SAMRecord]): String = {
    if (records.isEmpty) return ""
    val curIdx = yieldAndThen(readIdx)(readIdx += 1)
    val prefix = readNamePrefix.getOrElse(longestCommonPrefix(records.map(_.getReadName)).getOrElse("CONSENSUS"))
    s"$prefix:$curIdx"
  }

  /** Creates a `SAMRecord` for both ends of a pair.  If a consensus read is not given for one end of a pair, a dummy
    * record is created.  At least one consensus read must be given.
    */
  private def createAndEnqueueSamRecordPair(firstRecords: Seq[SAMRecord],
                                            firstRead: Option[ConsensusRead],
                                            secondRecords: Seq[SAMRecord],
                                            secondRead: Option[ConsensusRead]): Unit = {
    if (firstRead.isEmpty && secondRead.isEmpty) throw new IllegalArgumentException("Both consenus reads were empty.")
    val readName = nextReadName(firstRecords++secondRecords)
    // first end
    createAndEnqueueSamRecord(
      records  = firstRecords,
      read     = firstRead.getOrElse(dummyConsensusRead(secondRead.get)),
      readName = readName,
      readType = FirstOfPair
    )
    // second end
    createAndEnqueueSamRecord(
      records  = secondRecords,
      read     = secondRead.getOrElse(dummyConsensusRead(firstRead.get)),
      readName = readName,
      readType = SecondOfPair
    )
  }

  /** Creates a `SAMRecord` from the called consensus base and qualities. */
  private def createAndEnqueueSamRecord(records: Seq[SAMRecord],
                                        read: ConsensusRead,
                                        readName: String,
                                        readType: ReadType.Value): Unit = {
    val rec = new SAMRecord(header)
    rec.setReadName(readName)
    rec.setReadUnmappedFlag(true)
    readType match {
      case Fragment     => // do nothing
      case FirstOfPair  =>
        rec.setReadPairedFlag(true)
        rec.setFirstOfPairFlag(true)
        rec.setMateUnmappedFlag(true)
      case SecondOfPair =>
        rec.setReadPairedFlag(true)
        rec.setSecondOfPairFlag(true)
        rec.setMateUnmappedFlag(true)
    }
    rec.setReadString(read.bases)
    rec.setBaseQualityString(read.quals)
    rec.setAttribute(SAMTag.RG.name(), readGroupId)
    rec.setAttribute(options.tag, records.head.getStringAttribute(options.tag))
    // TODO: set custom SAM tags:
    // - # of reads contributing to this consensus

    // enqueue the record
    this.nextConsensusRecords.enqueue(rec)
  }

  /** Creates a dummy consensus read.  The read and quality strings will have the same length as the source, with
    * the read string being all Ns, and the quality string having zero base qualities. */
  private def dummyConsensusRead(source: ConsensusRead): ConsensusRead = {
    ConsensusRead(bases=source.bases.map(_ => 'N'), quals=source.quals.map(_ => PhredZeroChar))
  }
}
