/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

import com.fulcrumgenomics.umi.ConsensusCallerOptions._
import com.fulcrumgenomics.util.LogDouble._
import com.fulcrumgenomics.util.PhredScore._
import com.fulcrumgenomics.util.{LogDouble, PhredScore, ProgressLogger}
import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools._
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object ConsensusCallerOptions {
  /** Various default values for the consensus caller. */
  val DefaultTag: String                             = "MI"
  val DefaultErrorRatePreUmi: PhredScore             = 45
  val DefaultErrorRatePostUmi: PhredScore            = 40
  val DefaultMaxBaseQuality: PhredScore              = 40
  val DefaultBaseQualityShift: PhredScore            = 10
  val DefaultMinConsensusBaseQuality: PhredScore     = 13
  val DefaultMinReads: Int                           = 1
  val DefaultMinMeanConsensusBaseQuality: PhredScore = 13
  val DefaultRequireConsensusForBothPairs: Boolean   = true
}

/** Holds the parameters/options for consensus calling. */
case class ConsensusCallerOptions(tag: String                             = DefaultTag,
                                  errorRatePreUmi: PhredScore             = DefaultErrorRatePreUmi,
                                  errorRatePostUmi: PhredScore            = DefaultErrorRatePostUmi,
                                  maxBaseQuality: PhredScore              = DefaultMaxBaseQuality,
                                  baseQualityShift: PhredScore            = DefaultBaseQualityShift,
                                  minConsensusBaseQuality: PhredScore     = DefaultMinConsensusBaseQuality,
                                  minReads: Int                           = DefaultMinReads,
                                  minMeanConsensusBaseQuality: PhredScore = DefaultMinMeanConsensusBaseQuality,
                                  requireConsensusForBothPairs: Boolean   = DefaultRequireConsensusForBothPairs
                                 )

/** Stores information about a consensus read. */
case class ConsensusRead(bases: String, quals: String)

object ConsensusCaller {
  private val DnaBasesUpperCase = Array('A', 'C', 'G', 'T')
  private val ThreeLogDouble = 3.0.toLogDouble
  private val TwoThirdsLogDouble = 2.0.toLogDouble / ThreeLogDouble

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  def consensusFromSamRecords(records: Seq[SAMRecord],
                              options: ConsensusCallerOptions = new ConsensusCallerOptions()
                              ): Option[ConsensusRead] = {
    val bases = records.map { rec =>
      if (rec.getReadNegativeStrandFlag) {
        val newBases = new Array[Byte](rec.getReadLength)
        Array.copy(rec.getReadBases, 0, newBases, 0, rec.getReadLength)
        SequenceUtil.reverseComplement(newBases)
        newBases
      }
      else rec.getReadBases
    }
    val quals = records.map { rec =>
      if (rec.getReadNegativeStrandFlag) rec.getBaseQualities.reverse
      else rec.getBaseQualities
    }.map(qs => qs.map(_.toDouble))
    consensusFromBasesAndQualities(bases=bases, quals=quals, options=options)
  }

  /** Creates a consensus read from the given read and quality string tuples.  If no consensus read was created, None is returned.
    * Currently used for testing. */
  private[umi] def consensusFromStringBasesAndQualities(basesAndQualities: Seq[(String, String)],
                                                        options: ConsensusCallerOptions = new ConsensusCallerOptions()
                                                       ): Option[ConsensusRead] = {
    consensusFromBasesAndQualities(
      bases = basesAndQualities.map(_._1.getBytes),
      quals = basesAndQualities.map { _._2 }.map(SAMUtils.fastqToPhred).map(qs => qs.map(_.toDouble)),
      options           = options
    )
  }

  /** Creates a consensus read from the given read and qualities sequences.  If no consensus read was created, None
    * is returned.
    *
    * The same number of base sequences and quality sequences should be given.
    * */
  private[umi] def consensusFromBasesAndQualities(bases: Seq[Array[Byte]],
                                                  quals: Seq[Array[PhredScore]],
                                                  options: ConsensusCallerOptions = new ConsensusCallerOptions()
                                                 ): Option[ConsensusRead] = {
    // TODO: uncomment me
    if (bases.length != quals.length) throw new IllegalArgumentException("The # of base and quality sequences must be the same.")
    if (bases.zip(quals).exists { case (b, q) => b.length != q.length }) {
      throw new IllegalArgumentException("Found a read with different lengths for its bases and qualities.")
    }

    // check to see if we have enough reads.
    if (bases.length < options.minReads) return None

    // extract the bases and qualities, and adjust the qualities based on the given options.s
    val qualSeqs: Seq[Seq[LogDouble]] = quals.map { qs =>
        adjustBaseQualities(
          quals            = qs.map(_.fromPhredScore),
          maxBaseQuality   = options.maxBaseQuality.fromPhredScore,
          baseQualityShift = options.baseQualityShift,
          errorRatePostUmi = options.errorRatePostUmi.fromPhredScore
        )
      } // NB: binary Phred, not ASCII

    // get the most likely consensus bases and qualities
    val (consensusBases, consensusQualities) = consensusCall(
      baseSeqs                = bases,
      qualSeqs                = qualSeqs,
      errorRatePreUmi         = options.errorRatePreUmi.fromPhredScore,
      minConsensusBaseQuality = options.minConsensusBaseQuality.fromPhredScore,
      minReads                = options.minReads
    )

    // check that the mean base quality is high enough
    if (LogDouble.mean(consensusQualities:_*).toPhredScore < options.minMeanConsensusBaseQuality) None
    else Some(ConsensusRead(bases=consensusBases, quals=consensusQualities.map(_.toPhredScoreChar).mkString))
  }

  /** Get the most likely consensus bases and qualities. Used for testing. */
  private[umi] def consensusCallFromStringBasesAndQualities(baseStrings: Seq[String],
                                 qualSeqs: Seq[Seq[LogDouble]], // probability of an error
                                 errorRatePreUmi: LogDouble         = DefaultErrorRatePreUmi.fromPhredScore,
                                 minConsensusBaseQuality: LogDouble = DefaultMinConsensusBaseQuality.fromPhredScore,
                                 minReads: Int                      = DefaultMinReads): (String, Seq[LogDouble]) = {
    consensusCall(
      baseSeqs=baseStrings.map(_.getBytes()),
      qualSeqs=qualSeqs,
      errorRatePreUmi=errorRatePreUmi,
      minConsensusBaseQuality=minConsensusBaseQuality,
      minReads=minReads
    )
  }

  /** Get the most likely consensus bases and qualities. */
  private[umi] def consensusCall(baseSeqs: Seq[Array[Byte]],
                                 qualSeqs: Seq[Seq[LogDouble]], // probability of an error
                                 errorRatePreUmi: LogDouble         = DefaultErrorRatePreUmi.fromPhredScore,
                                 minConsensusBaseQuality: LogDouble = DefaultMinConsensusBaseQuality.fromPhredScore,
                                 minReads: Int                      = DefaultMinReads): (String, Seq[LogDouble]) = {
    val maxReadLength = baseSeqs.map(_.length).max


    // Array
    // by cycle, by candidate base
    val likelihoods = Array.ofDim[LogDouble](maxReadLength, DnaBasesUpperCase.length)
    val numReads = new Array[Int](maxReadLength)
    // init
    for (baseIdx <- 0 until maxReadLength) {
      for (i <- DnaBasesUpperCase.indices) {
        likelihoods(baseIdx)(i) = OneProbability.fromPhredScore
      }
      numReads(baseIdx) = 0
    }

    // Calculate the likelihoods
    for (readIdx <- baseSeqs.indices) { // for each read
      val bases = baseSeqs(readIdx)
      for (baseIdx <- bases.indices) { // for each base in the read
        val base = bases(baseIdx)
        if (base != 'N') {
          val pError = qualSeqs(readIdx)(baseIdx)
          for (i <- DnaBasesUpperCase.indices) {
            val candidateBase = DnaBasesUpperCase(i)
            val likelihood = {
              if (base == candidateBase) pError.oneMinus() // 1.0 - Pr(Error)
              else pError / ThreeLogDouble //  Pr(Error) for this specific base, assuming the error distributes uniformly across the other three bases
            }
            likelihoods(baseIdx)(i) *= likelihood
          }
          numReads(baseIdx) += 1
        }
      }
    }

    // Calculate the posteriors
    val consensusBases = new Array[Char](maxReadLength)
    val consensusErrorProbabilities = new Array[LogDouble](maxReadLength)
    for (baseIdx <- likelihoods.indices) { // for each base in the read
      if (numReads(baseIdx) < minReads) {
        consensusBases(baseIdx) = SequenceUtil.N.toChar
        consensusErrorProbabilities(baseIdx) = OneProbability.fromPhredScore
      }
      else {
        // get the sum of the likelihoods
        // pick the base with the maximum posterior
        var likelihoodSum = ZeroProbability.fromPhredScore
        var maxPosterior = LogDouble.toLogDouble(0)
        var maxPosteriorIdx = -1
        for (i <- DnaBasesUpperCase.indices) {
          val likelihood = likelihoods(baseIdx)(i)
          likelihoodSum = likelihoodSum + likelihood
          if (maxPosterior < likelihood) {
            maxPosterior = likelihood
            maxPosteriorIdx = i
          }
        }
        // normalize, assumes a uniform prior, so omits from the following calculation
        maxPosterior /= likelihoodSum
        val pConsensusError = maxPosterior.oneMinus() // convert to probability of the called consensus being wrong

        // Masks a base if the phred score would be too low
        consensusBases(baseIdx) = {
          if (pConsensusError.toPhredScoreInt < minConsensusBaseQuality.toPhredScoreInt) SequenceUtil.N.toChar
          else DnaBasesUpperCase(maxPosteriorIdx)
        }

        // Factor in the pre-UMI error rate.
        // Pr(error) = Pr(any pre-UMI error AND correct consensus) + Pr(no pre-UMI error AND any error in consensus)
        //               + Pr(pre-UMI error AND error in consensus, that do not give us the correct bases)
        // The last term tries to capture the case where a pre-UMI error modifies the base (ex. A->C) but a sequencing
        // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
        val p = probabilityOfErrorTwoTrials(errorRatePreUmi, pConsensusError)

        // Cap the quality
        consensusErrorProbabilities(baseIdx) = Math.min(PhredScore.MaxValue, p.toPhredScoreInt).fromPhredScore
      }
    }
    (consensusBases.mkString, consensusErrorProbabilities)
  }

  /** Adjusts the given base qualities.  The base qualities are first shifted by `baseQualityShift`, then capped using
    * `maxBaseQuality`, and finally the `errorRatePostUmi` is incorporated.
    */
  private[umi] def adjustBaseQualities(quals: Seq[LogDouble],
                                       maxBaseQuality: LogDouble    = DefaultMaxBaseQuality.fromPhredScore,
                                       baseQualityShift: PhredScore = DefaultBaseQualityShift,
                                       errorRatePostUmi: LogDouble  = DefaultErrorRatePostUmi.fromPhredScore
                                      ): Seq[LogDouble] = {
    quals.map { qual =>
      // shift the base qualities, then cap it.
      val newQual = if (qual.toPhredScore < baseQualityShift) ZeroProbability.fromPhredScore
      else {
        val shiftedQual = (qual.toPhredScore - baseQualityShift).fromPhredScore
        if (maxBaseQuality.toPhredScore < shiftedQual.toPhredScore) maxBaseQuality
        else shiftedQual
      }
      // Pr(err) = Pr(no post-UMI error AND any sequencing error) + Pr(any post-UMI error and no sequencing error) +
      //             Pr(post-UMI error AND sequencing error)
      // The last term tries to capture the case where a post-UMI error modifies the base (ex. A->C) but a sequencing
      // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
      probabilityOfErrorTwoTrials(errorRatePostUmi, newQual)
    }
  }

  /** Computes the probability of seeing an error in the base sequence if there are two independent error processes.
    * We sum three terms:
    * 1. the probability of an error in trial one and no error in trial two: Pr(A=Error, B=NoError).
    * 2. the probability of no error in trial one and an error in trial two: Pr(A=NoError, B=Error).
    * 3. the probability of an error in both trials, but when the second trial does not reverse the error in first one, which
    *    for DNA (4 bases) would only occur 2/3 times: Pr(A=x->y, B=y->z) * Pr(x!=z | x!=y, y!=z, x,y,z \in {A,C,G,T})
    */
  private[umi] def probabilityOfErrorTwoTrials(prErrorTrialOne: LogDouble, prErrorTrialTwo: LogDouble): LogDouble = {
    val pr1 = prErrorTrialOne            * prErrorTrialTwo.oneMinus()
    val pr2 = prErrorTrialOne.oneMinus() * prErrorTrialTwo
    val pr3 = prErrorTrialOne            * prErrorTrialTwo            * TwoThirdsLogDouble
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
      while (this.iter.hasNext && this.iter.head.getStringAttribute(options.tag) == tagToMatch) {
        val rec = this.iter.next()
        buffer += rec
      }
      progress.map(_.record(buffer:_*))
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
    val curIdx = readIdx
    readIdx += 1
    val prefix = readNamePrefix.getOrElse(ConsensusCaller.longestCommonPrefix(records.map(_.getReadName)).getOrElse("CONSENSUS"))
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
    if (records.isEmpty) return
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
