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
import com.fulcrumgenomics.util.NumericTypes._
import com.fulcrumgenomics.util.{MathUtil, ProgressLogger}
import htsjdk.samtools._
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object ConsensusCallerOptions {
  type PhredScore = Byte

  /** Various default values for the consensus caller. */
  val DefaultTag: String                             = "MI"
  val DefaultErrorRatePreUmi: PhredScore             = 45.toByte
  val DefaultErrorRatePostUmi: PhredScore            = 40.toByte
  val DefaultMaxBaseQuality: PhredScore              = 40.toByte
  val DefaultBaseQualityShift: PhredScore            = 10.toByte
  val DefaultMinConsensusBaseQuality: PhredScore     = 13.toByte
  val DefaultMinReads: Int                           = 1
  val DefaultMinMeanConsensusBaseQuality: PhredScore = 13.toByte
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
                                 ) {

  val errorRatePreUmiLn  = LogProbability.fromPhredScore(errorRatePreUmi)
  val errorRatePostUmiLn = LogProbability.fromPhredScore(errorRatePostUmi)
}

/** Stores all the information about a read going into a consensus. */
case class AdjustedRead(bases: Array[Byte], pError: Array[LogProbability], pCorrect: Array[LogProbability]) {
  assert(bases.length == pError.length,   "Bases and qualities not the same length.")
  assert(bases.length == pCorrect.length, "Bases and qualities not the same length.")

  def length: Int = bases.length
}

/** Stores information about a read to be fed into a consensus. */
case class SourceRead(bases: Array[Byte], quals: Array[Byte]) {
  assert(bases.length == quals.length,   "Bases and qualities not the same length.")
  def length = bases.length
}

/** Stores information about a consensus read. */
case class ConsensusRead(bases: Array[Byte], quals: Array[Byte]) {
  val meanQuality: Byte = MathUtil.mean(quals)

  /** Returns the consensus read a String - mostly useful for testing. */
  def baseString = new String(bases)
}

/** Contains methods to help perform perform various consensus calling operations efficiently. Requires some upfront
  * computation of various probabilities. */
trait ConsensusCallerMath {

  /** The Phred-scaled error rate for an error post the UMIs have been integrated. */
  protected def errorRatePostUmi: PhredScore

  // Some useful constants
  protected val LogThree      = LogProbability.toLogProbability(3.0)
  protected val LogFourThirds = LogProbability.toLogProbability(4.0) - LogThree

  /** Computes the probability of seeing an error in the base sequence if there are two independent error processes.
    * We sum three terms:
    * 1. the probability of an error in trial one and no error in trial two: Pr(A=Error, B=NoError).
    * 2. the probability of no error in trial one and an error in trial two: Pr(A=NoError, B=Error).
    * 3. the probability of an error in both trials, but when the second trial does not reverse the error in first one, which
    *    for DNA (4 bases) would only occur 2/3 times: Pr(A=x->y, B=y->z) * Pr(x!=z | x!=y, y!=z, x,y,z \in {A,C,G,T})
    */
  def probabilityOfErrorTwoTrials(prErrorTrialOne: LogProbability, prErrorTrialTwo: LogProbability): LogProbability = {
    if (prErrorTrialOne < prErrorTrialTwo) probabilityOfErrorTwoTrials(prErrorTrialTwo, prErrorTrialOne)
    else if (prErrorTrialOne - prErrorTrialTwo >= 6) prErrorTrialOne // a simple approximation since prErrorTrialOne will dominate
    else {
      // f(X, Y) = X(1-Y) + (1-X)Y + 2/3*XY
      //         = X - XY + Y - XY + 2/3*XY
      //         = X + Y - 2XY + 2/3*XY
      //         = X + Y + XY*(2/3 - 6/3)
      //         = X + Y - 4/3*XY
      val term1 = LogProbability.or(prErrorTrialOne, prErrorTrialTwo) // X + Y
      val term2 = LogFourThirds + prErrorTrialOne + prErrorTrialTwo // 4/3*XY
      LogProbability.aOrNotB(term1, term2)
    }
  }

  /** Pre-computes the the log-scale probabilities of an error for each a phred-scaled base quality from 0-100. */
  protected lazy val phredToLogProbError: Array[Double] = Range(0, 100).toArray.map(p => {
    val e1 = LogProbability.fromPhredScore(errorRatePostUmi)
    val e2 = LogProbability.fromPhredScore(p.toByte)
    probabilityOfErrorTwoTrials(e1, e2)
  })

  /** Pre-computes the the log-scale probabilities of an not an error for each a phred-scaled base quality from 0-100. */
  protected lazy val phredToLogProbCorrect: Array[Double] = phredToLogProbError.map(LogProbability.not)

  /**
    * Adjusts the given base qualities.  The base qualities are first shifted by `baseQualityShift`, then capped using
    *`maxBaseQuality`, and finally the `errorRatePostUmi` is incorporated.
    *
    * Implemented as a while loop with assignment into pre-created arrays as this function is called on every
    * single input read and unfortunately Array[Double].map() causes boxing on all the values which is a
    * significant performance drain.
    */
  def adjustBaseQualities(read: SourceRead, maxBaseQuality: PhredScore, baseQualityShift: PhredScore): AdjustedRead = {
    val len = read.length
    val pErrors   = new Array[LogProbability](len)
    val pCorrects = new Array[LogProbability](len)

    val qs = read.quals
    var i = 0
    while (i < len) {
      val q = qs(i)
      val newQ = Math.min(maxBaseQuality, Math.max(PhredScore.MinValue, q - baseQualityShift))
      pErrors(i)   = this.phredToLogProbError(newQ)
      pCorrects(i) = this.phredToLogProbCorrect(newQ)
      i += 1
    }

    new AdjustedRead(read.bases, pErrors, pCorrects)
  }

  /**
    * Adjusts the given base qualities.  The base qualities are first shifted by `baseQualityShift`, then capped using
    *`maxBaseQuality`, and finally the `errorRatePostUmi` is incorporated.
    *
    * Implemented as a while loop with assignment into pre-created arrays as this function is called on every
    * single input read and unfortunately Array[Double].map() causes boxing on all the values which is a
    * significant performance drain.
    */
  def adjustBaseQualities(read: SourceRead, options: ConsensusCallerOptions): AdjustedRead = {
    adjustBaseQualities(read=read, maxBaseQuality=options.maxBaseQuality, baseQualityShift=options.baseQualityShift)
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
) extends Iterator[SAMRecord] with ConsensusCallerMath {
  /** The type of consensus read to output. */
  private object ReadType extends Enumeration {
    val Fragment, FirstOfPair, SecondOfPair = Value
  }
  import ReadType._

  private val DnaBasesUpperCase: Array[Byte] = Array('A', 'C', 'G', 'T').map(_.toByte)

  private val iter = input.buffered
  private val nextConsensusRecords: mutable.Queue[SAMRecord] = mutable.Queue[SAMRecord]() // one per UMI group
  private var readIdx = 1

  protected val errorRatePostUmi = this.options.errorRatePostUmi

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  def consensusFromSamRecords(records: Seq[SAMRecord]): Option[ConsensusRead] = {
    val sourceReads = records.map { rec =>
      if (rec.getReadNegativeStrandFlag) {
        val newBases = rec.getReadBases.clone()
        val newQuals = rec.getBaseQualities.clone()
        SequenceUtil.reverseComplement(newBases)
        SequenceUtil.reverse(newQuals, 0, newQuals.length)
        SourceRead(newBases, newQuals)
      }
      else {
        SourceRead(rec.getReadBases.array, rec.getBaseQualities.array)
      }
    }

    consensusCall(sourceReads)
  }

  /** Creates a consensus read from the given read and qualities sequences.  If no consensus read was created, None
    * is returned.
    *
    * The same number of base sequences and quality sequences should be given.
    * */
  private[umi] def consensusCall(reads: Seq[SourceRead]): Option[ConsensusRead] = {
    // check to see if we have enough reads.
    if (reads.length < options.minReads) return None

    // extract the bases and qualities, and adjust the qualities based on the given options.s
    val adjustedReads = reads.map(read => adjustBaseQualities(read, options))

    // get the most likely consensus bases and qualities
    val consensusRead = consensusCallAdjustedReads(reads=adjustedReads)

    // check that the mean base quality is high enough
    if (consensusRead.meanQuality < options.minMeanConsensusBaseQuality) None
    else Some(consensusRead)
  }

  /** Get the most likely consensus bases and qualities. */
  private[umi] def consensusCallAdjustedReads(reads: Seq[AdjustedRead]): ConsensusRead = {
    type Base = Byte
    val maxReadLength = reads.map(_.length).max

    // Array
    // by cycle, by candidate base
    val likelihoods = Array.ofDim[LogProbability](maxReadLength, DnaBasesUpperCase.length) // NB: array is initialized to zero, which is toLogProbability(1.0)
    val numReads = new Array[Int](maxReadLength)

    // Calculate the likelihoods
    {
      val readCount = reads.length
      val DnaBasesUpperCaseCount = DnaBasesUpperCase.length

      var readIdx = 0
      while (readIdx < readCount) {
        // for each read
        val read = reads(readIdx)
        val bases = read.bases
        val baseCount = bases.length

        var baseIdx = 0
        while (baseIdx < baseCount) {
          // for each base in the read
          val base = bases(baseIdx)
          if (base != 'N') {

            var i = 0
            while (i < DnaBasesUpperCaseCount) {
              val candidateBase = DnaBasesUpperCase(i)
              val likelihood = {
                if (base == candidateBase) read.pCorrect(baseIdx)
                else LogProbability.normalizeByScalar(read.pError(baseIdx), 3) //  Pr(Error) for this specific base, assuming the error distributes uniformly across the other three bases
              }
              likelihoods(baseIdx)(i) = LogProbability.and( likelihoods(baseIdx)(i), likelihood)
              i += 1
            }
            numReads(baseIdx) += 1
          }

          baseIdx +=1
        }
        readIdx += 1
      }
    }

    // Calculate the posteriors
    val consensusBases     = new Array[Base](maxReadLength)
    val consensusQualities = new Array[PhredScore](maxReadLength)
    for (baseIdx <- likelihoods.indices) { // for each base in the read
      if (numReads(baseIdx) < this.options.minReads) {
        consensusBases(baseIdx) = SequenceUtil.N
        consensusQualities(baseIdx) = PhredScore.MinValue
      }
      else {
        // get the sum of the likelihoods
        // pick the base with the maximum posterior
        val lls  = likelihoods(baseIdx)
        val likelihoodSum   = LogProbability.or(lls)
        val (maxLikelihood, maxLlIndex) = MathUtil.maxWithIndex(lls)
        val maxPosterior    = LogProbability.normalizeByLogProbability(maxLikelihood, likelihoodSum)
        val pConsensusError = LogProbability.not(maxPosterior) // convert to probability of the called consensus being wrong

        // Masks a base if the phred score would be too low
        consensusBases(baseIdx) = {
          if (PhredScore.fromLogProbability(pConsensusError) < this.options.minConsensusBaseQuality) SequenceUtil.N
          else DnaBasesUpperCase(maxLlIndex)
        }

        // Factor in the pre-UMI error rate.
        // Pr(error) = Pr(any pre-UMI error AND correct consensus) + Pr(no pre-UMI error AND any error in consensus)
        //               + Pr(pre-UMI error AND error in consensus, that do not give us the correct bases)
        // The last term tries to capture the case where a pre-UMI error modifies the base (ex. A->C) but a sequencing
        // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
        val p = probabilityOfErrorTwoTrials(this.options.errorRatePreUmiLn, pConsensusError)

        // Cap the quality
        consensusQualities(baseIdx) = PhredScore.cap(PhredScore.fromLogProbability(p))
      }
    }

    ConsensusRead(consensusBases, consensusQualities)
  }

  private[umi] def adjustBaseQualities(read: SourceRead): AdjustedRead = adjustBaseQualities(read=read, options=this.options)

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
    consensusFromSamRecords(records=fragments) match {
      case None       => // reject
        rejectRecords(records=fragments);
      case Some(frag) => // output
        this.createAndEnqueueSamRecord(records=fragments, read=frag, readName=nextReadName(fragments), readType=Fragment)
        success = true
    }

    // pairs
    val needBothPairs = options.requireConsensusForBothPairs // for readability later
    val firstOfPairConsensus  = consensusFromSamRecords(records=firstOfPair)
    val secondOfPairConsensus = consensusFromSamRecords(records=secondOfPair)
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
    if (records.isEmpty) ""
    else this.options.tag + ":" + records.head.getStringAttribute(this.options.tag)
//    val curIdx = readIdx
//    readIdx += 1
//    val prefix = readNamePrefix.getOrElse(longestCommonPrefix(records.map(_.getReadName)).getOrElse("CONSENSUS"))
//    s"$prefix:$curIdx"
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
    rec.setReadBases(read.bases)
    rec.setBaseQualities(read.quals)
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
    ConsensusRead(bases=source.bases.map(_ => 'N'.toByte), quals=source.quals.map(_ => PhredScore.MinValue))
  }
}
