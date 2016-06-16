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
import com.fulcrumgenomics.util.PhredValue.{OneProbability, ZeroPhred, ZeroProbability}
import com.fulcrumgenomics.util.{PhredValue, ProgressLogger}
import dagr.commons.CommonsDef._
import htsjdk.samtools._
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object ConsensusCallerOptions {
  /** Various default values for the consensus caller. */
  val DefaultAttribute: String                       = "RX"
  val DefaultErrorRatePreUmi: PhredValue             = 50
  val DefaultErrorRatePostUmi: PhredValue            = 0
  val DefaultMaxBaseQuality: PhredValue              = 40
  val DefaultBaseQualityShift: PhredValue            = 10
  val DefaultMinConsensusBaseQuality: PhredValue     = 13
  val DefaultMinReads: Int                           = 1
  val DefaultMinMeanConsensusBaseQuality: PhredValue = 13
}

/** Holds the parameters/options for consensus calling. */
case class ConsensusCallerOptions(attribute: String                       = DefaultAttribute,
                                  errorRatePreUmi: PhredValue             = DefaultErrorRatePreUmi,
                                  errorRatePostUmi: PhredValue            = DefaultErrorRatePostUmi,
                                  maxBaseQuality: PhredValue              = DefaultMaxBaseQuality,
                                  baseQualityShift: PhredValue            = DefaultBaseQualityShift,
                                  minConsensusBaseQuality: PhredValue     = DefaultMinConsensusBaseQuality,
                                  minReads: Int                           = DefaultMinReads,
                                  minMeanConsensusBaseQuality: PhredValue = DefaultMinMeanConsensusBaseQuality
                                 )

object ConsensusCaller {
  /** Creates a consensus read from the given read and quality string tuples.  If no consensus read was created, None is returned. */
  def fromBasesAndQualities(basesAndQualities: Seq[(String, String)],
                            options: ConsensusCallerOptions = new ConsensusCallerOptions()
                           ): Option[(String, String)] = {
    if (basesAndQualities.exists { case (b, q) => b.length != q.length }) {
      throw new IllegalArgumentException("Found a read with the bases and qualities having different length")
    }
    // check to see if we have enough reads.
    if (basesAndQualities.length < options.minReads) return None

    // extract the bases and qualities, and adjust the qualities based on the given options.
    val bases: Seq[String] = basesAndQualities.map { _._1.toUpperCase }
    val quals: Seq[Seq[PhredValue]] = basesAndQualities.map(bq => SAMUtils.fastqToPhred(bq._2))
      .map { bq =>
        adjustBaseQualities(
          quals            = bq.map(q => PhredValue(q.toDouble)),
          maxBaseQuality   = options.maxBaseQuality,
          baseQualityShift = options.baseQualityShift,
          errorRatePostUmi = options.errorRatePostUmi
        )
      } // NB: binary Phred, not ASCII

    // get the most likely consensus bases and qualities
    val (consensusBases, consensusQualities) = consensusCalls(
      bases                   = bases,
      quals                   = quals,
      errorRatePreUmi         = options.errorRatePreUmi,
      minConsensusBaseQuality = options.minConsensusBaseQuality,
      minReads                = options.minReads
      )

    // check that the mean base quality is high enough
    if (PhredValue.mean(consensusQualities:_*) < options.minMeanConsensusBaseQuality) None
    else Some((consensusBases, SAMUtils.phredToFastq(consensusQualities.map(_.toByte).toArray[Byte])))
  }

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  def fromSamRecords(records: Seq[SAMRecord],
                     options: ConsensusCallerOptions = new ConsensusCallerOptions()
                    ): Option[(String, String)] = {
    // TODO: if mapped, this gets the bases in genomic order and not sequencing order.  Does it matter?
    val bases = records.map { rec => rec.getReadString.toUpperCase }
    val quals = records.map { rec => rec.getBaseQualityString }

    fromBasesAndQualities( basesAndQualities=bases.zip(quals), options = options)
  }

  /** Get the most likely consensus bases and qualities. */
  private[bam] def consensusCalls(bases: Seq[String],
                                  quals: Seq[Seq[PhredValue]],
                                  errorRatePreUmi: PhredValue         = DefaultErrorRatePreUmi,
                                  minConsensusBaseQuality: PhredValue = DefaultMinConsensusBaseQuality,
                                  minReads: Int                       = DefaultMinReads): (String, Seq[PhredValue]) = {
    val maxBaseLength = bases.map(_.length).max

    // go through each position in the read(s).
    val (consensusBases, consensusQualities) = Seq.range(0, maxBaseLength, 1).map { baseIdx =>
      val basesAtIdx = bases.filter { str => baseIdx < str.length}.map { str => str.charAt(baseIdx) }

      // check to see if we heave a enough reads to produce a consensus base.
      if (basesAtIdx.length < minReads) (SequenceUtil.N.toChar, ZeroPhred)
      else {
        val qualsAtIdx         = quals.filter { q => baseIdx < q.length}.map { q => q(baseIdx) }
        val basesAndQualsAtIdx = basesAtIdx.zip(qualsAtIdx)

        // TODO: handle Ns in the reads
        // get the phred-scaled probability of each potential DNA base
        var phreds = SequenceUtil.VALID_BASES_UPPER.map(_.toChar).map {
          candidateBase =>
            basesAndQualsAtIdx.map { case (base, qual) =>
              if (base == candidateBase) qual.inv() // probability of not an error
              else qual / 3.0 // probability of an error for this specific base
            }.foldLeft(OneProbability)((phred, qual) => phred * qual)
        }
        // normalize
        val phredsSum = phreds.foldLeft(ZeroProbability)((a, b) => a + b)
        phreds = phreds.map { phred => phred / phredsSum }

        // Find the base with the maximum probability, which means the smallest phred
        val minPhred    = phreds.min
        val minPhredIdx = phreds.indexOf(minPhred)
        val baseQual    = minPhred.inv() // convert the probability of the base to the phred-scaled probability of *not* the base
        // Masks a base if the phred quality is too low
        if (baseQual.toInt < minConsensusBaseQuality) (SequenceUtil.N.toChar, baseQual)
        else (SequenceUtil.VALID_BASES_UPPER(minPhredIdx).toChar, baseQual)
      }
    }.unzip

    // Factor in the pre-UMI error rate.
    val consensusQualitiesScaled = consensusQualities.map { phred =>
      // Pr(error) = Pr(any pre-UMI error AND correct consensus) + Pr(no pre-UMI error AND any error in consensus)
      //               + Pr(pre-UMI error AND error in consensus, that do not give us the correct bases)
      // The last term tries to capture the case where a pre-UMI error modifies the base (ex. A->C) but a sequencing
      // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
      val p = probabilityOfErrorTwoTrials(errorRatePreUmi, phred)
      PhredValue(Math.min(SAMUtils.MAX_PHRED_SCORE, p.value)) // cap the quality
    }
    (consensusBases.mkString, consensusQualitiesScaled)
  }

  /** Adjusts the given base qualities.  The base qualities are first shifted by `baseQualityShift`, then capped using
    * `maxBaseQuality`, and finally the `errorRatePostUmi` is incorporated.
    */
  private[bam] def adjustBaseQualities(quals: Seq[PhredValue],
                                       maxBaseQuality: PhredValue   = DefaultMaxBaseQuality,
                                       baseQualityShift: PhredValue = DefaultBaseQualityShift,
                                       errorRatePostUmi: PhredValue = DefaultErrorRatePostUmi
                                      ): Seq[PhredValue] = {
    quals.map { qual =>
      // shift the base qualities, then cap it
      PhredValue(Math.max(0.0, qual.value - baseQualityShift.value)) match {
        case q if maxBaseQuality < q => maxBaseQuality
        case q => q
      }
    }.map { phred =>
      // Pr(err) = Pr(no post-UMI error AND any sequencing error) + Pr(any post-UMI error and no sequencing error) +
      //             Pr(post-UMI error AND sequencing error)
      // The last term tries to capture the case where a post-UMI error modifies the base (ex. A->C) but a sequencing
      // error calls it the correct base (C->A).  Only 2/3 times will the two errors result in the incorrect base.
      probabilityOfErrorTwoTrials(errorRatePostUmi, phred)
    }
  }

  /** Computes the probability of seeing an error in the base sequence if there are two independent error processes.
    * We sum three terms:
    * 1. the probability of an error in trial one and no error in trial two: Pr(A=Error, B=NoError).
    * 2. the probability of no error in trial one and an error in trial two: Pr(A=NoError, B=Error).
    * 3. the probability of an error in both trials, but when the second trial does not reverse the error in first one, which
    *    for DNA (4 bases) would only occur 2/3 times: Pr(A=x->y, B=y->z) * Pr(x!=z | x!=y, y!=z, x,y,z \in {A,C,G,T})
    */
  private[bam] def probabilityOfErrorTwoTrials(prErrorTrialOne: PhredValue, prErrorTrialTwo: PhredValue): PhredValue = {
    val pr1 = prErrorTrialOne       * prErrorTrialTwo.inv()
    val pr2 = prErrorTrialOne.inv() * prErrorTrialTwo
    val pr3 = prErrorTrialOne       * prErrorTrialTwo       * 2.0 / 3.0
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
      val firsts = strs.map(_.head)
      if (firsts.toSet.size > 1) accu
      else longestCommonPrefix(strs.map(_.tail), Some(accu.getOrElse("") + firsts.head))
    }
  }
}

/** Calls consensus reads by grouping consecutive reads with the same attribute.
  *
  * Consecutive reads with the sam attribute are partitioned into fragments, first of pair, and
  * second of pair reads, and a consensus read is created for each partition.  A consensus read
  * for a given partition may not be returned if any of the conditions are not met (ex. minimum
  * number of reads, minimum mean consensus base quality, ...).
  * */
class ConsensusCaller
( input: Iterator[SAMRecord],
  val header: SAMFileHeader,
  val readNamePrefix: Option[String]   = None,
  val readGroup: String                = "A",
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
  def hasNext(): Boolean = {
    if (this.nextConsensusRecords.nonEmpty) true
    else this.advance()
  }

  /** Returns the next consensus read. */
  def next(): SAMRecord = {
    if (!this.hasNext()) throw new IllegalStateException("Calling next() when hasNext() is false.")
    this.nextConsensusRecords match {
      case mutable.Queue(x, _*) => this.nextConsensusRecords.dequeue()
      case mutable.Queue() => unreachable("hasNext() was true but nextConsensusRecords is empty")
    }
  }

  /** Consumes records until a consensus read can be created, or no more input records are available. Returns
    * true if a consensus read was created, false otherwise. */
  private def advance(): Boolean = {
    if (this.iter.isEmpty) return false

    // get the records to create the consensus read
    val buffer = getNextGroupOfRecords(attr=this.iter.head.toAttr)

    // partition the records to which end of a pair it belongs, or if it is a fragment read.
    val (fragments, firstOfPair, secondOfPair) = subGroupRecords(records = buffer)

    // build a consensus for each sub-group that has a read
    val fragmentName = nextReadName(fragments)
    val pairName = nextReadName(firstOfPair ++ secondOfPair) // must have the same read name
    var success = false
    success |= this.callConsensus(fragments,    readName=fragmentName, readType=Fragment)
    success |= this.callConsensus(firstOfPair,  readName=pairName,     readType=FirstOfPair)
    success |= this.callConsensus(secondOfPair, readName=pairName,     readType=SecondOfPair)

    if (success) true
    else if (this.iter.nonEmpty) this.advance() // try again
    else false // failed to build a consensus, and no more records in the input iterator
  }

  /** Adds a SAM record from the underlying iterator to the buffer if either the buffer is empty or the SAM attribute is
    * the same for the records in the buffer as the next record in the input iterator.  Returns true if a record was
    * added, false otherwise.
    */
  @annotation.tailrec
  private def getNextGroupOfRecords(buffer: ListBuffer[SAMRecord] = new ListBuffer[SAMRecord], attr: String): ListBuffer[SAMRecord] = {
    if (this.iter.isEmpty || this.iter.head.toAttr != attr) return buffer
    val rec = this.iter.next()
    progress.map(_.record(rec))
    buffer.append(rec)
    getNextGroupOfRecords(buffer=buffer, attr=attr)
  }

  /** Adds a method to the SAMRecord class to retrieve the attribute for the UMI identifier. */
  implicit private class SAMRecordToString(r: SAMRecord) {
    def toAttr: String = r.getStringAttribute(options.attribute)
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

  /** Attempts to create a consensus read from the given records adds it to `nextConsensusRecords`. Returns true if a
    * consensus read was created, false otherwise. */
  private def callConsensus(records: Seq[SAMRecord],
                   readName: String,
                   readType: ReadType.Value): Boolean = {
    ConsensusCaller.fromSamRecords(records=records, options=options) match {
      case Some((bases: String, quals: String)) =>
        // create the SAM record
        this.nextConsensusRecords.enqueue(
          this.createSamRecord(
            records  = records,
            bases    = bases,
            quals    = quals,
            readName = readName,
            readType = readType)
        )
        true
      case None =>
        this.rejects.foreach(rej => records.foreach(rej.addAlignment))
        false
    }
  }

  /** Creates a `SAMRecord` from the called consensus base and qualities. */
  private def createSamRecord(records: Seq[SAMRecord],
                              bases: String,
                              quals: String,
                              readName: String,
                              readType: ReadType.Value): SAMRecord = {
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
    rec.setReadString(bases)
    rec.setBaseQualityString(quals)
    rec.setAttribute(SAMTag.RG.name(), readGroup)
    rec.setAttribute(options.attribute, records.head.toAttr)
    // TODO: set custom SAM tags:
    // - # of reads contributing to this consensus
    rec
  }
}
