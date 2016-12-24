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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes._
import com.fulcrumgenomics.util.{MathUtil, ProgressLogger}
import dagr.commons.util.LazyLogging
import htsjdk.samtools._
import htsjdk.samtools.util.{Murmur3, SequenceUtil}

import scala.collection.JavaConversions.collectionAsScalaIterable
import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object VanillaUmiConsensusCallerOptions {
  type PhredScore = Byte

  /** Various default values for the consensus caller. */
  val DefaultTag: String                             = "MI"
  val DefaultErrorRatePreUmi: PhredScore             = 45.toByte
  val DefaultErrorRatePostUmi: PhredScore            = 40.toByte
  val DefaultMinInputBaseQuality: PhredScore         = 2.toByte
  @deprecated("Too redundant with DefaultErrorRatePostUmi.", since="0.1.2") val DefaultMaxBaseQuality: PhredScore = 40.toByte
  @deprecated("Too redundant with DefaultErrorRatePostUmi.", since="0.1.2") val DefaultBaseQualityShift: PhredScore = 10.toByte
  val DefaultMinConsensusBaseQuality: PhredScore     = 13.toByte
  val DefaultMinReads: Int                           = 1
  @deprecated("Too redundant with DefaultMinConsensusBaseQuality.", since="0.1.2") val DefaultMinMeanConsensusBaseQuality: PhredScore = 13.toByte
  val DefaultRequireConsensusForBothPairs: Boolean   = true
}

/** Holds the parameters/options for consensus calling. */
case class VanillaUmiConsensusCallerOptions
(
  tag: String                             = DefaultTag,
  errorRatePreUmi: PhredScore             = DefaultErrorRatePreUmi,
  errorRatePostUmi: PhredScore            = DefaultErrorRatePostUmi,
  minInputBaseQuality: PhredScore         = DefaultMinInputBaseQuality,
  @deprecated("Use errorRatePostUmi instead.", since="0.1.2") maxRawBaseQuality: PhredScore           = DefaultMaxBaseQuality,
  @deprecated("Use errorRatePostUmi instead.", since="0.1.2") rawBaseQualityShift: PhredScore         = DefaultBaseQualityShift,
  minConsensusBaseQuality: PhredScore     = DefaultMinConsensusBaseQuality,
  minReads: Int                           = DefaultMinReads,
  @deprecated("Use minConsensusBaseQuality instead.", since="0.1.2") minMeanConsensusBaseQuality: PhredScore = DefaultMinMeanConsensusBaseQuality,
  requireConsensusForBothPairs: Boolean   = DefaultRequireConsensusForBothPairs
)

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

object VanillaUmiConsensusCaller {
  /**
    * Takes in a seq of SAMRecords and filters them such that the returned seq only contains
    * those reads that share the most common alignment of the read sequence to the reference.
    * If two or more different alignments share equal numbers of reads, the 'most common' will
    * be an arbitrary pick amongst those alignments, and the group of reads with that alignment will
    * be returned.
    *
    * For the purposes of this method all that is implied by "same alignment" is that any
    * insertions or deletions are at the same position and of the same length.  This is done
    * to allow for differential read length (either due to sequencing or untracked hard-clipping
    * of adapters) and for differential soft-clipping at the starts and ends of reads.
    */
  def filterToMostCommonAlignment(recs: Seq[SAMRecord]): Seq[SAMRecord] = {
    val groups = recs.groupBy { r =>
      val builder = new mutable.StringBuilder
      val elems = r.getCigar.getCigarElements.iterator().bufferBetter

      // This loop constructs a pseudo-cigar for each read.  The result is a string that contains, for each indel
      // operator in the reads's cigar, the number of non-indel bases in the read since the start/last indel, then
      // the length of the indel.  E.g.:
      //     50M         => <empty> because there are no indels!
      //     10M2D10M    => 10M2D   because we don't care about anything after the last indel
      //     5H5S5M1I40M => 15M1I   because all the leading stuff gets collapsed
      while (elems.nonEmpty) {
        // Consume and merge all non-indel operators (clip, match/mismatch) at the head of the
        // iterator since we only care where the indels occur relative to the read
        val len = elems.takeWhile(e => !e.getOperator.isIndelOrSkippedRegion).map(_.getLength).sum

        // Now the iterator is either at the end of the read, in which case we don't need to do anything,
        // or it is before an I or D cigar entry, so we should output how much read it took to get here
        // followed by the I or D element
        if (elems.nonEmpty) {
          val indel = elems.next()
          if (len > 0) builder.append(len).append("M")
          builder.append(indel.toString)
        }
      }

      builder.toString()
    }

    val max = groups.values.map(_.size).max
    groups.values.find(_.size == max).getOrElse(unreachable("At least one group must have the max size!"))
  }
}


/** Calls consensus reads by grouping consecutive reads with the same SAM tag.
  *
  * Consecutive reads with the SAM tag are partitioned into fragments, first of pair, and
  * second of pair reads, and a consensus read is created for each partition.  A consensus read
  * for a given partition may not be returned if any of the conditions are not met (ex. minimum
  * number of reads, minimum mean consensus base quality, ...).
  * */
class VanillaUmiConsensusCaller
(input: Iterator[SAMRecord],
 val header: SAMFileHeader,
 val readNamePrefix: Option[String]   = None,
 val readGroupId: String              = "A",
 val options: VanillaUmiConsensusCallerOptions  = new VanillaUmiConsensusCallerOptions(),
 val rejects: Option[SAMFileWriter]   = None,
 val progress: Option[ProgressLogger] = None
) extends Iterator[SAMRecord] with LazyLogging {
  /** The type of consensus read to output. */
  private object ReadType extends Enumeration {
    val Fragment, FirstOfPair, SecondOfPair = Value
  }
  import ReadType._

  private val DnaBasesUpperCase: Array[Byte] = Array('A', 'C', 'G', 'T').map(_.toByte)
  private val LogThree = LogProbability.toLogProbability(3.0)
  private val caller = new ConsensusCaller(errorRatePreLabeling  = options.errorRatePreUmi,
                                           errorRatePostLabeling = options.errorRatePostUmi,
                                           maxRawBaseQuality     = options.maxRawBaseQuality,
                                           rawBaseQualityShift   = options.rawBaseQualityShift)

  private val iter = input.buffered
  private val nextConsensusRecords: mutable.Queue[SAMRecord] = mutable.Queue[SAMRecord]() // one per UMI group
  private var readIdx = 1

  private val NoCall: Base = 'N'.toByte
  private val NotEnoughReadsQual: PhredScore = 0.toByte // Score output when masking to N due to insufficient input reads
  private val TooLowQualityQual: PhredScore = 2.toByte  // Score output when masking to N due to too low consensus quality

  // var to track how many reads meet various fates
  private var _totalReads: Long = 0
  private var _filteredForMinReads: Long = 0
  private var _filteredForMinorityAlignment: Long = 0

  // Initializes the prefix that will be used to make sure read names are (hopefully) unique
  private val actualReadNamePrefix = readNamePrefix.getOrElse {
    val ids = header.getReadGroups.map(rg => Option(rg.getLibrary).getOrElse(rg.getReadGroupId)).toList.sorted.distinct
    // Read names have to fit into 255 bytes, and this is just the prefix
    if (ids.map(_.length+1).sum <= 200) ids.mkString("|")
    else Integer.toHexString(new Murmur3(1).hashUnencodedChars(ids.mkString("|")))
  }

  /** Returns the total number of input reads examined by the consensus caller so far. */
  def totalReads: Long = _totalReads

  /** Returns the number of raw reads filtered out due to the minReads threshold not being met. */
  def readsFilteredForMinReads: Long = _filteredForMinReads

  /** Returns the number of raw reads filtered out because their alignment disagreed with the majority alignment of
    * all raw reads for the same source molecule
    */
  def readsFilteredMinorityAlignment: Long = _filteredForMinorityAlignment

  /** Creates a consensus read from the given records.  If no consensus read was created, None is returned. */
  def consensusFromSamRecords(records: Seq[SAMRecord]): Option[ConsensusRead] = {
    this._totalReads += records.size

    if (records.size < this.options.minReads) {
      this._filteredForMinReads += records.size
      None
    }
    else {
      val filteredRecords = VanillaUmiConsensusCaller.filterToMostCommonAlignment(records)
      this._filteredForMinorityAlignment += (records.size - filteredRecords.size)

      if (filteredRecords.size < records.size) {
        val r = records.head
        val n = if (r.getReadPairedFlag && r.getSecondOfPairFlag) "/2" else "/1"
        val m = r.getAttribute(this.options.tag)
        val discards = records.size - filteredRecords.size
        logger.debug("Discarded ", discards, "/", records.size, " records due to mismatched alignments for ", m, n)
      }

      val sourceReads = filteredRecords.map { rec =>
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
  }

  /** Creates a consensus read from the given read and qualities sequences.  If no consensus read was created, None
    * is returned.
    *
    * The same number of base sequences and quality sequences should be given.
    * */
  private[umi] def consensusCall(reads: Seq[SourceRead]): Option[ConsensusRead] = {
    // check to see if we have enough reads.
    if (reads.size < this.options.minReads) {
      None
    }
    else {
      // get the most likely consensus bases and qualities
      val consensusLength = consensusReadLength(reads)
      val consensusBases  = new Array[Base](consensusLength)
      val consensusQuals  = new Array[PhredScore](consensusLength)

      var positionInRead = 0
      val builder = this.caller.builder()
      while (positionInRead < consensusLength) {
        // Add the evidence from all reads that are long enough to cover this base
        reads.foreach { read =>
          if (read.length > positionInRead) {
            val base = read.bases(positionInRead)
            val qual = read.quals(positionInRead)
            if (qual >= this.options.minInputBaseQuality) {
              builder.add(base=base, qual=qual)
            }
          }
        }

        // Call the consensus and do any additional filtering
        val (base, qual) = {
          if (builder.contributions < this.options.minReads) {
            (NoCall, NotEnoughReadsQual)
          }
          else {
            val (b,q) = builder.call()
            if (q < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual) else (b,q)

          }
        }
        consensusBases(positionInRead) = base
        consensusQuals(positionInRead) = qual

        // Get ready for the next pass
        builder.reset()
        positionInRead += 1
      }

      // check that the mean base quality is high enough
      if (MathUtil.mean(consensusQuals) >= options.minMeanConsensusBaseQuality)
        Some(ConsensusRead(consensusBases, consensusQuals))
      else
        None
    }
  }

  /**
    * Calculates the length of the consensus read that should be produced. The length is calculated
    * as the maximum length at which #options.minReads reads still have bases.
    *
    * @param reads the set of reads being fed into the consensus
    * @return the length of consensus read that should be created
    */
  private def consensusReadLength(reads: Seq[SourceRead]): Int = {
    val n = this.options.minReads
    if (reads.length < n) throw new IllegalArgumentException("Too few reads to create a consensus.")

    reads.map(_.length).sortBy(len => -len).drop(n-1).head
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

  /** Returns the next read name with format "<prefix>:<idx>", where "<prefix>" is either the supplied prefix or a
    * prefix composed by concatenating information from the input read groups.
    */
  private def nextReadName(records: Seq[SAMRecord]): String = {
    if (records.isEmpty) unreachable("Can't generate a consensus read name for zero input reads.")
    else this.actualReadNamePrefix + ":" + records.head.getStringAttribute(this.options.tag)
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
