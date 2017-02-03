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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.umi.UmiConsensusCaller._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import dagr.commons.util.Logger
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.{Murmur3, SequenceUtil}
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord, SAMRecord, SAMTag}

import scala.collection.mutable

/**
  * Contains shared types and functions used when writing UMI-driven consensus
  * callers that take in SAMRecords and emit SAMRecords.
  */
object UmiConsensusCaller {
  /** The type of consensus read to output. */
  object ReadType extends Enumeration {
    type ReadType = Value
    val Fragment, FirstOfPair, SecondOfPair = Value
  }

  /** A trait that consensus reads must implement. */
  trait SimpleRead {
    /** The ID of the molecule that generated the consensus read. */
    def id: String
    /** The bases of the consensus read. */
    def bases: Array[Byte]
    /** The quals of the consensus read. */
    def quals: Array[Byte]
    /** Gets the length of the consensus read. */
    def length: Int = bases.length
    /** Returns the consensus read a String - mostly useful for testing. */
    def baseString = new String(bases)
  }

  /** Stores information about a read to be fed into a consensus. */
  case class SourceRead(id: String, bases: Array[Byte], quals: Array[Byte]) extends SimpleRead {
    require(bases.length == quals.length, "Bases and qualities are not the same length.")
  }

  /**
    * Attempts to construct a String that can be used as a prefix for consensus read names based
    * on the contents of the incoming SAMFileHeader.
    */
  def makePrefixFromSamHeader(header: SAMFileHeader): String = {
    val ids = header.getReadGroups.map(rg => Option(rg.getLibrary).getOrElse(rg.getReadGroupId)).toList.sorted.distinct
    // Read names have to fit into 255 bytes, and this is just the prefix
    if (ids.map(_.length+1).sum <= 200) ids.mkString("|")
    else Integer.toHexString(new Murmur3(1).hashUnencodedChars(ids.mkString("|")))
  }

  /**
    * Constructs an output header with a single read group for the a BAM.
    */
  def outputHeader(in: SAMFileHeader, readGroupId: String, sortOrder: Option[SortOrder] = None): SAMFileHeader = {
    val oldRgs = in.getReadGroups.toSeq
    def collapse(f: SAMReadGroupRecord => String): String = oldRgs.map(f).filter(_ != null).distinct match {
      case Nil => null
      case xs  => xs.mkString(",")
    }

    val rg = new SAMReadGroupRecord(readGroupId)
    rg.setDescription     (collapse(_.getDescription))
    rg.setLibrary         (collapse(_.getLibrary))
    rg.setSample          (collapse(_.getSample))
    rg.setPlatform        (collapse(_.getPlatform))
    rg.setPlatformUnit    (collapse(_.getPlatformUnit))
    rg.setSequencingCenter(collapse(_.getSequencingCenter))

    val outHeader = new SAMFileHeader
    outHeader.addReadGroup(rg)
    outHeader.setSortOrder(sortOrder.getOrElse(SortOrder.unsorted))
    outHeader.setGroupOrder(GroupOrder.query)
    outHeader.addComment(s"Read group ${rg.getId} contains consensus reads generated from ${oldRgs.size} input read groups.")
    outHeader
  }
}


/**
  * A trait that can be mixed in by any consensus caller that works at the read level,
  * mapping incoming SAMRecords into consensus SAMRecords.
  *
  * @tparam C Internally, the type of lightweight consensus read that is used prior to
  *           rebuilding [[SAMRecord]]s.
  */
trait UmiConsensusCaller[C <: SimpleRead] {
  import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType._

  // vars to track how many reads meet various fates
  private var _totalReads: Long = 0
  private var _totalFiltered: Long = 0
  private var _filteredForMinorityAlignment: Long = 0
  private var _consensusReadsConstructed: Long = 0

  protected val NoCall: Byte = 'N'.toByte
  protected val NoCallQual: PhredScore = PhredScore.MinValue

  /** Returns the total number of input reads examined by the consensus caller so far. */
  def totalReads: Long = _totalReads

  /**
    * Returns the number of raw reads filtered out due to there being insufficient reads present
    * to build the necessary set of consensus reads.
    */
  def readsFilteredInsufficientSupport: Long = this._totalFiltered - this._filteredForMinorityAlignment

  /** Returns the number of raw reads filtered out because their alignment disagreed with the majority alignment of
    * all raw reads for the same source molecule.
    */
  def readsFilteredMinorityAlignment: Long = _filteredForMinorityAlignment

  /** Returns the number of consensus reads constructed by this caller. */
  def consensusReadsConstructed: Long = _consensusReadsConstructed

  /** Records that the supplied records were rejected, and not used to build a consensus read. */
  protected def rejectRecords(recs: Traversable[SAMRecord]) : Unit = {
    this._totalFiltered += recs.size
  }

  /** A RG.ID to apply to all generated reads. */
  protected def readGroupId: String

  /** A prefix to use on all read names.  If None then a suitable prefix will be synthesized. */
  protected val readNamePrefix: String

  /**
    * Needs to be implemented to return a value from a SAMRecord that represents the unit of
    * grouping, e.g. the MI tag for vanilla UMI data and the MI tag minus the /?? suffix for
    * duplex data.
    * @param rec a SAMRecord
    * @return an identified for the source molecule
    */
  protected[umi] def sourceMoleculeId(rec: SAMRecord): String

  /**
    * Converts from a SAMRecord into a SourceRead.
    */
  protected def toSourceRead(rec: SAMRecord, minBaseQuality: PhredScore): SourceRead = {
    val bases = rec.getReadBases.clone()
    val quals = rec.getBaseQualities.clone()

    if (rec.getReadNegativeStrandFlag) {
      SequenceUtil.reverseComplement(bases)
      SequenceUtil.reverse(quals, 0, quals.length)
    }

    forloop (from=0, until=bases.length) { i =>
      if (quals(i) < minBaseQuality) {
        bases(i) = NoCall
        quals(i) = NoCallQual
      }
    }

    SourceRead(sourceMoleculeId(rec), bases, quals)
  }

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SAMRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  final def consensusReadsFromSamRecords(recs: Seq[SAMRecord]): Seq[SAMRecord] = {
    this._totalReads += recs.size
    val result = consensusSamRecordsFromSamRecords(recs)
    this._consensusReadsConstructed += result.size
    result
  }

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SAMRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  protected def consensusSamRecordsFromSamRecords(recs: Seq[SAMRecord]): Seq[SAMRecord]

  /** Split records into those that should make a single-end consensus read, first of pair consensus read,
    * and second of pair consensus read, respectively.  The default method is to use the SAM flag to find
    * unpaired reads, first of pair reads, and second of pair reads.
    */
  protected def subGroupRecords(records: Seq[SAMRecord]): (Seq[SAMRecord], Seq[SAMRecord], Seq[SAMRecord]) = {
    val fragments    = records.filter { rec => !rec.getReadPairedFlag }
    val (firstOfPair, secondOfPair) = records.filter(_.getReadPairedFlag).partition(_. getFirstOfPairFlag)
    (fragments, firstOfPair, secondOfPair)
  }

  /**
    * Takes in a non-empty seq of SAMRecords and filters them such that the returned seq only contains
    * those reads that share the most common alignment of the read sequence to the reference.
    * If two or more different alignments share equal numbers of reads, the 'most common' will
    * be an arbitrary pick amongst those alignments, and the group of reads with that alignment will
    * be returned.
    *
    * For the purposes of this method all that is implied by "same alignment" is that any
    * insertions or deletions are at the same position and of the same length.  This is done
    * to allow for differential read length (either due to sequencing or untracked hard-clipping
    * of adapters) and for differential soft-clipping at the starts and ends of reads.
    *
    * NOTE: filtered out reads are sent to the [[rejectRecords()]] method and do not need further handling
    */
  protected[umi] def filterToMostCommonAlignment(recs: Seq[SAMRecord]): Seq[SAMRecord] = {
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
    val (keepers, rejects) = groups.values.partition(_.size == max)
    rejects.foreach { rejectGroup =>
      rejectRecords(rejectGroup)
      this._filteredForMinorityAlignment += rejectGroup.size
    }

    keepers.head
  }

  /** Creates a `SAMRecord` from the called consensus base and qualities. */
  protected def createSamRecord(read: C, readType: ReadType): SAMRecord = {
    val rec = new SAMRecord(null)
    rec.setReadName(this.readNamePrefix + ":" + read.id)
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
    rec.setAttribute(ConsensusTags.MolecularId, read.id)

    rec
  }

  /** Sums a short array into an Int to avoid overflow. */
  protected def sum(ss: Array[Short]): Int = {
    var total: Int = 0
    forloop(from=0, until=ss.length) { i => total += ss(i) }
    total
  }


  /**
    * Logs statistics about how many reads were seen, and how many were filtered/discarded due
    * to various filters.
    */
  def logStatistics(logger: Logger): Unit = {
    logger.info(f"Total Raw Reads Considered: ${totalReads}%,d.")
    logger.info(f"Raw Reads Filtered Due to Mismatching Alignments: ${readsFilteredMinorityAlignment}%,d.")
    logger.info(f"Raw Reads Filtered Due to Insufficient Raw Reads: ${readsFilteredInsufficientSupport}%,d.")
    logger.info(f"Consensus reads emitted: ${consensusReadsConstructed}%,d.")
  }
}
