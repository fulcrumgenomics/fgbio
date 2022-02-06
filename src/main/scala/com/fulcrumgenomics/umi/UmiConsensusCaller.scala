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
import com.fulcrumgenomics.alignment.{Cigar, CigarElem}
import com.fulcrumgenomics.bam.{ClippingMode, SamRecordClipper}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.commons.util.{Logger, SimpleCounter}
import com.fulcrumgenomics.umi.UmiConsensusCaller._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.{Murmur3, SequenceUtil, TrimmingUtil}

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.math.min

/**
  * Contains shared types and functions used when writing UMI-driven consensus
  * callers that take in SamRecords and emit SamRecords.
  */
object UmiConsensusCaller {
  /** The type of consensus read to output. */
  object ReadType extends Enumeration {
    type ReadType = Value
    val Fragment, FirstOfPair, SecondOfPair = Value
    val ReadTypeKey: String = "__read_type__"
  }

  /** Filter reason for when there are too few reads to form a consensus. */
  val FilterInsufficientSupport = "Insufficient Support"

  /** Filter reason for when reads are rejected for having a minority CIGAR. */
  val FilterMinorityAlignment   = "Mismatching Cigars"

  /** Filter reason for when reads are rejected due to low quality. */
  val FilterLowQuality = "Low Base Quality"

  /** Filter reason for when reads are rejected due creation of orphaned consensus (i.e. R1 or R2 failed). */
  val FilterOrphan = "Orphan Consensus Created"

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
    def baseString: String = new String(bases)
    /** Retrieves the quals as a phred+33/fastq ascii String. */
    def qualString: String = SAMUtils.phredToFastq(this.quals)
  }

  /** Stores information about a read to be fed into a consensus. */
  case class SourceRead(id: String, bases: Array[Byte], quals: Array[Byte], cigar: Cigar, sam: Option[SamRecord] = None) extends SimpleRead {
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
  def outputHeader(in: SAMFileHeader, readGroupId: String, sortOrder: Option[SamOrder] = None): SAMFileHeader = {
    val oldRgs = in.getReadGroups.toSeq
    def collapse(f: SAMReadGroupRecord => String): String = oldRgs.map(f).filter(_ != null).distinct match {
      case Nil => null
      case xs  => xs.mkString(",")
    }

    val rg = new SAMReadGroupRecord(readGroupId)
    rg.setDescription     (collapse(_.getDescription))
    rg.setLibrary         (collapse(_.getLibrary))
    rg.setSample          (collapse(_.getSample))
    rg.setPlatform        (collapse(x => Option(x.getPlatform).map(_.toUpperCase).orNull)) // NB: this to ensure that platforms are all upper-case; not all tools are modern or well-behaved
    rg.setPlatformUnit    (collapse(_.getPlatformUnit))
    rg.setSequencingCenter(collapse(_.getSequencingCenter))

    val outHeader = new SAMFileHeader
    outHeader.addReadGroup(rg)
    sortOrder match {
      case Some(so) => so.applyTo(outHeader)
      case None     => SamOrder.Unsorted.applyTo(outHeader); outHeader.setGroupOrder(GroupOrder.query)
    }

    outHeader.addComment(s"Read group ${rg.getId} contains consensus reads generated from ${oldRgs.size} input read groups.")
    outHeader
  }

  /**
    * Helper method to check that the input BAM is in the correct order for consensus calling. Will call
    * the warn function if the sort order looks like it's probably compatible but it's not 100% sure.
    * Will invoke the error function in cases where the sort order is definitely incompatible.
    *
    * @param header the header of the BAM file to be used for consensus calling
    * @param source a path or string representing the source of the header
    * @param warn a function to be called when any warnings are detected/emitted
    * @param error a function to be called when any errors are encountered; should probably throw an exception!
    */
  def checkSortOrder(header: SAMFileHeader, source: Any, warn: String => Unit, error: String => Unit): Unit = {
    // Check that the SAM file is sorted appropriately
    if (!SamOrder(header).contains(SamOrder.TemplateCoordinate)) {
      // Group reads used to output the header without the sub-sort, so allow this for now
      if (header.getSortOrder == SortOrder.unsorted && header.getGroupOrder == GroupOrder.query) {
        warn(s"File $source may not be sorted correctly for consensus read generation.")
      }
      else {
        error(s"File $source is not sorted correctly. Please sort with fgbio SortBam -s TemplateCoordinate.")
      }
    }
  }
}


/**
  * A trait that can be mixed in by any consensus caller that works at the read level,
  * mapping incoming SamRecords into consensus SamRecords.
  *
  * @tparam ConsensusRead Internally, the type of lightweight consensus read that is used prior to
  *           rebuilding [[com.fulcrumgenomics.bam.api.SamRecord]]s.
  */
trait UmiConsensusCaller[ConsensusRead <: SimpleRead] {
  import com.fulcrumgenomics.umi.UmiConsensusCaller.ReadType._

  // vars to track how many reads meet various fates
  private var _totalReads: Long = 0
  private val _filteredReads = new SimpleCounter[String]()
  private var _consensusReadsConstructed: Long = 0

  protected val NoCall: Byte = 'N'.toByte
  protected val NoCallQual: PhredScore = PhredScore.MinValue

  /** A consensus caller used to generate consensus UMI sequences */
  private val consensusBuilder = new SimpleConsensusCaller()

  /** Clipper utility used to _calculate_ clipping, but not do the actual clipping */
  private val clipper = new SamRecordClipper(mode=ClippingMode.Soft, autoClipAttributes=true)

  /** Returns a clone of this consensus caller in a state where no previous reads were processed.  I.e. all counters
    * are set to zero.*/
  def emptyClone(): UmiConsensusCaller[ConsensusRead]

  /** Returns the total number of input reads examined by the consensus caller so far. */
  def totalReads: Long = _totalReads

  /** Returns the total number of reads filtered for any reason. */
  def totalFiltered: Long = _filteredReads.total

  /**
    * Returns the number of raw reads filtered out due to there being insufficient reads present
    * to build the necessary set of consensus reads.
    */
  def readsFilteredInsufficientSupport: Long = this._filteredReads.countOf(FilterInsufficientSupport)

  /** Returns the number of raw reads filtered out because their alignment disagreed with the majority alignment of
    * all raw reads for the same source molecule.
    */
  def readsFilteredMinorityAlignment: Long = this._filteredReads.countOf(FilterMinorityAlignment)

  /** Returns the number of consensus reads constructed by this caller. */
  def consensusReadsConstructed: Long = _consensusReadsConstructed

  /** Records that the supplied records were rejected, and not used to build a consensus read. */
  protected def rejectRecords(recs: Iterable[SamRecord], reason: String) : Unit = this._filteredReads.count(reason, recs.size)

  /** Records that the supplied records were rejected, and not used to build a consensus read. */
  protected def rejectRecords(reason: String, rec: SamRecord*) : Unit = rejectRecords(rec, reason)

  /** A RG.ID to apply to all generated reads. */
  protected def readGroupId: String

  /** A prefix to use on all read names.  If None then a suitable prefix will be synthesized. */
  protected val readNamePrefix: String

  /**
    * Needs to be implemented to return a value from a SamRecord that represents the unit of
    * grouping, e.g. the MI tag for vanilla UMI data and the MI tag minus the /?? suffix for
    * duplex data.
    * @param rec a SamRecord
    * @return an identified for the source molecule
    */
  protected[umi] def sourceMoleculeId(rec: SamRecord): String

  /**
    * Converts from a SamRecord into a SourceRead.  During conversion the record is end-trimmed
    * to remove Ns and bases below the `minBaseQuality`.  Remaining bases that are below
    * `minBaseQuality` are then masked to Ns.  Also trims reads so that no mapped bases extend past their mate.
    *
    * @return Some(SourceRead) if there are any called bases with quality > minBaseQuality, else None
    */
  protected[umi] def toSourceRead(rec: SamRecord, minBaseQuality: PhredScore, trim: Boolean): Option[SourceRead] = {
    // Extract and possibly RC the source bases and quals from the SamRecord
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    val cigar = if (rec.positiveStrand) rec.cigar else {
      SequenceUtil.reverseComplement(bases)
      SequenceUtil.reverse(quals, 0, quals.length)
      rec.cigar.reverse
    }

    // Quality trim the reads if requested.
    val trimToLength = if (trim) TrimmingUtil.findQualityTrimPoint(quals, minBaseQuality) else bases.length

    // Mask remaining low quality bases to Ns
    forloop (from=0, until=trimToLength) { i =>
      if (quals(i) < minBaseQuality) {
        bases(i) = NoCall
        quals(i) = NoCallQual
      }
    }

    // Get the length of the read based on trimming bases that are beyond the mate's end (FR only) and then any
    // remaining trailing Ns.  This includes both mapped bases and soft-clipped bases past the mate's end.
    val len = {
      var index = if (!rec.isFrPair) trimToLength - 1 else {
        // Get the number of mapped bases to clip that maps beyond the mate's end, including any soft-clipped bases. Use
        // that to compute where in the read to keep.
        val clipPosition = rec.length - this.clipper.numBasesExtendingPastMate(rec=rec)
        min(clipPosition, trimToLength) - 1
      }
      // Find the last non-N base of sufficient quality in the record, starting from either the
      // end of the read, or the end of the insert, whichever is shorter!
      while (index >= 0 && (bases(index) == NoCall)) index -= 1
      index + 1
    }

    len match {
      case 0 =>
        rejectRecords(Iterable(rec), FilterLowQuality)
        None
      case n if n == bases.length =>
        Some(SourceRead(sourceMoleculeId(rec), bases, quals, cigar, Some(rec)))
      case n                         =>
        val trimmedBases = new Array[Byte](n)
        val trimmedQuals = new Array[Byte](n)
        System.arraycopy(bases, 0, trimmedBases, 0, n)
        System.arraycopy(quals, 0, trimmedQuals, 0, n)
        Some(SourceRead(sourceMoleculeId(rec), trimmedBases, trimmedQuals, cigar.truncateToQueryLength(n), Some(rec)))
    }
  }

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SamRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  final def consensusReadsFromSamRecords(recs: Seq[SamRecord]): Seq[SamRecord] = {
    this._totalReads += recs.size
    val result = consensusSamRecordsFromSamRecords(recs)
    this._consensusReadsConstructed += result.size
    result
  }

  /**
    * Takes in all the reads for a source molecule and, if possible, generates one or more
    * output consensus reads as SAM records.
    *
    * @param recs the full set of source SamRecords for a source molecule
    * @return a seq of consensus SAM records, may be empty
    */
  protected def consensusSamRecordsFromSamRecords(recs: Seq[SamRecord]): Seq[SamRecord]

  /** Split records into those that should make a single-end consensus read, first of pair consensus read,
    * and second of pair consensus read, respectively.  The default method is to use the SAM flag to find
    * unpaired reads, first of pair reads, and second of pair reads.
    */
  protected def subGroupRecords(records: Seq[SamRecord]): (Seq[SamRecord], Seq[SamRecord], Seq[SamRecord]) = {
    val fragments    = records.filter { rec => !rec.paired }
    val (firstOfPair, secondOfPair) = records.filter(_.paired).partition(_.firstOfPair)
    (fragments, firstOfPair, secondOfPair)
  }

  /** Used it [[filterToMostCommonAlignment()]] to store a cigar string and a set of flags for which reads match. */
  private final case class AlignmentGroup(cigar: Cigar, flags: mutable.BitSet, var size: Int = 0) {
    /** Adds the read at `idx` to the set included. */
    @inline def add(idx: Int): Unit = {
      flags(idx) = true
      size      += 1
    }

    @inline def contains(idx: Int): Boolean = flags(idx)
  }

  /**
    * Takes in a non-empty seq of SamRecords and filters them such that the returned seq only contains
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
    * NOTE: filtered out reads are sent to the [[rejectRecords]] method and do not need further handling
    */
  protected[umi] def filterToMostCommonAlignment(recs: Seq[SourceRead]): Seq[SourceRead] = if (recs.size < 2) recs else {
    val groups = new ArrayBuffer[AlignmentGroup]
    val sorted = recs.sortBy(r => -r.length).toIndexedSeq

    forloop (from=0, until=sorted.length) { i =>
      val simpleCigar = simplifyCigar(sorted(i).cigar)
      var found = false
      groups.foreach { g => if (simpleCigar.isPrefixOf(g.cigar)) { g.add(i); found = true } }

      if (!found) {
        val newGroup = AlignmentGroup(simpleCigar, new mutable.BitSet(sorted.size))
        newGroup.add(i)
        groups += newGroup
      }
    }

    if (groups.isEmpty) {
      Seq.empty
    }
    else {
      val bestGroup = groups.maxBy(_.size)
      val keepers = new ArrayBuffer[SourceRead](bestGroup.size)
      forloop (from=0, until=sorted.length) { i =>
        if (bestGroup.contains(i)) keepers += sorted(i)
        else sorted(i).sam.foreach(rejectRecords(FilterMinorityAlignment, _))
      }

      keepers.toIndexedSeq
    }
  }

  /** Simplifies the cigar by turning other operators into Ms if that's how we want to think of them. */
  private def simplifyCigar(cigar: Cigar) = {
    import CigarOperator._
    if (cigar.forall(e => e.operator == M || e.operator == I || e.operator == D)) {
      cigar
    }
    else {
      val newElems = cigar.elems.map {
        case CigarElem(S,  len) => CigarElem(M, len)
        case CigarElem(EQ, len) => CigarElem(M, len)
        case CigarElem(X,  len) => CigarElem(M, len)
        case CigarElem(H,  len) => CigarElem(M, len)
        case cig => cig
      }

      Cigar(newElems).coalesce
    }
  }

  /** Creates a `SamRecord` from the called consensus base and qualities. */
  protected def createSamRecord(read: ConsensusRead, readType: ReadType, umis: Seq[String] = Seq.empty): SamRecord = {
    val rec = SamRecord(null)
    rec.name = this.readNamePrefix + ":" + read.id
    rec.unmapped = true
    readType match {
      case Fragment     => // do nothing
      case FirstOfPair  =>
        rec.paired       = true
        rec.firstOfPair  = true
        rec.mateUnmapped = true
      case SecondOfPair =>
        rec.paired       = true
        rec.secondOfPair = true
        rec.mateUnmapped = true
    }
    rec.bases = read.bases
    rec.quals = read.quals
    rec(SAMTag.RG.name()) = readGroupId
    rec(ConsensusTags.MolecularId) = read.id
    if (umis.nonEmpty) rec(ConsensusTags.UmiBases) = this.consensusBuilder.callConsensus(umis)

    rec
  }

  /** Sums a short array into an Int to avoid overflow. */
  protected def sum(ss: Array[Short]): Int = {
    var total: Int = 0
    forloop(from=0, until=ss.length) { i => total += ss(i) }
    total
  }

  /** Adds the given caller's statistics (counts) to this caller. */
  def addStatistics(caller: UmiConsensusCaller[ConsensusRead]): Unit = {
    this._totalReads += caller.totalReads
    this._consensusReadsConstructed += caller.consensusReadsConstructed
    this._filteredReads += caller._filteredReads
  }

  /**
    * Logs statistics about how many reads were seen, and how many were filtered/discarded due
    * to various filters.
    */
  def logStatistics(logger: Logger): Unit = {
    logger.info(f"Total Raw Reads Considered: $totalReads%,d.")
    this._filteredReads.foreach { case (filter, count) =>
      logger.info(f"Raw Reads Filtered Due to $filter: $count%,d (${count/totalReads.toDouble}%.4f).")
    }
    logger.info(f"Consensus reads emitted: $consensusReadsConstructed%,d.")
  }
}
