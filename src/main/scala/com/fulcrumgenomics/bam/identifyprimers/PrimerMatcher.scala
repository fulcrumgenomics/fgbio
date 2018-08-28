/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.identifyprimers

import java.nio.ByteBuffer

import com.fulcrumgenomics.alignment.{Alignable, Aligner, Alignment, AlignmentTask}
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef.{forloop, unreachable}
import com.fulcrumgenomics.commons.util.LazyLogging
import htsjdk.samtools.util.{Interval, OverlapDetector, SequenceUtil}

import scala.collection.mutable

private[identifyprimers] sealed trait PrimerMatcher {
  /** The primers that are being matched against */
  def primers: Seq[Primer]

  /** For location-based matching, the maximum distance from the primer 5' start and the read's 5' start.  For gapped
    * alignment, the extra # of read bases to include. */
  def slop: Int

  /** Finds a primer match, if any. */
  def find(rec: SamRecord): Option[PrimerMatch]
}

/** A trait that can be mixed into [[PrimerMatcher]]s that provides a method `PrimerMatcherWithKmerFilter.getPrimersForAlignmentTasks()`
  * that returns only primers that share a kmer with the corresponding bases in the read. */
private[identifyprimers] sealed trait PrimerMatcherWithKmerFilter {

  /** The primers that are being matched against */
  def primers: Seq[Primer]

  /** The kmer-length for the cache. */
  def kmerLength: Option[Int]

  /** The maximum primer length across all primers.  This is the width of the search window in the read. */
  private val maxPrimerLength: Int = this.primers.map(_.length).max

  /** A map from a kmer to the list of primers that contain that kmer. */
  private val kmerToPrimers: Map[ByteBuffer, Set[Primer]] = kmerLength match {
    case None => Map.empty[ByteBuffer, Set[Primer]]
    case Some(len) =>
      val kmers = new mutable.TreeMap[ByteBuffer,mutable.HashSet[Primer]]()
      primers.foreach { primer =>
        // create all possible kmers from the primer
        val until = primer.length - len
        forloop (from = 0, until = until) { from =>
          // extract the kmer
          val kmer = ByteBuffer.wrap(primer.bases.slice(from = from, until = from + len))
          // add the primer to that the list for that kmer if the kmer is already in the map, otherwise, create a new entry
          kmers.get(kmer) match {
            case Some(b) => b.add(primer)
            case None    =>
              val buffer = new mutable.HashSet[Primer]()
              kmers.put(kmer, buffer)
              buffer.add(primer)

          }
        }
      }
      // convert to an immutable interface
      kmers.iterator.map { case (bytes, _primers) => bytes -> _primers.toSet }.toMap
  }


  /** Finds all [[Primer]]s that we should perform alignment against.  Does not use the kmer hash if it exists.
    *
    * If `ignorePrimerStrand` is true, then all primers and their reverse complements are returned.  Otherwise, if the
    * read is unmapped, try matching against all positive strand primers and the reverse complement of hte reverse strand
    * primers.  Otherwise, the read is mapped, so try only primers that map to the same strand.
    *
    * @param rec the records to check against
    * @param ignorePrimerStrand ignore the strand of the primers
    */
  private def getPrimersForAlignmentTasksNoHash(rec: SamRecord, ignorePrimerStrand: Boolean): Seq[Primer] = {
    (ignorePrimerStrand, rec.unmapped) match {
      case (true, _) =>
        // try all primers and the reverse complement of all the primers
        this.primers ++ this.primers.map(_.copy(reverseComplementBases = true))
      case (_, true) =>
        // We don't know if the read maps to the forward or reverse strand, so try primers with strands from both!
        // NB: the reverse primers have had their sequences reverse complemented already to be on the "top" genomic strand,
        // so we will need to reverse complement them again to get them back to sequencing order.
        this.primers.filter(_.positiveStrand) ++ this.primers.filter(_.negativeStrand).map(_.copy(reverseComplementBases = true))
      case _         => // (false, false)
        this.primers.filter(_.positiveStrand == rec.positiveStrand)
    }
  }

  /** Finds all [[Primer]]s that we should perform alignment against.  Uses the kmer hash.
    *
    * If `ignorePrimerStrand` is true:
    *   - If the read is unmapped or maps to the positive strand, match the read at the start against all primers and
    *   match against the end for the reverse complement of all primers.  If the read maps to the negative strand,
    *   match the read at the end against all primers and match against the start for the reverse complement of all
    *   primers.
    *
    * Otherwise, if the read is unmapped, match the read at the start against all forward strand primers and match
    * against the end for the reverse complement of all reverse strand primers
    *
    * Otherwise, the read is mapped, and so for positive strand reads we match the start of the read against forward
    * strand primers, and for negative strand reads we match the end of hte read against reverse strand primers.
    *
    * @param rec the records to check against
    * @param ignorePrimerStrand ignore the strand of the primers.s
    */
  private def getPrimersForAlignmentTasksWithHash(rec: SamRecord, ignorePrimerStrand: Boolean, kmerLength: Int): Seq[Primer] = {
    val _primers = mutable.HashSet[Primer]()
    val forward = rec.bases
    val revcomp = SequenceUtil.reverseComplement(rec.basesString).getBytes
    (ignorePrimerStrand, rec.unmapped) match {
      case (true, _) =>
        // If the read is unmapped or maps to the positive strand, match the read at the start against all primers and
        // match against the end for the reverse complement of all primers.  If the read maps to the negative strand,
        // match the read at the end against all primers and match against the start for the reverse complement of all primers.
        getPrimersFrom(bases = forward, matchStart = rec.unmapped || rec.positiveStrand, kmerLength = kmerLength)
          .foreach(_primers.add)
        getPrimersFrom(bases = revcomp, matchStart = rec.mapped && rec.negativeStrand, kmerLength = kmerLength)
          .foreach(_primers.add)
      case (_, true) =>
        // We don't know if the read maps to the forward or reverse strand, so try primers with strands from both!
        // NB: the reverse primers have had their sequences reverse complemented already to be on the "top" genomic strand,
        // so we will need to reverse complement them again to get them back to sequencing order.
        getPrimersFrom(bases = forward, matchStart = true, kmerLength = kmerLength)
          .filter(_.positiveStrand)
          .foreach(_primers.add)
        getPrimersFrom(bases = revcomp, matchStart = false, kmerLength = kmerLength)
          .filter(_.negativeStrand)
          .map(_.copy(reverseComplementBases = true)) // *** important *** reverse complement
          .foreach(_primers.add)
      case _ => // (false, false)
        // reads mapped to the forward should match on the start, and reads on the reverse should match the end
        getPrimersFrom(bases = forward, matchStart = rec.positiveStrand, kmerLength = kmerLength)
          .filter(_.positiveStrand == rec.positiveStrand)
          .foreach(_primers.add)
    }
    _primers.toSeq
  }

  /** Finds all [[Primer]]s containing any kmer from the 5' end of the read.
    *
    * See `getPrimersForAlignmentTasksNoHash` and `getPrimersForAlignmentTasksWithHash`.
    *
    * @param rec the records to check against
    * @param ignorePrimerStrand ignore the strand of the primers
    * @return
    */
  protected[identifyprimers] def getPrimersForAlignmentTasks(rec: SamRecord, ignorePrimerStrand: Boolean): Seq[Primer] = {
    kmerLength match {
      case None         => getPrimersForAlignmentTasksNoHash(rec, ignorePrimerStrand)
      case Some(length) => getPrimersForAlignmentTasksWithHash(rec, ignorePrimerStrand, length)
    }
  }

  /** Finds all primers that share a kmer within the range [start, end] of the bases. */
  private[identifyprimers] def getPrimersFrom(bases: Array[Byte], startOffset: Int, endOffset: Int, kmerLength: Int): Seq[Primer] = {
    require(bases.length >= kmerLength, s"bases.length: ${bases.length} kmerLength: $kmerLength")
    if (endOffset < startOffset) Seq.empty else {
      require(startOffset <= endOffset)
      Range.inclusive(start = startOffset, end = endOffset).flatMap { from =>
        val until = from + kmerLength
        val kmer  = bases.slice(from = from, until = until)
        require(kmer.length == kmerLength, s"from=$from until=$until startOffset=$startOffset endOffset=$endOffset")
        this.kmerToPrimers.getOrElse(ByteBuffer.wrap(kmer), Seq.empty)
      }
    }
  }

  /** Finds all primers that share a kmer with of the bases.
    *
    * @param bases the read bases
    * @param matchStart true if we are to search the start of the read, false to search the end
    * @param kmerLength the kmer length
    */
  private[identifyprimers] def getPrimersFrom(bases: Array[Byte], matchStart: Boolean, kmerLength: Int): Seq[Primer] = {
    if (bases.length < kmerLength) Seq.empty
    else if (matchStart) {
      val start = 0
      val end   = math.min(maxPrimerLength, bases.length) - kmerLength
      getPrimersFrom(bases, start, end, kmerLength)
    }
    else {
      val start = math.max(0, bases.length - maxPrimerLength)
      val end   = math.max(0, bases.length - kmerLength)
      getPrimersFrom(bases, start, end, kmerLength)
    }
  }
}

/** A little trait to facilitate ungapped alignment. */
private trait UngappedAlignment {
  /** The maximum per-base mismatch rate allowed for a [[PrimerMatch]]. */
  def maxMismatchRate: Double

  /** Returns a mismatch-based alignment match, otherwise None */
  protected[identifyprimers] def mismatchAlign(rec: SamRecord, primer: Primer): Option[UngappedAlignmentPrimerMatch] = {
    val bases          = rec.bases
    val maxMismatches  = math.min(bases.length, primer.length) * this.maxMismatchRate
    val matchReadStart = rec.unmapped || rec.positiveStrand
    val mm             = {
      if (matchReadStart) numMismatches(bases, primer.bases, 0, maxMismatches)
      else {
        val primerBases = if (bases.length < primer.length) primer.bases.drop(primer.length - bases.length) else primer.bases
        numMismatches(bases, primerBases, math.max(0, bases.length - primer.length), maxMismatches)
      }
    }
    if (mm <= maxMismatches) Some(UngappedAlignmentPrimerMatch(primer, mm, Int.MaxValue)) else None
  }

  /** Counts the # of mismatches, allowing the primer to have IUPAC bases. */
  private[identifyprimers] def numMismatches(bases: Array[Byte], primer: Array[Byte], basesStartIndex: Int, maxMismatches: Double): Int = {
    val length = math.min(bases.length - basesStartIndex, primer.length)
    var count = 0
    var primerIndex = 0
    while (primerIndex < length && count <= maxMismatches) { // empirically faster than a fgbio forloop
      if (!SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(bases(primerIndex + basesStartIndex), primer(primerIndex))) count += 1
      primerIndex += 1
    }
    count
  }
}


/** Finds primer matches based on mapping position.
  * @param primers the primers against which to match
  * @param slop the maximum distance from the primer 5' start and the read's 5' start.
  * @param ignorePrimerStrand ignores the strand of the primer.
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class LocationBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val ignorePrimerStrand: Boolean,
 val maxMismatchRate: Double) extends PrimerMatcher with UngappedAlignment with LazyLogging {
  import com.fulcrumgenomics.FgBioDef.javaIteratorAsScalaIterator

  private val fwdDetector: OverlapDetector[Primer] = {
    val d = new OverlapDetector[Primer](0, 0)
    this.primers.filter(_.positiveStrand).foreach { primer => d.addLhs(primer, primer) }
    d
  }

  private val revDetector: OverlapDetector[Primer] = {
    val d = new OverlapDetector[Primer](0, 0)
    this.primers.filter(_.negativeStrand).foreach { primer => d.addLhs(primer, primer) }
    d
  }

  private def matchStart(detector: OverlapDetector[Primer], interval: Interval, pos: Int): Iterator[Primer] = {
    detector.getOverlaps(interval).iterator()
      .filter(primer => math.abs(primer.start - pos) <= slop)
  }

  private def matchEnd(detector: OverlapDetector[Primer], interval: Interval, pos: Int): Iterator[Primer] = {
    detector.getOverlaps(interval).iterator()
      .filter(primer => math.abs(primer.end - pos) <= slop)
  }

  /** Gets forward (reverse) primers that overlap the read's start (end) position. The strand of the read is used
    * to determine against which primers to match (i.e. forward or reverse). Forward strand primers are searched for at
    * the start of the read, while reverse strand primers are searched for at the end of the read. */
  private[identifyprimers] def getOverlaps(rec: SamRecord): Iterator[Primer] = {
    require(rec.mapped)
    // NB: when `ignorePrimerStrand` is true, then when we are matching the start of the read against the reverse primers,
    // we must reverse complement the reverse primers, since they were reverse complemented when read in to make them
    // top genomic strand.  Similarly for matching the end of the read and forward primers.
    val fwdOverlaps: Iterator[Primer] = if (ignorePrimerStrand || rec.positiveStrand) {
      val pos      = rec.unclippedStart
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop)
      if (ignorePrimerStrand) matchStart(fwdDetector, interval, pos) ++ matchStart(revDetector, interval, pos).map(_.copy(reverseComplementBases = true))
      else matchStart(fwdDetector, interval, pos)
    } else Iterator.empty
    val revOverlaps: Iterator[Primer] = if (!ignorePrimerStrand || rec.negativeStrand) {
      val pos      = rec.unclippedEnd
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop)
      if (ignorePrimerStrand) matchEnd(fwdDetector, interval, pos).map(_.copy(reverseComplementBases = true)) ++ matchEnd(revDetector, interval, pos)
      else matchEnd(revDetector, interval, pos)
    } else Iterator.empty
    fwdOverlaps ++ revOverlaps
  }

  /** Returns a location-based match, otherwise None */
  def find(rec: SamRecord): Option[LocationBasedPrimerMatch] = if (rec.unmapped) None else {
    // Find primers that overlap the record, then prioritize by those matching the same strand.  Only one primer may be
    // returned.  In the future, it may be nice to either pick one randomly, keep them all, or prioritize by # of
    // mismatches.
    val primerMatches: Option[Primer] = getOverlaps(rec).toSeq match {
      case Seq()                       => None
      case Seq(primer)                 => Some(primer)
      case _primers                    =>
        _primers.filter(_.positiveStrand == rec.positiveStrand) match {
          case Seq()       => None
          case Seq(primer) => Some(primer)
          case ssPrimers   => unreachable(s"Found multiple primers on the same strand for $rec\n${ssPrimers.mkString("\n")}")
        }
    }

    primerMatches
      .flatMap { p => mismatchAlign(rec, p) }
      .map { m => // Convert to a [[LocationBasedPrimerMatch]]
        LocationBasedPrimerMatch(primer = m.primer, numMismatches = m.numMismatches)
      }
  }
}

/**
  *
  * @param primers the primers against which to match
  * @param slop the extra # of read bases to include when mapping gapped alignment tasks.
  * @param ignorePrimerStrand match against primers on the same strand as the read.
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class UngappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val ignorePrimerStrand: Boolean,
 val maxMismatchRate: Double,
 val kmerLength: Option[Int] = None
 ) extends PrimerMatcherWithKmerFilter with UngappedAlignment with LazyLogging {

  /** Examines all primers to find the best mismatch-based alignment match. */
  def find(rec: SamRecord): Option[UngappedAlignmentPrimerMatch] = {
    val alignments =  getPrimersForAlignmentTasks(rec, ignorePrimerStrand)
      .flatMap { p => mismatchAlign(rec, p) }
    getBestUngappedAlignment(alignments)
  }

  /** Examines all primers to find the best gapped-alignment-based match. */
  private[identifyprimers] def getBestUngappedAlignment(alignments: Seq[UngappedAlignmentPrimerMatch]): Option[UngappedAlignmentPrimerMatch] = {
    if (alignments.isEmpty) None
    else {
      val best: UngappedAlignmentPrimerMatch = alignments.minBy(_.numMismatches)
      alignments.filter(_ != best) match {
        case Seq() => Some(best.copy(nextNumMismatches = Int.MaxValue))
        case alns  => Some(best.copy(nextNumMismatches = alns.minBy(_.numMismatches).numMismatches))
      }
    }
  }
}

/** A little class that wraps a [[SamRecord]] to make it [[Alignable]] and also to defer extracting the bases to align.
  * NB: end is exclusive, and start and end are zero-based. */
private[identifyprimers] case class SamRecordAlignable(rec: SamRecord, offset: Int, override val length: Int) extends Alignable {
  def bases: Array[Byte] = rec.bases.slice(offset, offset + length)
}

/** A primer matcher.  Attempts to match based on location, then based on a mismatch-only alignment, then based on a
  * gapped-alignment.
  *
  * @param primers the primers against which to match
  * @param slop the extra # of read bases to include when mapping gapped alignment tasks.
  * @param ignorePrimerStrand match against primers on the same strand as the read.
  * @param aligner the aligner for gapped-alignment-based matching.
  * @param minAlignmentScoreRate the minimum per-base alignment score rate for a gapped alignment based match.  The
  *                              rate for alignment is the score divided by the primer length.
  * @param kmerLength require a read to share a k-mer of this length with a primer before performing gapped alignment
  */
private[identifyprimers] class GappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val ignorePrimerStrand: Boolean,
 val aligner: Aligner,
 val minAlignmentScoreRate: Double,
 val kmerLength: Option[Int] = None) extends PrimerMatcherWithKmerFilter with LazyLogging {

  /** Returns a gapped-alignment-based match, otherwise None */
  def find(rec: SamRecord): Option[GappedAlignmentPrimerMatch] = {
    toAlignmentTasks(rec) match {
      case Seq()          => None
      case alignmentTasks =>
        val alignments = alignmentTasks.map { task =>
          val alignment = this.aligner.align(query = task.queryBytes, target = task.targetBytes)
          AlignmentAndPrimer(alignment, task.query)
        }.toIndexedSeq
        getBestGappedAlignment(alignments)
    }
  }

  /** Creates a list of [[AlignmentTask]]s. */
  def toAlignmentTasks(rec: SamRecord): Seq[AlignmentTask[Primer, SamRecordAlignable]] = {
    // NB: we are aligning the full query to a sub-sequence of the record
    val samRecordAlignableCache = scala.collection.mutable.HashMap[(Int,Int), SamRecordAlignable]()
    val primersToCheck = getPrimersForAlignmentTasks(rec, ignorePrimerStrand)
    primersToCheck.map { primer =>
      val (targetOffset, targetLength) = {
        if (primer.positiveStrand) (0, math.min(primer.sequence.length + slop, rec.length))
        else {
          val offset = math.max(0, rec.length - primer.length - slop)
          val length = rec.length - offset
          (offset, length)
        }
      }
      val target = samRecordAlignableCache.getOrElseUpdate((targetOffset, targetLength), SamRecordAlignable(rec, targetOffset, targetLength))
      AlignmentTask(query = primer, target = target)
    }
  }

  // A little class to store an alignment to a given primer
  private case class AlignmentAndPrimer(alignment: Alignment, primer: Primer)

  /** A little helper method for [[find()]] to create a [[GappedAlignmentPrimerMatch]]. */
  private[identifyprimers] def getBestGappedAlignment(alignments: Seq[AlignmentAndPrimer]): Option[GappedAlignmentPrimerMatch] = {
    if (alignments.isEmpty) None
    else {
      val best: AlignmentAndPrimer = alignments.maxBy(_.alignment.score)
      alignments.filter(_ != best) match {
        case Seq() => Some(GappedAlignmentPrimerMatch(best.primer, best.alignment.score, (minAlignmentScoreRate * best.primer.length).toInt))
        case alns  => Some(GappedAlignmentPrimerMatch(best.primer, best.alignment.score, alns.maxBy(_.alignment.score).alignment.score))
      }
    }
  }
}