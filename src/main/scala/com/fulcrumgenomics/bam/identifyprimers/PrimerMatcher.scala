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

private[identifyprimers] sealed trait _PrimerMatcher {
  /** The primers that are being matched against */
  def primers: Seq[Primer]

  /** For location-based matching, the maximum distance from the primer 5' start and the read's 5' start.  For gapped
    * alignment, the extra # of read bases to include. */
  def slop: Int

  /** Finds a primer match, if any. */
  def find(rec: SamRecord): Option[PrimerMatch]
}

private[identifyprimers] sealed trait _PrimerMatcherWithCache extends _PrimerMatcher {

  private val maxPrimerLength: Int = this.primers.map(_.length).max

  /** A map from a kmer to the list of primers that contain that kmer. */
  private val kmerToPrimers: Map[ByteBuffer, Set[Primer]] = kmerLength match {
    case None => Map.empty[ByteBuffer, Set[Primer]]
    case Some(len) =>
      val kmers = new mutable.TreeMap[ByteBuffer,mutable.HashSet[Primer]]()
      primers.foreach { primer =>
        val until = primer.length - len
        forloop (from = 0, until = until) { from =>
          val kmer = ByteBuffer.wrap(primer.bases.slice(from = from, until = from + len))
          kmers.get(kmer) match {
            case Some(b) => b.add(primer)
            case None    =>
              val buffer = new mutable.HashSet[Primer]()
              kmers.put(kmer, buffer)
              buffer.add(primer)

          }
        }
      }
      kmers.iterator.map { case (bytes, _primers) => bytes -> _primers.toSet }.toMap
  }

  def kmerLength: Option[Int]

  def requireSameStrand: Boolean

  /** Creates a list of [[AlignmentTask]]s. */
  def toAlignmentTasks(rec: SamRecord): Seq[AlignmentTask[Primer, SamRecordAlignable]] = {
    // NB: we are aligning the full query to a sub-sequence of the record
    val samRecordAlignable = scala.collection.mutable.HashMap[(Int,Int), SamRecordAlignable]()
    val primersToCheck = getPrimersForAlignmentTasks(rec)
    primersToCheck.map { primer =>
      if (primer.positiveStrand) {
        val targetOffset = 0
        val targetLength = math.min(primer.sequence.length + slop, rec.length) - targetOffset
        val target       = samRecordAlignable.getOrElseUpdate((targetOffset, targetLength), SamRecordAlignable(rec, targetOffset, targetLength))
        AlignmentTask(query = primer, target = target)
      }
      else {
        val targetOffset = math.max(0, rec.length - primer.length - slop)
        val targetLength = rec.length - targetOffset
        val target       = samRecordAlignable.getOrElseUpdate((targetOffset, targetLength), SamRecordAlignable(rec, targetOffset, targetLength))
        AlignmentTask(query=primer, target=target)
      }
    }
  }

  /** Finds all [[Primer]]s containing any kmer from the 5' end of the read. */
  protected def getPrimersForAlignmentTasks(rec: SamRecord): Seq[Primer] = {
    kmerLength match {
      case None         => this.primers
      case Some(length) =>
        val _primers = mutable.HashSet[Primer]()
        if (rec.unmapped) {
          val forward = rec.bases
          val reverse = SequenceUtil.reverseComplement(rec.basesString).getBytes
          val both    = Seq(forward, reverse)
          both.foreach { bases =>
            getPrimersFrom(bases = bases, matchStart = true, matchEnd = true, kmerLength = length).foreach(_primers.add)
          }
        }
        else {
          // Reads mapped to the positive strand should match primers at the start of the read, while negative strand
          // reads should match at the end of their reads.  If we require to match only primers from the same strand,
          // then keep only forward primers if the read was mapped to the + strand etc.
          val primers = getPrimersFrom(
            bases      = rec.bases,
            matchStart = rec.positiveStrand,
            matchEnd   = rec.negativeStrand,
            kmerLength = length
          )
          (if (requireSameStrand) primers.filter(_.positiveStrand == rec.positiveStrand) else primers).foreach(_primers.add)
        }
        _primers.toSeq
    }
  }

  /** Finds all primers that share a kmer within the range [start, end] of the bases. */
  private def getPrimersFrom(bases: Array[Byte], start: Int, end: Int, kmerLength: Int): Seq[Primer] = if (end - start + 1 < kmerLength) Seq.empty else {
    Range.inclusive(start = start, end = end).flatMap { from =>
      val until = from + kmerLength
      val kmer  = bases.slice(from = from, until = until)
      this.kmerToPrimers.getOrElse(ByteBuffer.wrap(kmer), Seq.empty)
    }
  }

  /** Finds all primers that share a kmer with of the bases.
    *
    * @param bases the read bases
    * @param matchStart true if we are to search the start of the read, false otherwise
    * @param matchEnd true if we are to search the end of hte read, false otherwise
    * @param kmerLength the kmer length
    */
  private def getPrimersFrom(bases: Array[Byte], matchStart: Boolean, matchEnd: Boolean, kmerLength: Int): Seq[Primer] = {
    val startPrimers = if (matchStart) {
      val start = 0
      val end   = math.min(maxPrimerLength, bases.length) - kmerLength
      getPrimersFrom(bases, start, end, kmerLength)
    } else Seq.empty
    val endPrimers = if (matchEnd) {
      val start = math.max(0, bases.length - maxPrimerLength)
      val end   = math.max(0, bases.length - kmerLength)
      getPrimersFrom(bases, start, end, kmerLength)
    } else Seq.empty
    startPrimers ++ endPrimers
  }
}

trait _MismatchAlign {
  /** The maximum per-base mismatch rate allowed for a [[PrimerMatch]]. */
  def maxMismatchRate: Double

  /** Returns a mismatch-based alignment match, otherwise None */
  protected[identifyprimers] def mismatchAlign(bases: Array[Byte], primer: Primer): Option[UngappedAlignmentPrimerMatch] = {
    val maxMismatches = math.min(bases.length, primer.length) * this.maxMismatchRate
    val mm = {
      if (primer.positiveStrand) numMismatches(bases, primer.bases, 0, maxMismatches)
      else numMismatches(bases, primer.bases, math.max(0, bases.length - primer.length), maxMismatches)
    }
    if (mm <= maxMismatches) Some(UngappedAlignmentPrimerMatch(primer, mm, Int.MaxValue)) else None
  }

  /** Counts the # of mismatches, allowing the primer to have IUPAC bases. */
  private[identifyprimers] def numMismatches(bases: Array[Byte], primer: Array[Byte], basesStartIndex: Int, maxMismatches: Double): Int = {
    val length = math.min(bases.length, primer.length)
    var count = 0
    var primerIndex = 0
    while (primerIndex < length && count <= maxMismatches) { // empirically faster than a fgbio forloop
      if (!SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(bases(primerIndex + basesStartIndex), primer(primerIndex))) count += 1
      primerIndex += 1
    }
    count
    //bases.zip(primer).count { case (l, r) => !SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(l.toByte, r.toByte) }
  }
}


/** Finds primer matches based on mapping position.
  * @param primers the primers against which to match
  * @param slop the maximum distance from the primer 5' start and the read's 5' start.
  * @param requireSameStrand require that the read is mapped to the same strand as the primer.
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class LocationBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val requireSameStrand: Boolean,
 val maxMismatchRate: Double) extends _PrimerMatcher with _MismatchAlign with LazyLogging {
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

  /** Gets forward (reverse) primers that overlap the read's start (end) position. The strand of the read is used
    * to determine against which primers to match (i.e. forward or reverse), unless requireSameStrand is false, in
    * which case both strand primers are searched against. */
  private def getOverlaps(rec: SamRecord): Iterator[Primer] = {
    val fwdOverlaps: Iterator[Primer] = if (!requireSameStrand || rec.positiveStrand) {
      val pos = rec.unclippedStart
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop) // TODO: off the end of the refNameosome?
      fwdDetector.getOverlaps(interval).iterator()
        .filter(primer => math.abs(primer.start - pos) <= slop)
    } else Iterator.empty
    val revOverlaps: Iterator[Primer] = if (!requireSameStrand || rec.negativeStrand) {
      val pos = rec.unclippedEnd
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop) // TODO: off the end of the refNameosome?
      revDetector.getOverlaps(interval).iterator()
        .filter(primer => math.abs(primer.end - pos) <= slop)
    } else Iterator.empty
    fwdOverlaps ++ revOverlaps
  }

  /** Returns a location-based match, otherwise None */
  def find(rec: SamRecord): Option[LocationBasedPrimerMatch] = if (rec.unmapped) None else {
    // Find primers that overlap the record, then prioritize by those matching the same strand.  Only one primer may be
    // returned.  In the future, it may be nice to either pick one randomly, keep them all, or  prioritize by # of
    // mismatches.
    val primerMatches = getOverlaps(rec)
      .toSeq match {
      case Seq()                       => None
      case Seq(primer)                 => Some(primer)
      case _primers    if rec.unmapped => unreachable(s"Found multiple primers for $rec\n${_primers.mkString("\n")}")
      case _primers                    =>
        _primers.filter(_.positiveStrand == rec.positiveStrand) match {
          case Seq()       => None
          case Seq(primer) => Some(primer)
          case ssPrimers   => unreachable(s"Found multiple primers on the same strand for $rec\n${ssPrimers.mkString("\n")}")
        }
    }
    primerMatches
      .flatMap { p =>
        // Check the # of mismatches is OK
        mismatchAlign(rec.bases, p).map { m =>
          // Convert to a [[LocationBasedPrimerMatch]]
          LocationBasedPrimerMatch(primer = m.primer, numMismatches = m.numMismatches)
        }
      }
  }
}

/**
  *
  * @param primers the primers against which to match
  * @param slop the extra # of read bases to include when mapping gapped alignment tasks.
  * @param requireSameStrand match against primers on the same strand as the read.
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class UngappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val requireSameStrand: Boolean,
 val maxMismatchRate: Double,
 val kmerLength: Option[Int] = None
 ) extends _PrimerMatcherWithCache with _MismatchAlign with LazyLogging {

  /** Examines all primers to find the best mismatch-based alignment match. */
  def find(rec: SamRecord): Option[UngappedAlignmentPrimerMatch] = {
    val primersToCheck = getPrimersForAlignmentTasks(rec)
    val alignments = if (rec.unmapped) {
      val forward = rec.bases
      val reverse = SequenceUtil.reverseComplement(rec.basesString).getBytes
      val both    = Seq(forward, reverse)

      primersToCheck.flatMap { primer: Primer =>
        both.flatMap { bases => mismatchAlign(bases, primer) }
      }
    }
    else {
      val bases = rec.bases
      primersToCheck.flatMap { primer: Primer =>
        if (!requireSameStrand || rec.positiveStrand == primer.positiveStrand) mismatchAlign(bases, primer)
        else None
      }
    }
    getBestUngappedAlignment(alignments)
  }

  /** Examines all primers to find the best gapped-alignment-based match. */
  private[identifyprimers] def getBestUngappedAlignment(alignments: Seq[UngappedAlignmentPrimerMatch]): Option[UngappedAlignmentPrimerMatch] = {
    if (alignments.isEmpty) None
    else {
      val best: UngappedAlignmentPrimerMatch = alignments.minBy(_.numMismatches)
      alignments.filter(_.numMismatches != best.numMismatches) match {
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
  * @param requireSameStrand match against primers on the same strand as the read.
  * @param aligner the aligner for gapped-alignment-based matching.
  * @param minAlignmentScoreRate the minimum per-base alignment score rate for a gapped alignment based match.  The
  *                              rate for alignment is the score divided by the primer length.
  * @param kmerLength require a read to share a k-mer of this length with a primer before performing gapped alignment
  */
private[identifyprimers] class GappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val requireSameStrand: Boolean,
 val aligner: Aligner,
 val minAlignmentScoreRate: Double,
 val kmerLength: Option[Int] = None) extends _PrimerMatcherWithCache with LazyLogging {

  /** A map from a kmer to the list of primers that contain that kmer. */
  private val kmerToPrimers: Map[ByteBuffer, Set[Primer]] = kmerLength match {
    case None => Map.empty[ByteBuffer, Set[Primer]]
    case Some(len) =>
      val kmers = new mutable.TreeMap[ByteBuffer,mutable.HashSet[Primer]]()
      primers.foreach { primer =>
        val until = primer.length - len
        forloop (from = 0, until = until) { from =>
          val kmer = ByteBuffer.wrap(primer.bases.slice(from = from, until = from + len - 1))
          kmers.get(kmer) match {
            case Some(b) => b.add(primer)
            case None    =>
              val buffer = new mutable.HashSet[Primer]()
              kmers.put(kmer, buffer)
              buffer.add(primer)

          }
        }
      }
      kmers.iterator.map { case (bytes, _primers) => bytes -> _primers.toSet }.toMap
  }

  /** Returns a gapped-alignment-based match, otherwise None */
  def find(rec: SamRecord): Option[GappedAlignmentPrimerMatch] = {
    val alignmentTasks = toAlignmentTasks(rec)

    if (alignmentTasks.isEmpty) None else {
      val alignments = alignmentTasks.map { task =>
        val alignment = this.aligner.align(query = task.queryBytes, target = task.targetBytes)
        AlignmentAndPrimer(alignment, task.query)
      }.toIndexedSeq

      val best: AlignmentAndPrimer = alignments.maxBy(_.alignment.score)
      alignments.filter(_.alignment.score != best.alignment.score) match {
        case Seq() => toGappedAlignmentPrimerMatch(best, None)
        case alns => toGappedAlignmentPrimerMatch(best, Some(alns.maxBy(_.alignment.score).alignment.score))
      }
    }
  }

  // A little clss to store an alignment to a given primer
  private case class AlignmentAndPrimer(alignment: Alignment, primer: Primer)

  /** A little helper method for [[find()]] to create a [[GappedAlignmentPrimerMatch]]. */
  private def toGappedAlignmentPrimerMatch(best: AlignmentAndPrimer, nextBestScore: Option[Int]): Option[GappedAlignmentPrimerMatch] = {
    val primerMatch = GappedAlignmentPrimerMatch(
      primer          = best.primer,
      score           = best.alignment.score,
      secondBestScore = nextBestScore.getOrElse((minAlignmentScoreRate * best.primer.length).toInt)
    )
    Some(primerMatch)
  }
}