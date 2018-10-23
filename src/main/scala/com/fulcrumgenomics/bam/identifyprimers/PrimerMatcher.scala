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

private[identifyprimers] object PrimerMatcher {
  /** Gets the bases of the record in sequencing order. */
  implicit class SequencingOrderBases(rec: SamRecord) {
    def basesInSequencingOrder: Array[Byte] = {
      if (rec.unmapped || rec.positiveStrand) rec.bases else {
        SequenceUtil.reverseComplement(rec.basesString).getBytes
      }
    }
  }
}

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
  import PrimerMatcher.SequencingOrderBases

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
      val kmers = new mutable.TreeMap[ByteBuffer, mutable.HashSet[Primer]]()
      primers.foreach { primer =>
        // create all possible kmers from the primer
        val until = primer.length - len
        forloop (from = 0, until = until) { from =>
          // extract the kmer
          val kmer = ByteBuffer.wrap(primer.basesInSequencingOrder.slice(from = from, until = from + len))
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

  /** Gets the list of primers for ungapped and gapped alignments.
    *
    * If `kmerLength` is defined, uses the k-mer hash to filter the returned primers for common kmers at the start of
    * the read (in sequencing order) and primers.
    */
  protected[identifyprimers] def getPrimersForAlignmentTasks(rec: SamRecord): Seq[Primer] = {
    this.kmerLength match {
      case None             => this.primers
      case Some(kmerLength) => getPrimersForAlignmentTasksWithCommonKmer(rec.basesInSequencingOrder, kmerLength)
    }
  }

  /** Gets the list of primers for ungapped and gapped alignments.

    */
  private[identifyprimers] def getPrimersForAlignmentTasksWithCommonKmer(bases: Array[Byte], kmerLength: Int): Seq[Primer] = {
    if (bases.length < kmerLength) Seq.empty else {
      val end = math.min(maxPrimerLength, bases.length) - kmerLength
      Range.inclusive(start = 0, end = end).flatMap { from =>
        val until = from + kmerLength
        val kmer  = bases.slice(from = from, until = until)
        this.kmerToPrimers.getOrElse(ByteBuffer.wrap(kmer), Seq.empty)
      }.distinct
    }
  }
}


/** A little trait to facilitate ungapped alignment. */
private trait UngappedAlignment {
  import PrimerMatcher.SequencingOrderBases

  /** The maximum per-base mismatch rate allowed for a [[PrimerMatch]]. */
  def maxMismatchRate: Double

  /** Returns a mismatch-based alignment match, otherwise None */
  protected[identifyprimers] def mismatchAlign(rec: SamRecord, primer: Primer): Option[UngappedAlignmentPrimerMatch] = {
    val bases          = rec.basesInSequencingOrder
    val maxMismatches  = math.min(bases.length, primer.length) * this.maxMismatchRate
    val mm             = numMismatches(bases, primer.basesInSequencingOrder, 0, maxMismatches)
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
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class LocationBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val maxMismatchRate: Double) extends PrimerMatcher with UngappedAlignment with LazyLogging {
  import com.fulcrumgenomics.FgBioDef.javaIteratorAsScalaIterator

  private val positiveStrandDetector: OverlapDetector[Primer] = {
    val d = new OverlapDetector[Primer](0, 0)
    this.primers.filter(p => p.mapped && p.positiveStrand).foreach { primer => d.addLhs(primer, primer) }
    d
  }

  private val negativeStrandDetector: OverlapDetector[Primer] = {
    val d = new OverlapDetector[Primer](0, 0)
    this.primers.filter(p => p.mapped && p.negativeStrand).foreach { primer => d.addLhs(primer, primer) }
    d
  }

  /** Gets forward (reverse) primers that overlap the read's start (end) position. The strand of the read is used
    * to determine against which primers to match (i.e. forward or reverse). Forward strand primers are searched for at
    * the start of the read, while reverse strand primers are searched for at the end of the read. */
  private[identifyprimers] def getOverlaps(rec: SamRecord): Iterator[Primer] = {
    require(rec.mapped)
    if (rec.positiveStrand) {
      val pos      = rec.unclippedStart
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop)
      positiveStrandDetector.getOverlaps(interval).iterator()
        .filter(primer => math.abs(primer.start - pos) <= slop)
    }
    else {
      val pos      = rec.unclippedEnd
      val interval = new Interval(rec.refName, math.max(1, pos-slop), pos+slop)
      negativeStrandDetector.getOverlaps(interval).iterator()
        .filter(primer => math.abs(primer.end - pos) <= slop)
    }
  }

  /** Returns a [[LocationBasedPrimerMatch]], otherwise [[None]] */
  def find(rec: SamRecord): Option[LocationBasedPrimerMatch] = if (rec.unmapped) None else {
    // Find primers that overlap the record, then prioritize by those matching the same strand.  At most one primer will
    // be returned.  In the future, it may be nice to either pick one randomly, keep them all, or prioritize by # of
    // mismatches.
    val primerMatches: Option[Primer] = getOverlaps(rec).toSeq match {
      case Seq()                       => None
      case Seq(primer)                 => Some(primer)
      case _primers                    =>
        _primers.filter(_.positiveStrand == rec.positiveStrand) match {
          case Seq()       => None
          case Seq(primer) => Some(primer)
          case ssPrimers   => throw new IllegalStateException(s"Found multiple primers on the same strand for $rec\n\t${ssPrimers.mkString("\n\t")}")
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
  * @param maxMismatchRate the maximum per-base mismatch rate to accept a primer match.
  */
private[identifyprimers] class UngappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
 val maxMismatchRate: Double,
 val kmerLength: Option[Int] = None
 ) extends PrimerMatcherWithKmerFilter with UngappedAlignment with LazyLogging {

  /** Examines all primers to find the best mismatch-based alignment match. */
  def find(rec: SamRecord): Option[UngappedAlignmentPrimerMatch] = {
    val alignments = getPrimersForAlignmentTasks(rec).flatMap { p => mismatchAlign(rec, p) }
    getBestUngappedAlignment(alignments)
  }

  /** Examines all primers to find the best gapped-alignment-based match. */
  private[identifyprimers] def getBestUngappedAlignment(alignments: Seq[UngappedAlignmentPrimerMatch]): Option[UngappedAlignmentPrimerMatch] = {
    if (alignments.isEmpty) None
    else {
      alignments.sortBy(_.numMismatches).take(2) match {
        case Seq()               => unreachable("Should have found at least one.")
        case Seq(best)           => Some(best.copy(nextNumMismatches = Int.MaxValue))
        case Seq(best, nextBest) => Some(best.copy(nextNumMismatches=nextBest.numMismatches))
      }
    }
  }
}

/** A little class that wraps a [[SamRecord]] to make it [[Alignable]] and also to defer extracting the bases to align.
  * NB: end is exclusive, and start and end are zero-based. */
private[identifyprimers] case class SamRecordAlignable(rec: SamRecord, override val length: Int) extends Alignable {
  import PrimerMatcher.SequencingOrderBases
  def bases: Array[Byte] = rec.basesInSequencingOrder.slice(0, length)
}

/** A primer matcher.  Attempts to match based on location, then based on a mismatch-only alignment, then based on a
  * gapped-alignment.
  *
  * @param primers the primers against which to match
  * @param slop the extra # of read bases to include when mapping gapped alignment tasks.
  * @param aligner the aligner for gapped-alignment-based matching.
  * @param minAlignmentScoreRate the minimum per-base alignment score rate for a gapped alignment based match.  The
  *                              rate for alignment is the score divided by the primer length.
  * @param kmerLength require a read to share a k-mer of this length with a primer before performing gapped alignment
  */
private[identifyprimers] class GappedAlignmentBasedPrimerMatcher
(val primers: Seq[Primer],
 val slop: Int,
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
    val samRecordAlignableCache = scala.collection.mutable.HashMap[Int, SamRecordAlignable]()
    val primersToCheck = getPrimersForAlignmentTasks(rec)
    primersToCheck.map { primer =>
      val targetLength = math.min(primer.sequence.length + slop, rec.length)
      val target       = samRecordAlignableCache.getOrElseUpdate(targetLength, SamRecordAlignable(rec, targetLength))
      AlignmentTask(query = primer, target = target)
    }
  }

  // A little class to store an alignment to a given primer
  private[identifyprimers] case class AlignmentAndPrimer(alignment: Alignment, primer: Primer)

  /** A little helper method for [[find()]] to create a [[GappedAlignmentPrimerMatch]]. */
  private[identifyprimers] def getBestGappedAlignment(alignments: Seq[AlignmentAndPrimer]): Option[GappedAlignmentPrimerMatch] = {
    if (alignments.isEmpty) None
    else {
      alignments.sortBy(-_.alignment.score).take(2) match {
        case Seq()                =>
          unreachable("Should have found at least one.")
        case Seq(best)            =>
          val secondBestScore = (minAlignmentScoreRate * best.primer.length).toInt
          if (best.alignment.score < secondBestScore) None else {
            Some(GappedAlignmentPrimerMatch(best.primer, score = best.alignment.score, secondBestScore = (minAlignmentScoreRate * best.primer.length).toInt))
          }
         case Seq(best, nextBest) =>
           Some(GappedAlignmentPrimerMatch(best.primer, score=best.alignment.score, secondBestScore=nextBest.alignment.score))
      }
    }
  }
}