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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools.{CigarOperator, SAMRecord, SAMSequenceDictionary}
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

/**
  * Represents a pileup of reads/bases at a single genomic location.
  */
case class Pileup[A <: PileupEntry](refName: String, refIndex: Int, pos: Int, pile: Seq[A]) extends Iterable[A] {
  override def iterator: Iterator[A] = pile.iterator

  /**
    * Returns the depth of coverage at this position. This includes reads who have a base at
    * this position as well as reads with deletions that span this position, but does NOT
    * include insertions.
    */
  lazy val depth = pile.count {
    case i: InsertionEntry => false
    case _                 => true
  }

  /** Returns the [[depth]] of the pileup (NOT the number of records in the pileup. */
  override def size: Int = depth

  /** Returns a copy of this pileup with only a single observation from any template. */
  def withoutOverlaps: Pileup[A] = {
    val buffer = new ListBuffer[A]()
    val names  = new mutable.HashSet[String]()
    this.pile.foreach { p => if (names.add(p.rec.getReadName)) buffer += p }
    copy(pile=buffer)
  }

  /** Returns a copy of this pileup object without any indels. */
  def withoutIndels: Pileup[BaseEntry] = copy(pile=this.baseIterator.toSeq)

  /** Provides an iterator over the non-indel entries in the pile. */
  def baseIterator: Iterator[BaseEntry] = iterator.collect { case x: BaseEntry => x }
}

/**
  * Base trait for pileup entries that exposes the [[SAMRecord]] and the 0-based offset into
  * the record's bases and qualities that is relevant to the pileup.
  */
sealed trait PileupEntry {
  /** The SAMRecord that the pileup entry represents. */
  val rec: SAMRecord

  /**
    * The zero-based offset within the [[SAMRecord]] that is represented.
    * For matches and mismatches this is the offset of the base in question.
    * For deletions it is the base prior to the deleted sequence.
    * For insertions it is the first inserted base in the read.
    */
  val offset: Int

  require(offset >= 0 && offset < rec.getReadLength, s"Offset must be between 0 and read length. Offset=$offset, rlen=${rec.getReadLength}")
}

/** Pileup entry representing a match or mismatch. */
case class BaseEntry(override val rec: SAMRecord, override val offset: Int) extends PileupEntry {
  val base: Byte = SequenceUtil.upperCase(rec.getReadBases()(offset))
  val qual: Byte = rec.getBaseQualities()(offset)

  /** Returns the 1-based position within the read in the read's current orientation. */
  def positionInRead: Int = offset + 1

  /**
    * Returns the 1-based position within the read in the order the read was sequenced.
    * This can be thought of as the offset from the 5' end plus one.
    */
  def positionInReadInReadOrder: Int = {
    if (rec.getReadNegativeStrandFlag) rec.getReadLength - offset else positionInRead
  }

  /** Returns the base as it was read on the sequencer. */
  def baseInReadOrientation: Byte = if (rec.getReadNegativeStrandFlag) SequenceUtil.complement(base) else base

}

/** Pileup entry representing an insertion. */
case class InsertionEntry(override val rec: SAMRecord, override val offset: Int) extends PileupEntry

/** Pileup entry representing a deletion. */
case class DeletionEntry(override val rec: SAMRecord, override val offset: Int) extends PileupEntry

/**
  * Class that provides methods to build and filter [[Pileup]]s.
  *
  * @param dict the SAMSequenceDictionary for the reference sequence
  * @param minMapQ the minimum MAPQ to allow a read to contribute to the pileup
  * @param minBaseQ the minimum base quality to allow a base to contribute to the pileup
  * @param mappedPairsOnly if true only allow reads from read-pairs with both reads mapped to contribute to the pileup
  * @param includeDuplicates if true allow reads flagged as duplicates to contribute to the pileup
  * @param includeSecondaryAlignments if true allow secondary alignments to contribute to the pileup
  * @param includeSupplementalAlignments if true allow supplemental alignments to contribute to the pileup
  */
class PileupBuilder(dict: SAMSequenceDictionary,
                    minMapQ: Int = 20,
                    minBaseQ: Int = 20,
                    mappedPairsOnly: Boolean = true,
                    includeDuplicates: Boolean = false,
                    includeSecondaryAlignments: Boolean = false,
                    includeSupplementalAlignments: Boolean = false) {

  private val readFilters = new ListBuffer[SAMRecord => Boolean]()
  private val baseFilters = new ListBuffer[PileupEntry => Boolean]()

  // Setup the default read level filters
  if (minMapQ > 0)                    readFilters += { r => r.getMappingQuality >= minMapQ }
  if (mappedPairsOnly)                readFilters += { r => r.getReadPairedFlag && !r.getReadUnmappedFlag && !r.getMateUnmappedFlag }
  if (!includeDuplicates)             readFilters += { r => !r.getDuplicateReadFlag }
  if (!includeSecondaryAlignments)    readFilters += { r => !r.getNotPrimaryAlignmentFlag }
  if (!includeSupplementalAlignments) readFilters += { r => !r.getSupplementaryAlignmentFlag }

  // Setup the default base level filters
  if (minBaseQ > 0) baseFilters += {
    case p: BaseEntry => p.qual >= minBaseQ
    case _ => true
  }

  /** Adds a read level filter to the set of filters. */
  def withReadFilter(f: SAMRecord => Boolean): PileupBuilder = yieldAndThen(this) { this.readFilters += f }

  /** Adds a base level filter to the set of filters. */
  def withBaseFilter(f: PileupEntry => Boolean): PileupBuilder = yieldAndThen(this) { this.baseFilters += f }

  /** Checks to see if all read filters accept the read. */
  protected def accepts(rec: SAMRecord): Boolean = this.readFilters.forall(f => f(rec))

  /** Checks to see if all read filters accept the read. */
  protected def accepts(entry: PileupEntry): Boolean = this.baseFilters.forall(f => f(entry))

  /**
    * Constructs a pileup, at the single requested position, from the records returned by the iterator.
    *
    * @param iterator an Iterator of SAMRecords that may or may not overlap the position
    * @param refName the name of the reference sequence/contig/chrom on which the position resides
    * @param pos the 1-based position on the reference sequence at which to construct the pileup
    */
  def build(iterator: Iterator[SAMRecord], refName: String, pos: Int): Pileup[PileupEntry] = {
    val refIndex = dict.getSequenceIndex(refName)
    require(refIndex >= 0, s"Unknown reference sequence: ${refName}")

    val buffer = new ListBuffer[PileupEntry]()
    iterator.bufferBetter
      .dropWhile(r => r.getReferenceIndex < refIndex)
      .takeWhile(r => r.getReferenceIndex == refIndex && r.getAlignmentStart <= pos + 1)
      .filterNot(r => r.getReadUnmappedFlag)
      .filter(r => r.getAlignmentStart <= pos || startsWithInsertion(r))
      .filter(r => r.getAlignmentEnd >= pos)
      .filter(r => this.readFilters.forall(f => f(r)))
      .foreach { rec =>
        if (rec.getAlignmentStart == pos + 1) {
          // Must be an insertion at the start of the read
          val entry = new InsertionEntry(rec, 0)
          if (accepts(entry)) buffer += entry
        }
        else {
          rec.getReadPositionAtReferencePosition(pos) match {
            case 0 =>
              // Must be deleted in the read
              val delPos = rec.getReadPositionAtReferencePosition(pos, true)
              val entry = DeletionEntry(rec, delPos-1)
              if (accepts(entry)) buffer += entry
            case p =>
              val entry = BaseEntry(rec, p-1)
              if (accepts(entry)) buffer += entry

              // Also check to see if the subsequent base represents an insertion
              if (p < rec.getReadLength - 1 && rec.getReferencePositionAtReadPosition(p+1) == 0) {
                val entry = InsertionEntry(rec, p)
                if (accepts(entry)) buffer += entry
              }
          }
        }
      }

    Pileup(refName=refName, refIndex=refIndex, pos=pos, pile=buffer)
  }

  /** Returns true if the read is mapped and the first non-clipping operator is an insertion. */
  private def startsWithInsertion(r: SAMRecord) = {
    if (r.getReadUnmappedFlag) {
      false
    }
    else {
      val cigarOps = r.getCigar.getCigarElements.iterator().bufferBetter.dropWhile(_.getOperator.isClipping)
      cigarOps.hasNext && cigarOps.head.getOperator == CigarOperator.INSERTION
    }
  }
}

