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

package com.fulcrumgenomics.bam.pileup

import com.fulcrumgenomics.bam.api.SamRecord
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable

/** Represents a pileup of SAM records and their alignments at a single genomic location. */
case class Pileup[A <: PileupEntry](refName: String, refIndex: Int, pos: Int, pile: Seq[A]) extends Iterable[A] {

  /** Returns the depth of coverage at this position. This includes reads who have a base at this position as well as
    * reads with deletions that span this position, but does NOT include insertions.
    */
  lazy val depth: Int = pile.count {
    case _: InsertionEntry => false
    case _                 => true
  }

  /** An iterator over the elements in this pile. */
  override def iterator: Iterator[A] = pile.iterator

  /** Returns the [[depth]] of the pileup (NOT the number of records in the pileup). */
  override def size: Int = depth

  /** Returns a copy of this pileup with only a single observation from any template. */
  def withoutOverlaps: Pileup[A] = {
    val builder = IndexedSeq.newBuilder[A]
    val names   = mutable.HashSet.empty[String]
    names.sizeHint(pile.length)
    pile.foreach { p => if (names.add(p.rec.name)) builder += p }
    copy(pile = builder.result())
  }

  /** Returns a copy of this pileup object without any indels and only base entries. */
  def withoutIndels: Pileup[BaseEntry] = copy(pile = this.baseIterator.toSeq)

  /** Provides an iterator over the non-indel entries in the pileup. */
  def baseIterator: Iterator[BaseEntry] = this.iterator.collect { case x: BaseEntry => x }
}

/** Base trait for pileup entries that exposes the [[SamRecord]] and the 0-based offset into the record's bases and
  * qualities that is relevant to the pileup.
  */
sealed trait PileupEntry {

  /** The SamRecord that the pileup entry represents. */
  val rec: SamRecord

  /** The zero-based offset within the [[SamRecord]] that is represented.
    *
    *   - For matches and mismatches this is the offset of the base in question.
    *   - For deletions it is the base prior to the deleted sequence.
    *   - For insertions it is the first inserted base in the read.
    */
  val offset: Int

  require(
    offset >= 0 && offset < rec.length,
    s"Offset must be between 0 and read length. Offset: $offset, Read Length: ${rec.length}"
  )
}

/** Pileup entry representing a match or mismatch. */
case class BaseEntry(override val rec: SamRecord, override val offset: Int) extends PileupEntry {
  val base: Byte = SequenceUtil.upperCase(rec.bases(offset))
  val qual: Byte = rec.quals(offset)

  /** Returns the 1-based position within the read in the read's current orientation. */
  def positionInRead: Int = offset + 1

  /** Returns the 1-based position within the read in the order the read was sequenced. This can be thought of as the
    * offset from the 5-prime end plus one.
    */
  def positionInReadInReadOrder: Int = {
    if (rec.negativeStrand) rec.length - offset else positionInRead
  }

  /** Returns the base as it was read on the sequencer. */
  def baseInReadOrientation: Byte = if (rec.negativeStrand) SequenceUtil.complement(base) else base
}

/** Pileup entry representing an insertion. */
case class InsertionEntry(override val rec: SamRecord, override val offset: Int) extends PileupEntry

/** Pileup entry representing a deletion. */
case class DeletionEntry(override val rec: SamRecord, override val offset: Int) extends PileupEntry
