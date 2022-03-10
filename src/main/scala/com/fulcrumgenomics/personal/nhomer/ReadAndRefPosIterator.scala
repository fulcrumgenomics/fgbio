/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

package com.fulcrumgenomics.personal.nhomer

import com.fulcrumgenomics.FgBioDef.yieldAndThen
import com.fulcrumgenomics.alignment.CigarElem
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.personal.nhomer.ReadAndRefPosIterator.ReadAndRefPos
import htsjdk.samtools.util.CoordMath


object ReadAndRefPosIterator {
  /** Stores a 1-based position in the read and reference
    *
    * @param read 1-based position in the read
    * @param ref 1-based position in the reference
    */
  case class ReadAndRefPos(read: Int, ref: Int)

  /** Builds a [[ReadAndRefPosIterator]] over a [[SamRecord]].
    *
    * The start/end read/rec positions default to the full alignment unless given.
    * */
  def apply(rec: SamRecord,
            startReadPos: Option[Int] = None,
            endReadPos: Option[Int] = None,
            startRefPos: Option[Int] = None,
            endRefPos: Option[Int] = None,
           ): ReadAndRefPosIterator = new ReadAndRefPosIterator(
    rec          = rec,
    startReadPos = startReadPos.getOrElse(1),
    endReadPos   = endReadPos.getOrElse(rec.length),
    startRefPos  = startRefPos.getOrElse(rec.start),
    endRefPos    = endRefPos.getOrElse(rec.end),
  )
}

/** An iterator for each _mapped_ read base (i.e. has a corresponding reference base), with each value being the
  * 1-based position in the read and reference respectively.
  *
  * Insertions and deletions will not be returned.  For example, an insertion consumes a read base, but has no
  * corresponding reference base.  A deletion consumes a reference base, but has no corresponding read base.
  *
  * @param rec the record over which read and reference positions should be returned
  * @param startReadPos the 1-based start position in the read
  * @param endReadPos the 1-based inclusive end position in the read
  * @param startRefPos the 1-based start position in the reference
  * @param endRefPos the 1-based inclusive end position in the reference
  * */
class ReadAndRefPosIterator (rec: SamRecord,
                             startReadPos: Int,
                             endReadPos: Int,
                             startRefPos: Int,
                             endRefPos: Int,
                            ) extends Iterator[ReadAndRefPos] {
  private val elems        = rec.cigar.elems // the cigar elements
  private var curRefPos    = rec.start // the current 1-based position in the reference
  private var curReadPos   = 1 // the current 1-based position in the read
  private var elementIndex = 0 // the current element index (0-based)
  private var inElemOffset = 0 // the current offset in the current element (0-based)

  // the next read and reference position
  private var nextValue: Option[ReadAndRefPos] = None

  // Initialize: iterate through the cigar elements to find the first aligned (M/X/EQ) element that
  // that has bases that map at or past both the start of the ref and read window respectively
  {
    // The current 1-based inclusive reference end
    def curRefEnd: Int  = CoordMath.getEnd(curRefPos, curElem.lengthOnTarget)
    // The current 1-based inclusive read end
    def curReadEnd: Int = CoordMath.getEnd(curReadPos, curElem.lengthOnQuery)

    // skip until we have an element is at or past both the read/ref starts
    while (elementIndex < elems.length
      && (curRefEnd < startRefPos || curReadEnd < startReadPos)
    ) {
      curRefPos += curElem.lengthOnTarget
      curReadPos += curElem.lengthOnQuery
      elementIndex += 1
    }

    // skip to the next mapped element
    skipNonAlignedBasesAndAdjustOffset()

    // initialize the next value to return
    this.nextValue = yieldNextValue()
  }

  override def hasNext: Boolean = nextValue.nonEmpty

  override def next(): ReadAndRefPos = {
    yieldAndThen(nextValue.get) {
      nextValue = this.yieldNextValue()
    }
  }

  /** Gets the current element; no OOB is performed */
  @inline private def curElem: CigarElem = this.elems(this.elementIndex)

  /** Gets the next value to return, or None if no more elements exist that fall within the read/ref window */
  private def yieldNextValue(): Option[ReadAndRefPos] = {
    // if we've consumed all the mapped bases in the current element, move the the next mapped element
    if (this.elementIndex < this.elems.length && curElem.length <= inElemOffset) {
      elementIndex += 1
      skipNonAlignedBasesAndAdjustOffset()
    }

    if (this.elementIndex == this.elems.length || endReadPos < curReadPos || endRefPos < curRefPos) None
    else {
      // return the current value and increment the element offset
      yieldAndThen(Some(ReadAndRefPos(read=curReadPos, ref=curRefPos))) {
        curReadPos += 1
        curRefPos += 1
        inElemOffset += 1
      }
    }
  }

  /** Consumes non-aligned (non-M/X/EQ) elements, then adjusts the read/ref position and offset into the next mapped
    * element if necessary. */
  private def skipNonAlignedBasesAndAdjustOffset(): Unit = {
    inElemOffset = 0

    // skip over any non-aligned bases
    while (elementIndex < elems.length && !curElem.operator.isAlignment) {
      curRefPos += curElem.lengthOnTarget
      curReadPos += curElem.lengthOnQuery
      elementIndex += 1
    }

    // update the offset in the current element only if the current element spans the read-start or reference-start
    if (elementIndex < elems.length && (curRefPos < startRefPos || curReadPos < startReadPos)) {
      inElemOffset = Math.max(startRefPos - curRefPos, startReadPos - curReadPos)
      curRefPos += inElemOffset
      curReadPos += inElemOffset
    }
  }
}

object MateOverlappingReadAndRefPosIterator {
  /** Builds a [[MateOverlappingReadAndRefPosIterator]] over a [[SamRecord]] and its mate.
    *
    * The start/end read positions default to the full overlap unless given.
    * */
  def apply(rec: SamRecord,
            mate: SamRecord,
            startReadPos: Option[Int] = None,
            endReadPos: Option[Int] = None
           ): MateOverlappingReadAndRefPosIterator = new MateOverlappingReadAndRefPosIterator(
    rec          = rec,
    mate         = mate,
    startReadPos = startReadPos.getOrElse(1),
    endReadPos   = endReadPos.getOrElse(rec.length),
  )
}


/** An iterator for each _mapped_ read base (i.e. has a corresponding reference base) that overlaps the alignment of
  * its mate, with each value being the 1-based position in the read and reference respectively.
  *
  * Insertions and deletions will not be returned.  For example, an insertion consumes a read base, but has no
  * corresponding reference base.  A deletion consumes a reference base, but has no corresponding read base.
  *
  * @param rec the record to iterate
  * @param mate the mate of the record
  */
class MateOverlappingReadAndRefPosIterator(rec: SamRecord,
                                           mate: SamRecord,
                                           startReadPos: Int,
                                           endReadPos: Int)
  extends ReadAndRefPosIterator(
    rec          = rec,
    startReadPos = startReadPos,
    endReadPos   = endReadPos,
    startRefPos  = Math.max(rec.start, mate.start),
    endRefPos    = Math.min(rec.end, mate.end)
  ) {
  require(rec.paired && rec.mapped && mate.paired && mate.mapped && rec.name == mate.name)
  require(rec.refIndex == mate.refIndex)
}
