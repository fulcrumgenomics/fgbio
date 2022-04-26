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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.yieldAndThen
import com.fulcrumgenomics.alignment.CigarElem
import com.fulcrumgenomics.bam.ReadAndRefPosIterator.ReadAndRefPos
import com.fulcrumgenomics.bam.ReadMateAndRefPosIterator.ReadMateAndRefPos
import com.fulcrumgenomics.bam.api.SamRecord
import htsjdk.samtools.util.CoordMath

/** An iterator for each _mapped_ read base (i.e. has a corresponding reference base), with each value being the
  * 1-based position in the read and reference respectively.
  *
  * Insertions and deletions will not be returned.  For example, an insertion consumes a read base, but has no
  * corresponding reference base.  A deletion consumes a reference base, but has no corresponding read base.
  *
  * The _mapped_ read bases returned can be limited using the arguments for the minimum and maximum position in the read
  * and reference respectively.
  *
  * @param rec the record over which read and reference positions should be returned
  * @param minReadPos the minimum 1-based position in the read to return
  * @param maxReadPos the maximum 1-based inclusive end position in the read to return
  * @param minRefPos the minimum 1-based start position in the reference to return
  * @param maxRefPos the maximum 1-based inclusive end position in the reference to return
  * */
class ReadAndRefPosIterator (rec: SamRecord,
                             minReadPos: Int = 1,
                             maxReadPos: Int = Int.MaxValue,
                             minRefPos: Int = 1,
                             maxRefPos: Int = Int.MaxValue,
                            ) extends Iterator[ReadAndRefPos] {
  private val elems        = rec.cigar.elems // the cigar elements
  private var curRefPos    = rec.start // the current 1-based position in the reference
  private var curReadPos   = 1 // the current 1-based position in the read
  private var elementIndex = 0 // the current element index (0-based)
  private var inElemOffset = 0 // the current offset in the current element (0-based)

  // the next read and reference position
  private var nextValue: Option[ReadAndRefPos] = None

  // the final bounds in the read and reference to return
  private val startReadPos = math.max(1, minReadPos)
  private val endReadPos   = math.min(rec.length, maxReadPos)
  private val startRefPos  = math.max(rec.start, minRefPos)
  private val endRefPos    = math.min(rec.end, maxRefPos)

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
      yieldAndThen(Some(ReadAndRefPos(readPos=curReadPos, refPos=curRefPos))) {
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

object ReadAndRefPosIterator {
  /** Stores a 1-based position in the read and reference
    *
    * @param readPos 1-based position in the read
    * @param refPos 1-based position in the reference
    */
  case class ReadAndRefPos(readPos: Int, refPos: Int)

  /** An iterator for each _mapped_ read base (i.e. has a corresponding reference base) that overlaps the alignment of
    * its mate, with each value being the 1-based position in the read and reference respectively.
    *
    * Insertions and deletions will not be returned.  For example, an insertion consumes a read base, but has no
    * corresponding reference base.  A deletion consumes a reference base, but has no corresponding read base.
    *
    * The _mapped_ read bases returned can be limited using the arguments for the minimum and maximum position in the read.
    *
    * @param rec the record over which read and reference positions should be returned
    * @param mate the mate of the record
    * @param minReadPos the minimum 1-based position in the read to return
    * @param maxReadPos the maximum 1-based inclusive end position in the read to return
    */
  def apply(rec: SamRecord,
            mate: SamRecord,
            minReadPos: Int = 1,
            maxReadPos: Int = Int.MaxValue): ReadAndRefPosIterator = {
    require(rec.paired && rec.mapped && mate.paired && mate.mapped && rec.name == mate.name)
    require(rec.refIndex == mate.refIndex)

    new ReadAndRefPosIterator(
      rec        = rec,
      minReadPos = minReadPos,
      maxReadPos = maxReadPos,
      minRefPos  = Math.max(rec.start, mate.start),
      maxRefPos  = Math.min(rec.end, mate.end)
    )
  }
}


/** An iterator for each _mapped_ read base (i.e. has a corresponding reference base) that also has a _mapped_ base in
  * its mate, with each value being the 1-based position in the read and reference respectively.
  *
  * Insertions and deletions will not be returned.  For example, an insertion consumes a read base, but has no
  * corresponding reference base.  A deletion consumes a reference base, but has no corresponding read base.
  *
  * The _mapped_ read bases returned can be limited using the arguments for the minimum and maximum position in the read.
  *
  * @param rec the record over which read and reference positions should be returned
  * @param mate the mate of the record
  * @param minReadPos the minimum 1-based position in the read to return
  * @param maxReadPos the maximum 1-based inclusive end position in the read to return
  */
class ReadMateAndRefPosIterator(rec: SamRecord,
                                mate: SamRecord,
                                minReadPos: Int = 1,
                                maxReadPos: Int = Int.MaxValue) extends Iterator[ReadMateAndRefPos] {
  // Use an iterators that returns the position for the read and mate respectively, then just find the reference positions
  // where both have a read position defined.
  private val recIter  = ReadAndRefPosIterator(rec=rec, mate=mate, minReadPos=minReadPos, maxReadPos=maxReadPos).buffered
  private val mateIter = ReadAndRefPosIterator(rec=mate, mate=rec).buffered

  private var nextItem: Option[ReadMateAndRefPos] = None

  override def hasNext(): Boolean = {
    while (nextItem.isEmpty && recIter.hasNext && mateIter.hasNext) {
      val nextRec  = recIter.head
      val nextMate = mateIter.head
      if (nextRec.refPos < nextMate.refPos) recIter.next()
      else if (nextMate.refPos < nextRec.refPos) mateIter.next()
      else {
        nextItem = Some(ReadMateAndRefPos(readPos=nextRec.readPos, matePos=nextMate.readPos, refPos=nextRec.refPos))
        recIter.next()
        mateIter.next()
      }
    }
    nextItem.isDefined
  }

  def next(): ReadMateAndRefPos = {
   if (!hasNext()) throw new NoSuchElementException()
   val retval = this.nextItem.get
   this.nextItem = None
   retval
  }
}

object ReadMateAndRefPosIterator {
  /** Stores a 1-based position in the read, mate and reference
    *
    * @param readPos 1-based position in the read
    * @param matePos 1-based position in the mate
    * @param refPos 1-based position in the reference
    */
  case class ReadMateAndRefPos(readPos: Int, matePos: Int, refPos: Int)
}
