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
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Cigar, CigarElem}
import com.fulcrumgenomics.bam.api.SamRecord
import enumeratum.EnumEntry
import htsjdk.samtools.{SAMUtils, CigarOperator => Op}

import scala.collection.mutable.ArrayBuffer

/** The base trait for all clipping modes. */
sealed trait ClippingMode extends EnumEntry

/** An enumeration representing the various ways to clip bases within a read. */
object ClippingMode extends FgBioEnum[ClippingMode] {
  def values: scala.collection.immutable.IndexedSeq[ClippingMode] = findValues
  /** Soft-clip the read. */
  case object Soft extends ClippingMode
  /** Soft-clip the read and mask the bases and qualities. */
  case object SoftWithMask extends ClippingMode
  /** Hard-clip the read. */
  case object Hard extends ClippingMode
}

/**
  * Holds types and constants related to SamRecord clipping.
  */
object SamRecordClipper {

  /** The set of tags that should be invalidated if a read undergoes clipping. */
  val TagsToInvalidate: Seq[String] = Bams.AlignmentTags

  private val NoCallBase = 'N'.toByte
  private val NoCallQual = 2.toByte
}


/**
  * Provides a suite of methods for clipping (soft- and hard-) bases from the beginnings
  * and ends of reads in various ways.  Cigar strings, bases, qualities and alignment
  * positions are all correctly adjusted post-clipping.
  *
  * Note that there are several "flavours" of method:
  *   - Ones that work on the start and end of the read in whatever orientation the read
  *     is in, vs. working on the 5' or 3' end of the read
  *   - Ones that attempt to clip an additional N bases beyond any clipping already
  *     provided (clip[Start|End|5PrimeEnd|3PrimeEnd]OfAlignment) and ones that
  *     attempt to make it so that the read has N bases clipped, including any existing
  *     clipping (clip[Start|End|5PrimeEnd|3PrimeEnd]OfRead)
  *
  * @param mode how should clipping be performed (hard, soft, soft with masking)
  * @param autoClipAttributes if true attributes that are the same length as the bases
  *                           and qualities be automatically clipped in the same way as
  *                           the bases and qualities, otherwise attributes are not
  *                           touched.
  */
class SamRecordClipper(val mode: ClippingMode, val autoClipAttributes: Boolean) {
  import SamRecordClipper._

  /**
    * Adds clipping of _at least_ numberOfBasesToClip to the start (left hand end)
    * of an alignment. If clipping already exists at the start of the read it is
    * preserved, and numberOfBasesToClip more clipping is added.
    *
    * If is unmapped or the numberOfBasesToClip is < 1 nothing is done.
    *
    * If the read has fewer clippable bases than requested clipping, the read is
    * unmapped.
    *
    * If hard-clipping is requested and the read has existing soft-clipping at the start
    * it is converted to hard-clipping before adding further clipping.
    *
    * If soft-clipping-with-masking is requested and the read already contains soft
    * clipping at the start of the read, both the existing and new soft clipped bases
    * are masked.
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to clip beyond any clipping
    *                            already applied to the read
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clipStartOfAlignment(rec: SamRecord, numberOfBasesToClip: Int) : Int = {
    val numClippable = numberOfClippableBases(rec)
    if (rec.unmapped || numberOfBasesToClip < 1) {
      0
    }
    else if (numClippable <= numberOfBasesToClip) {
      SAMUtils.makeReadUnmapped(rec.asSam)
      numClippable
    }
    else {
      val oldElems = rec.cigar.elems
      val oldBases = rec.bases
      val oldQuals = rec.quals
      val (newElems, readBasesClipped, refBasesClipped, bases, quals) = clip(oldElems, numberOfBasesToClip, oldBases, oldQuals)
      rec.start = rec.start + refBasesClipped
      rec.cigar = Cigar(newElems)
      rec.bases = bases
      rec.quals = quals
      cleanupClippedRecord(rec)
      clipExtendedAttributes(rec, oldBases.length - rec.bases.length, fromStart=true)
      readBasesClipped
    }
  }

  /**
    * Adds clipping of _at least_ numberOfBasesToClip to the end (right hand end)
    * of an alignment. If clipping already exists at the end of the read it is
    * preserved, and numberOfBasesToClip more clipping is added.
    *
    * If is unmapped or the numberOfBasesToClip is < 1 nothing is done.
    *
    * If the read has fewer clippable bases than requested clipping, the read is
    * unmapped.
    *
    * If hard-clipping is requested and the read has existing soft-clipping at the end
    * it is converted to hard-clipping before adding further clipping.
    *
    * If soft-clipping-with-masking is requested and the read already contains soft
    * clipping at the end of the read, both the existing and new soft clipped bases
    * are masked.
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to clip beyond any clipping
    *                            already applied to the read
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clipEndOfAlignment(rec: SamRecord, numberOfBasesToClip: Int) : Int = {
    val numClippable = numberOfClippableBases(rec)
    if (rec.unmapped || numberOfBasesToClip < 1) {
      0
    }
    else if (numClippable <= numberOfBasesToClip) {
      SAMUtils.makeReadUnmapped(rec.asSam)
      numClippable
    }
    else {
      val oldElems = rec.cigar.elems.reverse
      val oldBases = rec.bases.reverse
      val oldQuals = rec.quals.reverse
      val (newElems, readBasesClipped, _, bases, quals) = clip(oldElems, numberOfBasesToClip, oldBases, oldQuals)
      rec.cigar = Cigar(newElems.reverse)
      rec.bases = bases.reverse
      rec.quals = quals.reverse
      cleanupClippedRecord(rec)
      clipExtendedAttributes(rec, oldBases.length - rec.bases.length, fromStart=false)
      readBasesClipped
    }
  }

  /**
    * Attempts to clip an additional numberOfBasesToClip from the 5' end of the read. For
    * details see [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]] and
    * [[com.fulcrumgenomics.bam.SamRecordClipper.clipEndOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clip5PrimeEndOfAlignment(rec: SamRecord, numberOfBasesToClip: Int): Int = {
    if (rec.negativeStrand) clipEndOfAlignment(rec, numberOfBasesToClip)
    else clipStartOfAlignment(rec, numberOfBasesToClip)
  }

  /**
    * Attempts to clip an additional numberOfBasesToClip from the 3' end of the read. For
    * details see [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]] and
    * [[com.fulcrumgenomics.bam.SamRecordClipper.clipEndOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clip3PrimeEndOfAlignment(rec: SamRecord, numberOfBasesToClip: Int): Int = {
    if (rec.negativeStrand) clipStartOfAlignment(rec, numberOfBasesToClip)
    else clipEndOfAlignment(rec, numberOfBasesToClip)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the start (left-hand end)
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param clipLength the total amount of clipping desired, including any existing clipping
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clipStartOfRead(rec: SamRecord, clipLength: Int): Int = {
    val existingClipping = rec.cigar.takeWhile(_.operator.isClipping).map(_.length).sum
    if (clipLength > existingClipping) clipStartOfAlignment(rec, clipLength - existingClipping)
    else { upgradeClipping(rec, clipLength, fromStart=true); 0 }
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the end (right-hand end)
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param clipLength the total amount of clipping desired, including any existing clipping
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clipEndOfRead(rec: SamRecord, clipLength: Int): Int = {
    val existingClipping = rec.cigar.reverseIterator.takeWhile(_.operator.isClipping).map(_.length).sum
    if (clipLength > existingClipping) clipEndOfAlignment(rec, clipLength - existingClipping)
    else { upgradeClipping(rec, clipLength, fromStart=false); 0 }
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the 5' end
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clip5PrimeEndOfRead(rec: SamRecord, numberOfBasesToClip: Int): Int = {
    if (rec.mapped && rec.negativeStrand) clipEndOfRead(rec, numberOfBasesToClip)
    else clipStartOfRead(rec, numberOfBasesToClip)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the 3' end
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfAlignment]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    * @return the additional number bases clipped, not including bases already clipped
    */
  def clip3PrimeEndOfRead(rec: SamRecord, numberOfBasesToClip: Int): Int = {
    if (rec.mapped && rec.negativeStrand) clipStartOfRead(rec, numberOfBasesToClip)
    else clipEndOfRead(rec, numberOfBasesToClip)
  }

  /** Converts all clipping to the current mode.
    * Can convert:
    * 1. From [[ClippingMode.Soft]] to [[ClippingMode.SoftWithMask]]
    * 2. From [[ClippingMode.Soft]] to [[ClippingMode.Hard]]
    * 3. From [[ClippingMode.SoftWithMask]] to [[ClippingMode.Hard]]
    * In all other cases, clipping remains the same.
    *
    * Calculates any clipping required and delegates to [[com.fulcrumgenomics.bam.SamRecordClipper.clipStartOfRead]] and
    * [[com.fulcrumgenomics.bam.SamRecordClipper.clipEndOfRead]].
    * @param rec the record to be clipped
    * @return the number of bases converted at the start and end of the read respectively.  For [[ClippingMode.SoftWithMask]],
    *         any existing masked bases that would be converted will be counted.
    */
  def upgradeAllClipping(rec: SamRecord): (Int, Int) = {
    if (rec.unmapped || this.mode == ClippingMode.Soft) (0, 0) else {
      val numBasesClippedStart = {
        val clippingElems = rec.cigar.elems.takeWhile(_.operator.isClipping)
        val numSoft = clippingElems.filter(_.operator == Op.S).map(_.length).sum
        if (numSoft > 0) this.clipStartOfRead(rec, clippingElems.map(_.length).sum)
        numSoft
      }
      val numBasesClippedEnd = {
        val clippingElems = rec.cigar.reverseIterator.takeWhile(_.operator.isClipping).toSeq
        val numSoft       = clippingElems.filter(_.operator == Op.S).map(_.length).sum
        if (numSoft > 0) this.clipEndOfRead(rec, clippingElems.map(_.length).sum)
        numSoft
      }
      (numBasesClippedStart, numBasesClippedEnd)
    }
  }

  /** Clips overlapping read pairs, where both ends of the read pair are mapped to the same chromosome, and in FR orientation.
    *
    * If the reads do not overlap, or are not an FR pair, return (0, 0).
    *
    * Clips at the reference midpoint between the two reads.
    *
    * @param rec the read
    * @param mate the mate
    * @return the number of overlapping bases to that were clipped on the record and mate respectively (3' end in sequencing order)
    */
  def clipOverlappingReads(rec: SamRecord, mate: SamRecord): (Int, Int) = {
    if (rec.matesOverlap.contains(false)) (0, 0) // do not overlap, don't clip
    else if (!rec.isFrPair) (0, 0)
    else if (rec.negativeStrand) clipOverlappingReads(rec=mate, mate=rec).swap // don't for get to swap the results since we swapped inputs
    else {
      // Pick the mid point in the reference window between the 5' ends of the record and its mate.  The read ends at the
      // mid point, or if the mid point is in a deletion, the base prior to the deletion.  The mate ends at the mid
      // point, or if the mid point is in a deletion, the base after the deletion.
      val midPoint  = { // Need to ensure the midpoint between 5' ends actually occurs in the overlapped region
        val retval = (rec.start + mate.end) / 2
        if (retval > rec.end) rec.end // trim only the minus strand read
        else if (retval < mate.start) mate.start - 1 // trim only positive-strand read
        else retval
      }
      val readEnd   = rec.readPosAtRefPos(pos=midPoint, returnLastBaseIfDeleted=true)
      val mateStart = { // NB: need to be careful if the midpoint falls in a deletion
        val retval = mate.readPosAtRefPos(pos=midPoint + 1, returnLastBaseIfDeleted=false)
        if (retval != 0) retval else mate.readPosAtRefPos(pos=midPoint + 1, returnLastBaseIfDeleted=true) + 1
      }
      val numOverlappingBasesRead = this.clip3PrimeEndOfRead(rec, rec.cigar.trailingHardClippedBases + rec.length - readEnd)
      val numOverlappingBasesMate = this.clip3PrimeEndOfRead(mate, mate.cigar.leadingHardClippedBases + mateStart - 1)
      (numOverlappingBasesRead, numOverlappingBasesMate)
    }
  }

  /** Returns the number of query bases at the 3' end of a record that extend past the 3' end of its mate, for FR pairs,
    * and zero otherwise.
    *
    * The two reads of an FR pair sequence the same insert from opposite ends, so beyond the last reference position
    * that both alignments cover they must have sequenced the same number of insert bases.  Any query bases a record has
    * beyond that count are read-through past the end of the insert (i.e. into the adapter).  Both counts are taken in
    * _query_ space so that an indel near the end of either alignment does not skew the comparison; the shared reference
    * position is used only as the landmark from which to count.
    *
    * Soft-clipped bases are counted, since they are query bases that may still belong to the insert.  Hard-clipped
    * bases are not, since they are not present in the record.
    *
    * When the two alignments share no reference position there is no landmark to count from, and the count falls back
    * to `numBasesExtendingPastDisjointMate` below.
    *
    * @param rec the record to examine
    * @param mateStart the alignment start of the mate
    * @param mateCigar the alignment cigar of the mate
    * @return the number of query bases to clip from the 3' end of the record (in sequencing order)
    */
  def numBasesExtendingPastMate(rec: SamRecord, mateStart: Int, mateCigar: Cigar): Int = {
    if (!rec.isFrPair) 0 else {
      val mateEnd      = mateStart + mateCigar.lengthOnTarget - 1
      val overlapStart = Math.max(rec.start, mateStart)
      val overlapEnd   = Math.min(rec.end, mateEnd)
      if (overlapEnd < overlapStart) numBasesExtendingPastDisjointMate(rec, mateStart, mateCigar, mateEnd)
      else if (rec.positiveStrand) {
        // The 3' end of a positive strand record is its right-hand end, so count from the last shared position.
        val recBases  = numQueryBasesAfter(start=rec.start, cigar=rec.cigar, refPos=overlapEnd)
        val mateBases = numQueryBasesAfter(start=mateStart, cigar=mateCigar, refPos=overlapEnd)
        Math.max(0, recBases - mateBases)
      }
      else {
        // The 3' end of a negative strand record is its left-hand end, so count from the first shared position.
        val recBases  = numQueryBasesBefore(start=rec.start, cigar=rec.cigar, refPos=overlapStart)
        val mateBases = numQueryBasesBefore(start=mateStart, cigar=mateCigar, refPos=overlapStart)
        Math.max(0, recBases - mateBases)
      }
    }
  }

  /** Returns the number of query bases at the 3' end of a record that extend past the 3' end of its mate, when the two
    * alignments cover no reference position in common.
    *
    * With no shared position there is no landmark to count query bases from, so the record's own soft-clipped tail is
    * projected into reference space at one query base per reference base and compared against the mate's
    * un-soft-clipped far end — the estimate this function has applied since #842.  Two reads can absolutely have read
    * through into each other's adapter here: both alignments lie within the insert, so a short insert read from both
    * ends leaves them disjoint whenever each alignment stops short of the other's, which soft clipping at a variant or
    * a mis-called base is enough to produce.  The extrapolation is exact whenever the region between the two
    * alignments is ungapped, and it is continuous with the query-space count above, which reduces to exactly this
    * expression as the shared span shrinks to a single reference position.
    *
    * The reference/query conflation that the query-space count exists to avoid cannot bite here: neither alignment
    * places a base between the two, so neither can place an indel between them either.  An indel in the sample that
    * falls in that gap is invisible to both reads and shifts the estimate by its own length, in either direction; the
    * result is capped by the record's soft clipping regardless, so no aligned base is ever removed on the strength of
    * it.
    *
    * A record whose 3' end points away from its mate — a positive strand record lying entirely to the right of its
    * mate, or a negative strand record entirely to its left — describes an outward-facing template rather than
    * read-through, and nothing is clipped.
    *
    * @param rec the record to examine
    * @param mateStart the alignment start of the mate
    * @param mateCigar the alignment cigar of the mate
    * @param mateEnd the alignment end of the mate, which the only caller has already computed
    * @return the number of query bases to clip from the 3' end of the record (in sequencing order)
    */
  private def numBasesExtendingPastDisjointMate(rec: SamRecord, mateStart: Int, mateCigar: Cigar, mateEnd: Int): Int = {
    if (rec.positiveStrand && rec.end < mateStart) {
      val mateUnSoftClippedEnd = mateEnd + mateCigar.trailingSoftClippedBases
      Math.max(0, rec.cigar.trailingSoftClippedBases - (mateUnSoftClippedEnd - rec.end))
    }
    else if (rec.negativeStrand && rec.start > mateEnd) {
      val mateUnSoftClippedStart = mateStart - mateCigar.leadingSoftClippedBases
      Math.max(0, rec.cigar.leadingSoftClippedBases - (rec.start - mateUnSoftClippedStart))
    }
    else 0
  }

  /** Returns the number of query bases in an alignment that lie to the right of the given reference position.
    *
    * Inserted and soft-clipped bases are counted when they follow the last query base at or before the given position.
    * If the position falls within a deletion or skipped region then no query base is at the position, and every query
    * base after the gap is counted.  Hard-clipped bases are never counted since they are not present in the record.
    *
    * @param start the alignment start of the alignment
    * @param cigar the cigar of the alignment
    * @param refPos the reference position, which must fall within the aligned span of the alignment
    */
  private def numQueryBasesAfter(start: Int, cigar: Cigar, refPos: Int): Int = {
    require(start <= refPos && refPos <= start + cigar.lengthOnTarget - 1,
      s"Reference position $refPos is outside the aligned span of $cigar at $start.")
    var numAtOrBefore = 0
    var refCur        = start
    var done          = false
    val iter          = cigar.iterator
    while (!done && iter.hasNext) {
      val elem = iter.next()
      if (elem.operator.consumesReferenceBases()) {
        val refLast = refCur + elem.length - 1
        if (refLast <= refPos) { numAtOrBefore += elem.lengthOnQuery; refCur = refLast + 1 }
        else { // the position falls within this element; count only the query bases up to and including it
          if (elem.operator.consumesReadBases()) numAtOrBefore += refPos - refCur + 1
          done = true
        }
      }
      else if (refCur > refPos) done = true // the element lies to the right of the position
      else numAtOrBefore += elem.lengthOnQuery
    }
    cigar.lengthOnQuery - numAtOrBefore
  }

  /** Returns the number of query bases in an alignment that lie to the left of the given reference position.
    *
    * Inserted and soft-clipped bases are counted when they precede the first query base at or after the given position.
    * If the position falls within a deletion or skipped region then no query base is at the position, and every query
    * base before the gap is counted.  Hard-clipped bases are never counted since they are not present in the record.
    *
    * @param start the alignment start of the alignment
    * @param cigar the cigar of the alignment
    * @param refPos the reference position, which must fall within the aligned span of the alignment
    */
  private def numQueryBasesBefore(start: Int, cigar: Cigar, refPos: Int): Int = {
    require(start <= refPos && refPos <= start + cigar.lengthOnTarget - 1,
      s"Reference position $refPos is outside the aligned span of $cigar at $start.")
    var numBefore = 0
    var refCur    = start
    var done      = false
    val iter      = cigar.iterator
    while (!done && iter.hasNext) {
      val elem = iter.next()
      if (elem.operator.consumesReferenceBases()) {
        val refLast = refCur + elem.length - 1
        if (refLast < refPos) { numBefore += elem.lengthOnQuery; refCur = refLast + 1 }
        else { // the position falls within this element; count only the query bases strictly before it
          if (elem.operator.consumesReadBases()) numBefore += refPos - refCur
          done = true
        }
      }
      else numBefore += elem.lengthOnQuery // refCur is always <= refPos here, so the element lies to the left
    }
    numBefore
  }

  /** Clips the reads in FR read pairs whose alignments extend beyond the far end of their mate's alignment.
    *
    * @param rec the read
    * @param mate the mate
    * @return the additional number of bases clipped (3' end in sequencing order) for the read and mate respectively
    */
  def clipExtendingPastMateEnds(rec: SamRecord, mate: SamRecord): (Int, Int) = {
    if (!rec.isFrPair) (0, 0)
    else {
      // NB: both amounts must be computed before either read is clipped, otherwise the second read is compared against
      // an alignment that has already been shortened (or, if the first read was clipped away entirely, un-mapped).
      val numBases1 = numBasesExtendingPastMate(rec=rec, mateStart=mate.start, mateCigar=mate.cigar)
      val numBases2 = numBasesExtendingPastMate(rec=mate, mateStart=rec.start, mateCigar=rec.cigar)
      (clipExtendingPastMateEnd(rec, numBases1), clipExtendingPastMateEnd(mate, numBases2))
    }
  }

  /** Clips the given number of query bases from the 3' end of a read.
    *
    * @param rec the record to clip
    * @param numBases the number of query bases present in the record to clip from its 3' end
    * @return the additional number of bases clipped (3' end in sequencing order)
    */
  private def clipExtendingPastMateEnd(rec: SamRecord, numBases: Int): Int = {
    if (numBases <= 0) 0 else {
      // NB: clip3PrimeEndOfRead takes the total clipping desired _including_ any existing clipping, while numBases
      // counts only bases present in the record, so any existing hard-clipping at that end must be added back in.
      val cigar                = rec.cigar
      val existingHardClipping = if (rec.negativeStrand) cigar.leadingHardClippedBases else cigar.trailingHardClippedBases
      this.clip3PrimeEndOfRead(rec, numBases + existingHardClipping)
    }
  }

  /** Computes the number of bases that are available to be clipped in a mapped SamRecord. */
  private def numberOfClippableBases(rec: SamRecord): Int = {
    rec.cigar.filter(e => e.operator.isAlignment || e.operator == Op.INSERTION).map(_.length).sum
  }

  /**
    * Private method that is the workhorse of all the clipping methods. _Always_ works at the start of
    * the Cigar/bases/quals, and _always_ attempts to clip >= numberOfBasesToClip.
    *
    * @param originalElems the elements in the cigar string to be clipped
    * @param numberOfBasesToClip the number of bases to attempt to clip
    * @param bases the array of bases for the read to be clipped
    * @param quals the array of quality scores for the read to be clipped
    * @return a tuple of (newCigarElems, readBasesClipped, refBasesClipped, newBases, newQuals)
    */
  private def clip(originalElems: Seq[CigarElem],
                   numberOfBasesToClip: Int,
                   bases: Array[Byte],
                   quals: Array[Byte]): (IndexedSeq[CigarElem], Int, Int, Array[Byte], Array[Byte]) = {

    if (originalElems.exists(_.operator == Op.PADDING)) throw new IllegalArgumentException("Can't handle reads with padding.")

    val existingHardClip = originalElems.takeWhile(_.operator == Op.HARD_CLIP).map(_.length).sum
    val existingSoftClip = originalElems.dropWhile(_.operator == Op.HARD_CLIP).takeWhile(_.operator == Op.SOFT_CLIP).map(_.length).sum
    val postClipElems = originalElems.dropWhile(e => e.operator == Op.HARD_CLIP || e.operator == Op.SOFT_CLIP).iterator.bufferBetter
    var readBasesClipped = 0
    var refBasesClipped  = 0
    val newElems = ArrayBuffer[CigarElem]()

    // The loop skips over all operators that are getting turned into clipping, while keeping track of
    // how many reference bases and how many read bases are skipped over.  If the clipping point falls
    // between existing operators then the `newElems` buffer is empty at the end of the while loop. If
    // the clip point falls within:
    //    a) a match/mismatch operator then the operator is split and the remainder added to the buffer
    //    b) an insertion: the remainder of the insertion is also clipped
    // If the operator immediately after the clip is a deletion, it is also discarded.
    //
    // At the end of the while loop newElems is either:
    //   a) Empty
    //   b) Contains a single element which is the remainder of an element that had to be split
    while (readBasesClipped < numberOfBasesToClip ||
      (readBasesClipped == numberOfBasesToClip && newElems.isEmpty && postClipElems.hasNext && postClipElems.head.operator == Op.DELETION)) {
      val elem = postClipElems.next()
      val op   = elem.operator
      val len  = elem.length

      if (op.consumesReadBases() && len > (numberOfBasesToClip - readBasesClipped)) {
        op match {
          case Op.INSERTION => readBasesClipped += len
          case _ =>
            val remainingClip = numberOfBasesToClip - readBasesClipped
            val remainingLength = len - remainingClip
            readBasesClipped += remainingClip
            refBasesClipped  += remainingClip
            newElems += CigarElem(op, remainingLength)
        }
      }
      else {
        if (op.consumesReadBases()) readBasesClipped     += len
        if (op.consumesReferenceBases()) refBasesClipped += len
      }
    }

    // Add the remainder of the post-clipping elements that haven't been turned into clips
    newElems ++= postClipElems

    // Prepend the appropriate clipping elements & fix up the read
    if (mode == ClippingMode.Hard) {
      val addedHardClip = existingSoftClip + readBasesClipped
      val totalHardClip = existingHardClip + addedHardClip
      newElems.prepend(CigarElem(Op.HARD_CLIP, totalHardClip))
      (newElems.toIndexedSeq, readBasesClipped, refBasesClipped, bases.drop(addedHardClip), quals.drop(addedHardClip))
    }
    else {
      val addedSoftClip = readBasesClipped
      val totalSoftClip = existingSoftClip + addedSoftClip
      newElems.prepend(CigarElem(Op.SOFT_CLIP, totalSoftClip))
      if (existingHardClip > 0) newElems.prepend(CigarElem(Op.HARD_CLIP, existingHardClip))
      if (mode == ClippingMode.SoftWithMask) hardMaskStartOfRead(bases, quals, totalSoftClip)
      (newElems.toIndexedSeq, readBasesClipped, refBasesClipped, bases, quals)
    }
  }

  /** Hard masks (to N) a number of bases starting from the start of the read. */
  private def hardMaskStartOfRead(bases: Array[Byte], quals: Array[Byte], maskBases: Int): Unit = {
    forloop (from=0, until=maskBases) { i =>
      bases(i) = NoCallBase
      quals(i) = NoCallQual
    }
  }

  /** Clips extended attributes that are the same length as the bases were prior to clipping. */
  protected def clipExtendedAttributes(rec: SamRecord, remove: Int, fromStart: Boolean): Unit = {
    if (mode == ClippingMode.Hard && remove > 0 && autoClipAttributes) {
      val newLength = rec.length
      val oldLength = newLength + remove
      rec.attributes.foreach {
        case (tag, s: String)   if s.length == oldLength => rec(tag) = if (fromStart) s.drop(remove) else s.take(newLength)
        case (tag, a: Array[_]) if a.length == oldLength => rec(tag) = if (fromStart) a.drop(remove) else a.take(newLength)
        case _  => ()
      }
    }
  }

  /**
    * Ensures sufficient masking or hard clipping exists on reads that may have soft-clipping. The read will
    * only be altered if:
    *   1. ClippingMode is Hard or SoftWithMask
    *   2. The read is already clipped at the appropriate end
    *   3. Soft-clipping exist within the first `length` bases of the appropriate end
    *
    * If all of those conditions are met, the soft-clipping will be upgraded to masking or hard clipping
    * up until length bases are clipped or masked. E.g. if the incoming cigar is `10H10S80M` and `length` is
    * 15, the resulting cigar will be `15H5S80M`.
    */
  protected def upgradeClipping(rec: SamRecord, length: Int, fromStart: Boolean): Unit = if (mode != ClippingMode.Soft && length > 0) {
    def iter = if (fromStart) rec.cigar.iterator else rec.cigar.reverseIterator
    val hardClipped = iter.takeWhile(_.operator == Op.H).map(_.length).sum
    val softClipped = iter.dropWhile(_.operator == Op.H).takeWhile(_.operator == Op.S).map(_.length).sum

    // If the requested length isn't all hard-clipped, and some of it's soft-clipped, then do the thing!
    if (hardClipped < length && softClipped > 0) {
      val lengthToUpgrade = math.min(softClipped, length - hardClipped)
      var (elems, bases, quals) = {
        if (fromStart) (rec.cigar.coalesce.elems, rec.bases, rec.quals)
        else (rec.cigar.coalesce.elems.reverse, rec.bases.reverse, rec.quals.reverse)
      }

      mode match {
        case ClippingMode.Soft =>
          unreachable("Should never reach here when mode == Soft")
        case ClippingMode.SoftWithMask =>
          hardMaskStartOfRead(bases, quals, lengthToUpgrade)
        case ClippingMode.Hard =>
          bases = bases.drop(lengthToUpgrade)
          quals = quals.drop(lengthToUpgrade)
          val newElems = IndexedSeq.newBuilder[CigarElem]

          elems match {
            case CigarElem(Op.H, h) +: CigarElem(Op.S, s) +: remaining =>
              newElems += CigarElem(Op.H, h+lengthToUpgrade)
              if (s > lengthToUpgrade) newElems += CigarElem(Op.S, s-lengthToUpgrade)
              newElems ++= remaining
            case CigarElem(Op.S, s) +: remaining =>
              newElems += CigarElem(Op.H, lengthToUpgrade)
              if (s > lengthToUpgrade) newElems += CigarElem(Op.S, s-lengthToUpgrade)
              newElems ++= remaining
            case es =>
              unreachable(s"Cigar $es doesn't contain soft clipping at the start")
          }

          elems = newElems.result()
      }

      if (fromStart) {
        rec.cigar = Cigar(elems)
        rec.bases = bases
        rec.quals = quals
      }
      else {
        rec.cigar = Cigar(elems.reverse)
        rec.bases = bases.reverse
        rec.quals = quals.reverse
      }

      clipExtendedAttributes(rec, lengthToUpgrade, fromStart=fromStart)
    }
  }

  /** Invalidates the set of tags that cannot be trusted if clipping is applied to a read. */
  protected def cleanupClippedRecord(rec: SamRecord): Unit = TagsToInvalidate.foreach(tag => rec(tag) = null)
}
