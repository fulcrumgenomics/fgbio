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

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.SamRecordClipper.ClippingMode
import htsjdk.samtools.{Cigar, CigarElement, SAMRecord, SAMUtils, CigarOperator => Op}

import scala.collection.mutable.ListBuffer

/**
  * Holds types and constants related to SAMRecord clipping.
  */
object SamRecordClipper {
  /** An enumeration representing the various ways to clip bases within a read. */
  object ClippingMode extends Enumeration {
    val Soft, Hard, SoftWithMask = Value
  }

  type ClippingMode = ClippingMode.Value

  /** The set of tags that should be invalidated if a read undergoes clipping. */
  val TagsToInvalidate = Seq("MD", "NM", "UQ")

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
    */
  def clipStartOfAlignment(rec: SAMRecord, numberOfBasesToClip: Int) : Unit = {
    if (rec.getReadUnmappedFlag || numberOfBasesToClip < 1) {
      Unit
    }
    else if (numberOfClippableBases(rec) <= numberOfBasesToClip) {
      SAMUtils.makeReadUnmapped(rec)
    }
    else {
      val oldElems = rec.getCigar.toSeq
      val oldBases = rec.getReadBases
      val oldQuals = rec.getBaseQualities
      val (newElems, readBasesClipped, refBasesClipped, bases, quals) = clip(oldElems, numberOfBasesToClip, oldBases, oldQuals)
      rec.setAlignmentStart(rec.getAlignmentStart + refBasesClipped)
      rec.setCigar(new Cigar(util.Arrays.asList(newElems:_*)))
      rec.setReadBases(bases)
      rec.setBaseQualities(quals)
      cleanupClippedRecord(rec)

      if (mode == ClippingMode.Hard && readBasesClipped > 0 && autoClipAttributes) {
        val oldLength = oldBases.length
        val shortenBy = oldLength - bases.length
        rec.getAttributes.foreach { kv =>
          kv.value match {
            case s: String   if s.length == oldLength => rec.setAttribute(kv.tag, s.drop(shortenBy))
            case a: Array[_] if a.length == oldLength => rec.setAttribute(kv.tag, a.drop(shortenBy))
            case _  => Unit
          }
        }
      }
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
    */
  def clipEndOfAlignment(rec: SAMRecord, numberOfBasesToClip: Int) : Unit = {
    if (rec.getReadUnmappedFlag || numberOfBasesToClip < 1) {
      Unit
    }
    else if (numberOfClippableBases(rec) <= numberOfBasesToClip) {
      SAMUtils.makeReadUnmapped(rec)
    }
    else {
      val oldElems = rec.getCigar.toSeq.reverse
      val oldBases = rec.getReadBases.reverse
      val oldQuals = rec.getBaseQualities.reverse
      val (newElems, readBasesClipped, refBasesClipped, bases, quals) = clip(oldElems, numberOfBasesToClip, oldBases, oldQuals)
      rec.setCigar(new Cigar(util.Arrays.asList(newElems.reverse:_*)))
      rec.setReadBases(bases.reverse)
      rec.setBaseQualities(quals.reverse)
      cleanupClippedRecord(rec)

      if (mode == ClippingMode.Hard && readBasesClipped > 0 && autoClipAttributes) {
        val oldLength = oldBases.length
        val newLength = bases.length
        rec.getAttributes.foreach { kv =>
          kv.value match {
            case s: String   if s.length == oldLength => rec.setAttribute(kv.tag, s.take(newLength))
            case a: Array[_] if a.length == oldLength => rec.setAttribute(kv.tag, a.take(newLength))
            case _  => Unit
          }
        }
      }
    }
  }

  /**
    * Attempts to clip an additional numberOfBasesToClip from the 5' end of the read. For
    * details see [[clipStartOfAlignment()]] and [[clipEndOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    */
  def clip5PrimeEndOfAlignment(rec: SAMRecord, numberOfBasesToClip: Int): Unit = {
    if (rec.getReadNegativeStrandFlag) clipEndOfAlignment(rec, numberOfBasesToClip)
    else clipStartOfAlignment(rec, numberOfBasesToClip)
  }

  /**
    * Attempts to clip an additional numberOfBasesToClip from the 5' end of the read. For
    * details see [[clipStartOfAlignment()]] and [[clipEndOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    */
  def clip3PrimeEndOfAlignment(rec: SAMRecord, numberOfBasesToClip: Int): Unit = {
    if (rec.getReadNegativeStrandFlag) clipStartOfAlignment(rec, numberOfBasesToClip)
    else clipEndOfAlignment(rec, numberOfBasesToClip)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the start (left-hand end)
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[clipStartOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param clipLength the total amount of clipping desired, including any existing clipping
    */
  def clipStartOfRead(rec: SAMRecord, clipLength: Int): Unit = {
    val existingClipping = rec.getCigar.getCigarElements.takeWhile(_.getOperator.isClipping).map(_.getLength).sum
    if (clipLength > existingClipping) clipStartOfAlignment(rec, clipLength - existingClipping)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the end (right-hand end)
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[clipStartOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param clipLength the total amount of clipping desired, including any existing clipping
    */
  def clipEndOfRead(rec: SAMRecord, clipLength: Int): Unit = {
    val existingClipping = rec.getCigar.getCigarElements.toSeq.reverse.takeWhile(_.getOperator.isClipping).map(_.getLength).sum
    if (clipLength > existingClipping) clipEndOfAlignment(rec, clipLength - existingClipping)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the 5' end
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[clipStartOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    */
  def clip5PrimeEndOfRead(rec: SAMRecord, numberOfBasesToClip: Int): Unit = {
    if (!rec.getReadUnmappedFlag && rec.getReadNegativeStrandFlag) clipEndOfRead(rec, numberOfBasesToClip)
    else clipStartOfRead(rec, numberOfBasesToClip)
  }

  /**
    * Ensures that there are at least clipLength bases clipped at the 3' end
    * of the read, _including_ any existing soft and hard clipping.  Calculates any
    * additional clipping and delegates to [[clipStartOfAlignment()]].
    *
    * @param rec the record to be clipped
    * @param numberOfBasesToClip the number of additional bases to be clipped
    */
  def clip3PrimeEndOfRead(rec: SAMRecord, numberOfBasesToClip: Int): Unit = {
    if (!rec.getReadUnmappedFlag && rec.getReadNegativeStrandFlag) clipStartOfRead(rec, numberOfBasesToClip)
    else clipEndOfRead(rec, numberOfBasesToClip)
  }

  /** Computes the number of bases that are available to be clipped in a mapped SAMRecord. */
  private def numberOfClippableBases(rec: SAMRecord): Int = {
    rec.getCigar.filter(e => e.getOperator.isAlignment || e.getOperator == Op.INSERTION).map(_.getLength).sum
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
  private def clip(originalElems: Seq[CigarElement],
                   numberOfBasesToClip: Int,
                   bases: Array[Byte],
                   quals: Array[Byte]): (Seq[CigarElement], Int, Int, Array[Byte], Array[Byte]) = {

    if (originalElems.exists(_.getOperator == Op.PADDING)) throw new IllegalArgumentException("Can't handle reads with padding.")

    val existingHardClip = originalElems.takeWhile(_.getOperator == Op.HARD_CLIP).map(_.getLength).sum
    val existingSoftClip = originalElems.dropWhile(_.getOperator == Op.HARD_CLIP).takeWhile(_.getOperator == Op.SOFT_CLIP).map(_.getLength).sum
    val postClipElems = originalElems.dropWhile(e => e.getOperator == Op.HARD_CLIP || e.getOperator == Op.SOFT_CLIP).iterator.bufferBetter
    var readBasesClipped = 0
    var refBasesClipped  = 0
    val newElems = ListBuffer[CigarElement]()

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
      (readBasesClipped == numberOfBasesToClip && newElems.isEmpty && postClipElems.hasNext && postClipElems.head.getOperator == Op.DELETION)) {
      val elem = postClipElems.next()
      val op   = elem.getOperator
      val len  = elem.getLength

      if (op.consumesReadBases() && len > (numberOfBasesToClip - readBasesClipped)) {
        op match {
          case Op.INSERTION => readBasesClipped += len
          case _ =>
            val remainingClip = numberOfBasesToClip - readBasesClipped
            val remainingLength = len - remainingClip
            readBasesClipped += remainingClip
            refBasesClipped  += remainingClip
            newElems += new CigarElement(remainingLength, op)
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
      newElems.prepend(new CigarElement(totalHardClip, Op.HARD_CLIP))
      (newElems, readBasesClipped, refBasesClipped, bases.drop(addedHardClip), quals.drop(addedHardClip))
    }
    else {
      val addedSoftClip = readBasesClipped
      val totalSoftClip = existingSoftClip + readBasesClipped
      newElems.prepend(new CigarElement(totalSoftClip, Op.SOFT_CLIP))
      if (existingHardClip > 0) newElems.prepend(new CigarElement(existingHardClip, Op.HARD_CLIP))
      if (mode == ClippingMode.SoftWithMask) hardMaskStartOfRead(bases, quals, totalSoftClip)
      (newElems, readBasesClipped, refBasesClipped, bases, quals)
    }
  }

  /** Hard masks (to N) a number of bases starting from the start of the read. */
  private def hardMaskStartOfRead(bases: Array[Byte], quals: Array[Byte], maskBases: Int): Unit = {
    forloop (from=0, until=maskBases) { i =>
      bases(i) = NoCallBase
      quals(i) = NoCallQual
    }
  }

  /** Invalidates the set of tags that cannot be trusted if clipping is applied to a read. */
  protected def cleanupClippedRecord(rec: SAMRecord): Unit = TagsToInvalidate.foreach(tag => rec.setAttribute(tag, null))
}
