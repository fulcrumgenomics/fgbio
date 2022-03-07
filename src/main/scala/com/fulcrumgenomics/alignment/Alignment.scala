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

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools.{CigarElement, CigarOperator, SAMRecord, Cigar => HtsJdkCigar}

import scala.collection.mutable.ArrayBuffer

/**
  * Represents an element in a Cigar.
  *
  * @param operator the type of element (e.g. match, insertion, etc.)
  * @param length the length of the element in bases (must be greater than 0).
  */
case class CigarElem(operator: CigarOperator, length: Int) {
  require(length > 0, s"Cigar element must have length >= 0. Operator: $operator, Length: $length")

  /** Returns how this element should be represented in a cigar string. */
  override def toString: String = s"${length}${operator}"

  /** Returns the length of the element on the query sequence. */
  def lengthOnQuery: Int = if (operator.consumesReadBases()) length else 0

  /** Returns the length of the element on the query sequence. */
  def lengthOnTarget: Int = if (operator.consumesReferenceBases()) length else 0
}

/** Companion object for Cigar that offers alternative constructors. */
object Cigar {
  private val ZeroCharAsInt = '0'.toInt

  /** An empty Cigar object. */
  val empty = Cigar(IndexedSeq.empty)

  /** Constructs a Cigar object from a Cigar string. */
  def apply(cigar: String): Cigar = {
    require(cigar != null && cigar.nonEmpty, "Can't parse a null or empty cigar string.")
    val elements = new ArrayBuffer[CigarElem]
    val chars    = cigar.toCharArray

    var i = 0
    while (i < chars.length) {
      // Parse out the length of the operator
      val expectedNumberAtPosition = i
      var length = 0
      while (Character.isDigit(chars(i))) { length = (length * 10) + chars(i).toInt - ZeroCharAsInt; i += 1 }
      require(i > expectedNumberAtPosition, s"Malformed CIGAR string $cigar, found non-digit ${chars(i)} where number expected.")

      // Grab the operator (characterToEnum will throw an exception if the operator is invalid)
      val operator = CigarOperator.characterToEnum(chars(i))
      elements += CigarElem(operator, length)
      i += 1
    }

    Cigar(elements.toIndexedSeq)
  }

  /**
    * Constructs a Cigar from a Cigar string, allowing the string to be either `null` or
    * `*`, in which case Cigar.empty is returned.  Useful for parsing Cigars from SAM files.
    */
  def fromSam(cigar: String): Cigar = {
    if (cigar == null || cigar == SAMRecord.NO_ALIGNMENT_CIGAR) Cigar.empty else apply(cigar)
  }

  /** Constructs a Cigar objects from the HTSJDK class of the same name. */
  def apply(ciggy: HtsJdkCigar): Cigar = Cigar(ciggy.iterator().map(e => CigarElem(e.getOperator, e.getLength)).toIndexedSeq)
}

/** A data class for holding statistics about a [[Cigar]]. */
private[alignment] case class CigarStats(lengthOnQuery: Int,
                                         lengthOnTarget: Int,
                                         alignedBases: Int,
                                         clippedBases: Int,
                                         leadingSoftClippedBases: Int,
                                         trailingSoftClippedBases: Int)

/** Companion object for [[Cigar]]. */
private[alignment] object CigarStats {

  /** Build a [[CigarStats]] from a [[Cigar]]. */
  private[alignment] def apply(cigar: Cigar): CigarStats = {
    var lengthOnQuery            = 0
    var lengthOnTarget           = 0
    var alignedBases             = 0
    var clippedBases             = 0
    var leadingSoftClippedBases  = 0
    var trailingSoftClippedBases = 0
    var index                    = 0

    // leading clipping
    while (index < cigar.elems.length && cigar.elems(index).operator.isClipping) {
      val elem = cigar.elems(index)
      lengthOnQuery  += elem.lengthOnQuery
      lengthOnTarget += elem.lengthOnTarget
      clippedBases += elem.length
      if (elem.operator == CigarOperator.S) leadingSoftClippedBases += elem.length
      index += 1
    }

    // aligned bases
    while (index < cigar.elems.length && !cigar.elems(index).operator.isClipping) {
      val elem = cigar.elems(index)
      lengthOnQuery  += elem.lengthOnQuery
      lengthOnTarget += elem.lengthOnTarget
      if (elem.operator.isAlignment) alignedBases += elem.length
      index += 1
    }

    // trailing clipping
    while (index < cigar.elems.length && cigar.elems(index).operator.isClipping) {
      val elem = cigar.elems(index)
      lengthOnQuery  += elem.lengthOnQuery
      lengthOnTarget += elem.lengthOnTarget
      clippedBases += elem.length
      if (elem.operator == CigarOperator.S) trailingSoftClippedBases += elem.length
      index += 1
    }

    require(index == cigar.elems.length,
      s"Invalid cigar, did not consume all the elements: ${Cigar(cigar.elems.take(index))}[${Cigar(cigar.elems.drop(index))}]"
    )

    CigarStats(
      lengthOnQuery            = lengthOnQuery,
      lengthOnTarget           = lengthOnTarget,
      alignedBases             = alignedBases,
      clippedBases             = clippedBases,
      leadingSoftClippedBases  = leadingSoftClippedBases,
      trailingSoftClippedBases = trailingSoftClippedBases
    )
  }

  /** Build a [[CigarStats]] from a string representation of the cigar. */
  private[alignment] def apply(cigar: String): CigarStats = this.apply(Cigar(cigar))
}

/**
  * Object representation of a Cigar string representing an alignment between two sequences.
  * @param elems the ordered sequence of elements in the Cigar
  */
case class Cigar(elems: IndexedSeq[CigarElem]) extends Iterable[CigarElem] {
  // Cache whether or not the Cigar is coalesced already (i.e. has no pair of adjacent elements with the same operator)
  private lazy val isCoalesced: Boolean = {
    var itIs = true
    var index = 0
    while (index < elems.length-1 && itIs) {
      itIs = elems(index).operator != elems(index+1).operator
      index += 1
    }
    itIs
  }

  /** Statistics about this cigar which are computed once on first access and then cached. */
  private lazy val stats: CigarStats = CigarStats(this)

  /** Provides an iterator over the elements in the cigar. */
  override def iterator: Iterator[CigarElem] = elems.iterator

  /** Provides an iterator over the elements in the reverse order. */
  def reverseIterator: Iterator[CigarElem] = elems.reverseIterator

  /** Returns the length of the alignment on the query sequence. */
  def lengthOnQuery: Int = stats.lengthOnQuery

  /** Returns the length of the alignment on the target sequence. */
  def lengthOnTarget: Int = stats.lengthOnTarget

  /** Yields a new cigar that is truncated to the given length on the query. */
  def truncateToQueryLength(len: Int): Cigar = truncate(len, e => e.operator.consumesReadBases())

  /** Yields a new cigar that is truncated to the given length on the target. */
  def truncateToTargetLength(len: Int): Cigar = truncate(len, e => e.operator.consumesReferenceBases())

  /** Returns the number of bases that are directly aligned between the two sequences. */
  def alignedBases: Int = stats.alignedBases

  /** Returns the number of bases that are clipped between the two sequences. */
  def clippedBases: Int = stats.clippedBases

  /** Returns the number of bases that are soft-clipped at the start of the sequence  Ignores hard-clips. */
  def leadingSoftClippedBases: Int = stats.leadingSoftClippedBases

  /** Returns the number of bases that are soft-clipped at the end of this sequence.  Ignores hard-clips. */
  def trailingSoftClippedBases: Int = stats.trailingSoftClippedBases

  /** Returns the number of bases that are hard-clipped at the start of the sequence. */
  def leadingHardClippedBases = this.headOption.map { elem =>
    if (elem.operator == CigarOperator.H) elem.length else 0
  }.getOrElse(0)

  /** Returns the number of bases that are clipped at the start of the sequence. */
  def leadingClippedBases: Int = leadingHardClippedBases + leadingSoftClippedBases

  /** Returns the number of bases that are hard-clipped at the end of the sequence. */
  def trailingHardClippedBases = this.lastOption.map { elem =>
    if (elem.operator == CigarOperator.H) elem.length else 0
  }.getOrElse(0)

  /** Returns the number of bases that are clipped at the end of the sequence. */
  def trailingClippedBases: Int = trailingHardClippedBases + trailingSoftClippedBases

  /** Truncates the cigar based on either query or target length cutoff. */
  private def truncate(len: Int, shouldCount: CigarElem => Boolean): Cigar = {
    var pos = 1
    val iter = iterator
    val builder = IndexedSeq.newBuilder[CigarElem]
    while (pos <= len && iter.hasNext) {
      val elem = iter.next()
      if (shouldCount(elem)) {
        val maxElemLength = len - pos + 1
        builder += (if (elem.length <= maxElemLength) elem else elem.copy(length=maxElemLength))
        pos += elem.length
      }
      else {
        builder += elem
      }
    }

    Cigar(builder.result())
  }

  /**
    * Returns true if this Cigar is a prefix of the other cigar. Prefix is defined as having the same
    * operators for the same bases through the length of the cigar.  E.g.:
    *   10M is a prefix of 10M2S
    *   10M is a prefix of 10M
    *   10M is a prefix of 20M
    *   10M is not a prefix of 2S10M
    */
  def isPrefixOf(that: Cigar): Boolean = {
    if (that.elems.size < this.elems.size) {
      false
    }
    else {
      val lastIndex = this.elems.size - 1
      var isPrefix  = true
      var i = 0
      while (isPrefix && i <= lastIndex) {
        val lhs = this.elems(i)
        val rhs = that.elems(i)
        isPrefix = lhs.operator == rhs.operator && (if (i == lastIndex) lhs.length <= rhs.length else lhs.length == rhs.length)
        i += 1
      }

      isPrefix
    }
  }

  /** Returns a new Cigar that contains the same elements in the reverse order of this cigar. */
  def reverse: Cigar = Cigar(this.elems.reverse)

  /** Coalesces adjacent operators of the same type into single operators. */
  def coalesce: Cigar = {
    if (isCoalesced) {
      this
    }
    else {
      val builder = IndexedSeq.newBuilder[CigarElem]
      val iter   = iterator.bufferBetter

      while (iter.hasNext) {
        val elem = iter.next()
        val same = iter.takeWhile(_.operator == elem.operator).foldLeft(elem)((a,b) => CigarElem(a.operator, a.length + b.length))
        builder += same
      }

      Cigar(builder.result())
    }
  }

  /** Returns the canonical Cigar string. */
  override def toString(): String = elems.mkString

  /** Converts the Cigar into an HTSJDK Cigar object. */
  def toHtsjdkCigar: htsjdk.samtools.Cigar = {
    new htsjdk.samtools.Cigar(elems.iterator.map(e => new CigarElement(e.length, e.operator)).toJavaList)
  }
}


/** Companion object for Alignment. */
object Alignment {
  /** Construct an alignment using Strings for sequences instead of byte arrays. */
  def apply(query: String, target: String, queryStart: Int, targetStart: Int, cigar: Cigar, score: Int): Alignment = {
    new Alignment(query.getBytes, target.getBytes, queryStart, targetStart, cigar, score)
  }
}

/**
  * A general class to describe the alignment between two sequences or partial ranges thereof
  * @param query the query sequence
  * @param target the target sequence
  * @param queryStart the 1-based position in the query sequence where the alignment begins
  * @param targetStart the 1-based position in the target sequence where the alignment begins
  * @param cigar a [[Cigar]] object describing the alignment of the two sequences
  * @param score the alignment score
  */
case class Alignment(query: Array[Byte],
                     target: Array[Byte],
                     queryStart: Int,
                     targetStart: Int,
                     cigar: Cigar,
                     score: Int) {

  /** One based closed coordinate of the end of the alignment on the query sequence. */
  def queryEnd: Int = queryStart + cigar.lengthOnQuery - 1

  /** One based closed coordinate of the end of the alignment on the query sequence. */
  def targetEnd: Int = targetStart + cigar.lengthOnTarget - 1

  /**
    * Generates a padded text representation of the alignment for visualization. The returned
    * sequence will consist of three lines as follows (minus the labels on the left):
    *
    * query : ACGTGAACTGACT-ACTGTATGCG
    * align : |||||  |||||| ||||||||.|
    * target: ACGTG--CTGACTGACTGTATGGG
    *
    * @param matchChar the character to use in the alignment line when the bases match
    * @param mismatchChar the character to use in the alignment line when the bases mismatch
    * @param gapChar the character to use in the alignment line when one of the sequences is gapped
    * @param padChar the character to use in the query or target sequence lines when padding is required
    * @return Three lines representing the alignment
    */
  def paddedString(matchChar: Char = '|', mismatchChar: Char = '.', gapChar: Char = ' ', padChar: Char = '-'): Seq[String] = {
    val buffers = Seq(new StringBuilder, new StringBuilder, new StringBuilder)
    val Seq(queryBuffer, alignBuffer, targetBuffer) = buffers

    var qOffset = queryStart - 1
    var tOffset = targetStart - 1
    cigar.foreach { cig =>
      val len = cig.length

      cig.operator match {
        case CigarOperator.INSERTION =>
          forloop (from=0, until=len) { _ =>
            queryBuffer.append(this.query(qOffset).toChar)
            alignBuffer.append(gapChar)
            targetBuffer.append(padChar)
            qOffset += 1
          }
        case CigarOperator.DELETION =>
          forloop (from=0, until=len) { _ =>
            queryBuffer.append(padChar)
            alignBuffer.append(gapChar)
            targetBuffer.append(this.target(tOffset).toChar)
            tOffset += 1
          }
        case op =>
          forloop (from=0, until=len) { _ =>
            val q = this.query(qOffset).toChar
            val t = this.target(tOffset).toChar
            queryBuffer.append(q)
            targetBuffer.append(t)
            op match {
              case CigarOperator.EQ => alignBuffer.append(matchChar)
              case CigarOperator.X  => alignBuffer.append(mismatchChar)
              case CigarOperator.M  => alignBuffer.append(if (q == t) matchChar else mismatchChar)
              case other            => throw new IllegalStateException(s"Unsupported cigar operator: $other")
            }

            qOffset += 1
            tOffset += 1
          }
        }
      }

    buffers.map(_.toString())
  }

  /**
    * Returns a subset of the [[Alignment]] representing the region defined by `start` and `end` on the query
    * sequence. The returned alignment will contain the entire query and target sequences, but will have
    * adjusted `queryStart` and `targetStart` positions and an updated `cigar`.  The score is set to 0.
    *
    * @param start the 1-based inclusive position of the first base on the query sequence to include
    * @param end   the 1-based inclusive position of the last base on the query sequence to include
    * @return a new [[Alignment]] with updated coordinates and cigar
    */
  def subByQuery(start: Int, end: Int): Alignment = {
    require(start >= queryStart && start <= queryEnd, "start is outside of aligned region of target sequence")
    require(end   >= queryStart && end   <= queryEnd, "end is outside of aligned region of target sequence")
    sub(start, end, this.queryStart, _.operator.consumesReadBases())
  }

  /**
    * Returns a subset of the Alignment representing the region defined by `start` and `end` on the target
    * sequence. The returned alignment will contain the entire query and target sequences, but will have
    * adjusted `queryStart` and `targetStart` positions and an updated `cigar`.  The score is set to 0.
    *
    * @param start the 1-based inclusive position of the first base on the target sequence to include
    * @param end   the 1-based inclusive position of the last base on the target sequence to include
    * @return a new [[Alignment]] with updated coordinates and cigar
    */
  def subByTarget(start: Int, end: Int): Alignment = {
    require(start >= targetStart && start <= targetEnd, "start is outside of aligned region of target sequence")
    require(end   >= targetStart && end   <= targetEnd, "end is outside of aligned region of target sequence")
    sub(start, end, this.targetStart, _.operator.consumesReferenceBases())
  }

  /**
    * Private helper method that helps generate an Alignment that is a subset of the current alignment.
    *
    * @param start the start (on either the query or target sequence) of the desired window
    * @param end the end (on either the query or target sequence) of the desired window
    * @param initialStart the start position on the relevant sequence (either target OR query)
    * @param consumes a function that returns true if the cigar element passed as a parameter consumes
    *                 bases on the sequence on which `start` and `end` are expresssed
    * @return a new [[Alignment]] with adjusted start, end and cigar, and with score=0
    */
  private def sub(start: Int, end: Int, initialStart: Int, consumes: CigarElem => Boolean): Alignment = {
    val elems = IndexedSeq.newBuilder[CigarElem]
    var (qStart, tStart, currStart) = (queryStart, targetStart, initialStart) // currStart = start of current element

    this.cigar.foreach { elem =>
      // If the current element consumes the chosen sequence, calculate the end, else set to the same as start

      val elementConsumes = consumes(elem)
      val currEnd = if (elementConsumes) currStart + elem.length - 1 else currStart-1

      if (currEnd < start) {
        // Element before the desired window, need to bump start positions
        qStart += elem.lengthOnQuery
        tStart += elem.lengthOnTarget
        if (elementConsumes) currStart += elem.length
      }
      else if (currStart > end) {
        // Don't include, beyond the range we're interested
      }
      else if (currStart >= start && currEnd <= end) {
        // Contained within the target region
        // Only add the element if it a) consumes bases or b) isn't right at the boundary (to avoid trailing gaps)
        if (elementConsumes || currStart != start) {
          elems += elem
          if (elementConsumes) currStart += elem.length
        }
      }
      else {
        // Element overlaps the desired region, may enclose, be enclosed by or straddle
        var len = elem.length // <- the length of the element to be added, after any clipping

        if (currStart < start) {
          // Element is split over the start of the desired region
          val diff = start - currStart
          len -= diff
          if (elem.operator.consumesReadBases())      qStart += diff
          if (elem.operator.consumesReferenceBases()) tStart += diff
          currStart += diff
        }

        if (currEnd > end) {
          // Element is split over the end of the desired region
          len -= (currEnd - end)
        }

        elems += CigarElem(elem.operator, len)
        currStart += len
      }
    }

    copy(queryStart=qStart, targetStart=tStart, cigar=Cigar(elems.result()), score=0)
  }
}
