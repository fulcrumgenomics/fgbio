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
import htsjdk.samtools.CigarOperator

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
  override def toString: String = length + operator.toString
}

/** Companion object for Cigar that offers alternative constructors. */
object Cigar {
  private val ZeroCharAsInt = '0'.toInt

  /** Constructs a Cigar object from a Cigar string. */
  def apply(cigar: String): Cigar = {
    require(cigar != null && cigar.length > 0, "Can't parse a null or empty cigar string.")
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
}

/**
  * Object representation of a Cigar string representing an alignment between two sequences.
  * @param elems the ordered sequence of elements in the Cigar
  */
case class Cigar(elems: IndexedSeq[CigarElem]) extends Iterable[CigarElem] {
  require(elems.nonEmpty, s"Cigar must have at least one element.")

  /** Provides an iterator over the elements in the cigar. */
  override def iterator: Iterator[CigarElem] = elems.iterator

  /** Returns the length of the alignment on the query sequence. */
  def lengthOnQuery: Int = elems.filter(_.operator.consumesReadBases()).map(_.length).sum

  /** Returns the length of the alignment on the query sequence. */
  def lengthOnTarget: Int = elems.filter(_.operator.consumesReferenceBases()).map(_.length).sum

  /** Returns the canonical Cigar string. */
  override def toString(): String = elems.mkString
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
    cigar.foreach {
      case CigarElem(CigarOperator.M | CigarOperator.EQ | CigarOperator.X, len) =>
        forloop (from=0, until=len) { i =>
          val q = this.query(qOffset).toChar
          val t = this.target(tOffset).toChar
          queryBuffer.append(q)
          alignBuffer.append(if (q == t) matchChar else mismatchChar)
          targetBuffer.append(t)
          qOffset += 1
          tOffset += 1
        }
      case CigarElem(CigarOperator.INSERTION, len) =>
        forloop (from=0, until=len) { i =>
          queryBuffer.append(this.query(qOffset).toChar)
          alignBuffer.append(gapChar)
          targetBuffer.append(padChar)
          qOffset += 1
        }
      case CigarElem(CigarOperator.DELETION, len) =>
        forloop (from=0, until=len) { i =>
          queryBuffer.append(padChar)
          alignBuffer.append(gapChar)
          targetBuffer.append(this.target(tOffset).toChar)
          tOffset += 1
        }
      case other =>
        throw new IllegalStateException(s"Padding string cannot support cigar operator $other")
    }

    buffers.map(_.toString())
  }
}
