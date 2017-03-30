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
}
