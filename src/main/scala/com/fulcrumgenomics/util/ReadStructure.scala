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

package com.fulcrumgenomics.util

import scala.collection.immutable
import scala.reflect.ClassTag

// Future improvement: I think we really need to support having a final segment of indeterminate length. Probably it
// would be good to have support for * and + to indicate zero or more and one or more bases.
object ReadStructure {
  /** Creates a read structure from a string */
  def apply(readStructure: String): ReadStructure = new ReadStructure(segs=segments(readStructure=readStructure))

  /** Creates a sequence of read segments from a string. */
  private def segments(readStructure: String): Seq[ReadSegment] = {
    var segmentOffset = 0 // for each segment we keep track of the offset (# of bases) from the start of the read
    var startIndexOfDigits = 0 // the index into the `readStructure` string where the first digit character starts for the current segment (inclusive)
    // Parse each character and build up segments.
    val segments = readStructure.zipWithIndex.flatMap { case (c, idx) =>
      if (!c.isDigit) { // no more digits, so we are now at the start of the segment type
      val endIndexOfDigits  = idx // the index into the `readStructure` string where the first digit character ends for the current segment (exclusive)
        if (endIndexOfDigits <= startIndexOfDigits) { // if we had no digit characters, then throw an exception
          throw new IllegalArgumentException("Read structure missing a segment length: " + formatReadStructureError(readStructure, idx, idx+1))
        }
        val segmentLength = Integer.parseInt(readStructure.substring(startIndexOfDigits, endIndexOfDigits)) // parse the length of the segment
        if (segmentLength == 0) { // make sure we have a non-zero segment
          throw new IllegalArgumentException("Read structure had a segment with zero length: " + formatReadStructureError(readStructure, startIndexOfDigits, endIndexOfDigits))
        }
        val currentSegmentOffset = segmentOffset // store the offset for hte current segment we are building
        startIndexOfDigits = idx + 1 // update the index of the first digit character for the next segment
        segmentOffset += segmentLength // update the total # of bases consumed
        try { // this relies on the fact that the segment type is only a single character (`c`).
          Some(ReadSegment(offset = currentSegmentOffset, length = segmentLength, c = c)) // try building a read segment.
        } catch {
          case ex: IllegalArgumentException =>
            throw new IllegalArgumentException("Read structure had a segment of unknown type: " + formatReadStructureError(readStructure, idx, idx+1))
        }
      }
      else { // not ready yet
        None
      }
    }
    // make sure we do not have any trailing characters
    if (startIndexOfDigits < readStructure.length) throw new IllegalArgumentException("Read structure had trailing digits: " + formatReadStructureError(readStructure, startIndexOfDigits, readStructure.length))
    segments
  }

  /** Contains the bases and optionally base qualities that correspond to the given read segment. */
  case class SubRead(bases: String, quals: Option[String] = None, segment: ReadSegment)

  /**
    * Inserts square brackets around the characers in the read structure that are causing the error.
    *
    * @param s the read structure string
    * @param start the start of the error in the string (inclusive)
    * @param end the ned of the error in the string (exclusive)
    */
  private def formatReadStructureError(s: String, start: Int, end: Int) = {
    val prefix = s.substring(0, start)
    val error  = s.substring(start, end)
    val suffix = if (end == s.length) "" else s.substring(end, s.length)
    s"$prefix[$error]$suffix"
  }
}

/**
  * Describes the structure of a give read.  A read contains one or more read segments. A read segment describes
  * a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
  * of the read.
  *
  * @param segs the segments composing this read structure.
  * @param resetOffsets true if we are to update the offsets of each segment to be relative to each other, false
  *                     otherwise.  This may be useful to perform when segments are extracted from another read
  *                     structure.  New segment objects will be created in this case.
  */
class ReadStructure(segs: Seq[ReadSegment], resetOffsets: Boolean = false) extends immutable.Seq[ReadSegment] {
  import ReadStructure.SubRead

  /** Creates a read structure from a string representation (ex. 8M150T8B). */
  def this(readStructure: String) = this(ReadStructure.segments(readStructure=readStructure))

  private val segments: Seq[ReadSegment] = if (resetOffsets) {
    var idx = 0
    var off = 0
    val segments = new Array[ReadSegment](segs.length)
    while (idx < segments.length) {
      val segment = segs(idx)
      segments(idx) = ReadSegment(segment, off, segment.length)
      off += segment.length
      idx += 1
    }
    segments
  }
  else segs

  def length: Int = segments.length
  def apply(idx: Int): ReadSegment = segments(idx)
  def iterator: Iterator[ReadSegment] = segments.iterator
  override def toString: String = segments.map(_.toString).mkString
  lazy val totalBases: Int = segments.map(_.length).sum

  /** Splits the given bases into tuples with its associated read segment.  If strict is false then only return
    * tuples for which we have bases in `bases`, otherwise throw an exception.
    **/
  def structureRead(bases: String, strict: Boolean = true): Seq[SubRead] = {
    // TODO: TF: strict == false allows us to skip segments (in the filter)
    this.filter(s => strict || s.offset < bases.length).map {
      segment =>
        val (start, end) = segment.range(bases.length, strict=strict)
        SubRead(bases=bases.substring(start, end), segment=ReadSegment(segment, end - start))
    }
  }

  /** Splits the given bases and qualities into triples with its associated read segment.  If strict is false then only
    * return tuples for which we have bases in `bases`, otherwise throw an exception.
    **/
  def structureReadWithQualities(bases: String, qualities: String, strict: Boolean = true): Seq[SubRead] = {
    assert(bases.length == qualities.length)
    this.filter(s => strict || s.offset < bases.length).map {
      segment =>
        val (start, end) = segment.range(Math.min(bases.length, qualities.length), strict=strict)
        SubRead(bases=bases.substring(start, end), quals=Some(qualities.substring(start,end)), segment=ReadSegment(segment, end - start))
    }
  }

  private def collectSegments[T <: ReadSegment: ClassTag]: Seq[T] = this.collect { case t: T => t }
  def template: Seq[Template] = this.collectSegments[Template]
  def sampleBarcode: Seq[SampleBarcode] = this.collectSegments[SampleBarcode]
  def molecularBarcode: Seq[MolecularBarcode] = this.collectSegments[MolecularBarcode]
  def skip: Seq[Skip] = this.collectSegments[Skip]
}

object ReadSegment {
  val Types = Seq('T', 'B', 'M', 'S')

  def apply(offset: Int, length: Int, s: String): ReadSegment = {
    if (s.isEmpty) throw new IllegalArgumentException("Empty string given")
    apply(offset, length, s.head)
  }

  def apply(offset: Int, length: Int, c: Char): ReadSegment = {
    c match {
      case 'T' => Template(offset, length)
      case 'B' => SampleBarcode(offset, length)
      case 'M' => MolecularBarcode(offset, length)
      case 'S' => Skip(offset, length)
      case _ => throw new IllegalArgumentException(s"Unknown read segment type: $c")
    }
  }

  /** Creates a new read segment of the same type but with a new length */
  def apply(segment: ReadSegment, length: Int): ReadSegment = segment match {
    case seg: Template         => Template(seg.offset, length)
    case seg: SampleBarcode    => SampleBarcode(seg.offset, length)
    case seg: MolecularBarcode => MolecularBarcode(seg.offset, length)
    case seg: Skip             => Skip(seg.offset, length)
  }

  /** Creates a new read segment of the same type but with a new offset and length */
  def apply(segment: ReadSegment, offset: Int, length: Int): ReadSegment = segment match {
    case seg: Template         => Template(offset, length)
    case seg: SampleBarcode    => SampleBarcode(offset, length)
    case seg: MolecularBarcode => MolecularBarcode(offset, length)
    case seg: Skip             => Skip(offset, length)
  }
}

/** Base for all types of segments of a read */
sealed trait ReadSegment {
  /** The number of bases in the read preceding this segment. */
  def offset: Int
  /** The length of this segment. */
  def length: Int
  /** The single-character symbol representing this segment. */
  def symbol: Char

  /** Gets first and last bases (exclusive, 0-based) associated with this read segment.  If strict is false then only
    * bound the last base to `basesLen`, otherwise if the end is too big, throw an exception.
    **/
  private[util] def range(basesLen: Int, strict: Boolean): (Int, Int) = {
    val start = this.offset
    val end = length + start
    if (basesLen < end) {
      if (strict) throw new IllegalStateException("Not enough bases for read segment: " + this)
      else (start, basesLen)
    }
    else (start, end)
  }

  /** Gets the bases associated with this read segment.  If strict is false then only return
    * the sub-sequence for which we have bases in `bases`, otherwise throw an exception.
    **/
  private[util] def extractBases(bases: String, strict: Boolean = true): String = {
    val (start, end) = range(bases.length, strict=strict)
    bases.substring(start, end)
  }

  /** Gets the bases and qualities associated with this read segment.  If strict is false then only return
    * the sub-sequence for which we have bases in `bases`, otherwise throw an exception.
    **/
  private[util] def extractBasesAndQualities(bases: String, qualities: String, strict: Boolean = true): (String, String) = {
    val (start, end) = range(Math.min(bases.length, qualities.length), strict=strict)
    (bases.substring(start, end), qualities.substring(start, end))
  }

  /** Provides a string representation of this segment (ex. "23T" or "4M"). */
  override def toString: String = s"$length$symbol"
}

/** Bases from the template */
case class Template(offset: Int, length: Int, symbol: Char = 'T') extends ReadSegment

/** Bases used to identify the sample */
case class SampleBarcode(offset: Int, length: Int, symbol: Char = 'B') extends ReadSegment

/** Bases used to identify the source molecule */
case class MolecularBarcode(offset: Int, length: Int, symbol: Char = 'M') extends ReadSegment

/** Bases to ignore */
case class Skip(offset: Int, length: Int, symbol: Char = 'S') extends ReadSegment
