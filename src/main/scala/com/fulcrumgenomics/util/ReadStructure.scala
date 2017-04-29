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

import com.fulcrumgenomics.util.ReadStructure.{SubReadWithQuals, SubReadWithoutQuals}
import com.fulcrumgenomics.FgBioDef._

import scala.collection.immutable
import scala.collection.mutable.ArrayBuffer

/**
  * Companion object for ReadStructure that provides factory methods.
  */
object ReadStructure {
  /** A character that can be put in place of a number in a read structure to mean "0 or more bases". */
  val AnyLengthChar: Char = '+'

  /** Can be subtracted from any digit character to get it's integer value. */
  private val DigitOffset = '0'.toInt

  /** Contains the bases and optionally base qualities that correspond to the given read segment. */
  sealed trait SubRead {
    def bases: String
    def segment: ReadSegment
    def kind: SegmentType = this.segment.kind
  }
  case class SubReadWithoutQuals(bases: String, segment: ReadSegment) extends SubRead
  case class SubReadWithQuals(bases: String, quals: String, segment: ReadSegment) extends SubRead

  /** Creates a new ReadStructure, optionally resetting the offsets on each of the segments. */
  def apply(segments: Seq[ReadSegment], resetOffsets: Boolean = false): ReadStructure = {
    // Check that none but the last segment has an indefinite length
    require(segments.dropRight(1).forall(_.hasFixedLength), s"Variable length ($AnyLengthChar) can only be used in the last segment: ${segments.mkString}")

    val segs = if (!resetOffsets) segments else {
      var idx, off = 0
      segments.map { s => yieldAndThen(s.copy(offset=off)) { off += s.length.getOrElse(0) } }
    }

    segments.filter(_.length.exists(_ <= 0)) match {
      case Seq() => new ReadStructure(segs)
      case zeros => throw new IllegalArgumentException(s"Read structure contained zero length segments: ${segments.mkString}")
    }
  }

  /** Creates a read structure from a string */
  def apply(readStructure: String): ReadStructure = {
    ReadStructure(segments=segments(rs=readStructure.trim.toUpperCase), resetOffsets=true)
  }

  /** Creates a sequence of read segments from a string. */
  private def segments(rs: String): Seq[ReadSegment] = {
    var i = 0
    val segs = ArrayBuffer[ReadSegment]()
    while (i < rs.length) {
      // Stash the beginning position of our parsing so we can highlight what we're having trouble with
      val parsePosition = i

      // Parse out the length segment which many be 1 or more digits or the AnyLengthChar
      val segLength = rs.charAt(i) match {
        case AnyLengthChar  =>
          i += 1
          None
        case d if d.isDigit =>
          var len = 0
          while (i < rs.length && rs.charAt(i).isDigit) { len = (len*10) + rs.charAt(i) - DigitOffset; i += 1 }
          Some(len)
        case x =>
          invalid("Read structure missing length information", rs, parsePosition, parsePosition+1)
      }

      // Parse out the operator and make a segment
      if (i == rs.length) {
        invalid("Read structure with invalid segment", rs, parsePosition, i)
      }
      else {
        val code = rs.charAt(i)
        i += 1
        try   { segs += ReadSegment(0, segLength, SegmentType(code)) }
        catch { case ex: Exception => invalid("Read structure segment had unknown type", rs, parsePosition, i) }
      }
    }

    segs
  }

  /**
    * Inserts square brackets around the characters in the read structure that are causing the error.
    *
    * @param rs the read structure string
    * @param start the start of the error in the string (inclusive)
    * @param end the ned of the error in the string (exclusive)
    */
  private def invalid(msg: String, rs: String, start: Int, end: Int): Nothing = {
    val prefix = rs.substring(0, start)
    val error  = rs.substring(start, end)
    val suffix = if (end == rs.length) "" else rs.substring(end, rs.length)
    throw new IllegalArgumentException(s"$msg: $prefix[$error]$suffix")
  }
}

/**
  * Describes the structure of a give read.  A read contains one or more read segments. A read segment describes
  * a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
  * of the read.
  *
  * @param segments the segments composing this read structure.
  */
class ReadStructure private(val segments: Seq[ReadSegment]) extends immutable.Seq[ReadSegment] {
  require(segments.nonEmpty, "Can't create a read structure with zero segments.")

  /** The minimum length read that this read structure can process. */
  private val minLength = segments.flatMap(_.length).sum

  /** Returns true if the ReadStructure has a fixed (i.e. non-variable) length. */
  def hasFixedLength: Boolean = segments.last.hasFixedLength

  /** Returns the fixed length if there is one. Throws an exception on segments without fixed lengths! */
  def fixedLength: Int = {
    require(hasFixedLength, s"fixedLength called on variable length segment: $this")
    this.minLength
  }

  /** Length is defined as the number of segments (not bases!) in the read structure. */
  override def length: Int = segments.length

  /** Fetches the segment at the given index. */
  override def apply(idx: Int): ReadSegment = segments(idx)

  /** Provides an iterator over the segments. */
  override def iterator: Iterator[ReadSegment] = segments.iterator

  /** Generates the String format of the ReadStructure that can be re-parsed. */
  override def toString: String = segments.map(_.toString).mkString

  /** Generates a new ReadStucture that is the same as this one except that the last segment has undefined length. */
  def withVariableLastSegment: ReadStructure =
    if (this.segments.last.hasFixedLength) new ReadStructure(segments.dropRight(1) :+ segments.last.copy(length=None)) else this

  /** Splits the given bases into tuples with its associated read segment.  If strict is false then only return
    * tuples for which we have bases in `bases`, otherwise throw an exception.
    **/
  def extract(bases: String): Seq[SubReadWithoutQuals] = segments.map(_.extract(bases))

  /** Splits the given bases and qualities into triples with its associated read segment.  If strict is false then only
    * return tuples for which we have bases in `bases`, otherwise throw an exception.
    **/
  def extract(bases: String, quals: String): Seq[SubReadWithQuals] = {
    assert(bases.length == quals.length)
    segments.map(s => s.extract(bases, quals))
  }

  /** Returns just the segments of a given kind. */
  def segments(kind: SegmentType): Seq[ReadSegment] = this.segments.filter(_.kind == kind)

  def templateSegments:         Seq[ReadSegment] = segments(SegmentType.Template)
  def sampleBarcodeSegments:    Seq[ReadSegment] = segments(SegmentType.SampleBarcode)
  def molecularBarcodeSegments: Seq[ReadSegment] = segments(SegmentType.MolecularBarcode)
  def skipSegments:             Seq[ReadSegment] = segments(SegmentType.Skip)
}

/** Sealed class hierarchy for the types of segments that can show up in a read structure. */
sealed abstract class SegmentType(val code: Char) { final override def toString: String = String.valueOf(code) }
object SegmentType {
  case object Template extends SegmentType('T')
  case object SampleBarcode extends SegmentType('B')
  case object MolecularBarcode extends SegmentType('M')
  case object Skip extends SegmentType('S')

  /** All the possible types. */
  val values: Seq[SegmentType] = Seq(Template, SampleBarcode, MolecularBarcode, Skip)

  /** Returns the [[SegmentType]] for the given code/letter. */
  def apply(code: Char): SegmentType = code match {
    case 'T' => Template
    case 'B' => SampleBarcode
    case 'M' => MolecularBarcode
    case 'S' => Skip
    case _   => throw new IllegalArgumentException(s"Invalid read segment type: $code")
  }
}

object ReadSegment {
  val Types = Seq('T', 'B', 'M', 'S')

  /** Creates a new ReadSegment with undefined length (i.e. length=0 or more). */
  def apply(offset: Int, c: Char): ReadSegment = new ReadSegment(offset, None, kind=SegmentType(c))

  /** Constructs a ReadSegment with a definite length using Char c to find the segment type. */
  def apply(offset: Int, length: Int, c: Char): ReadSegment = new ReadSegment(offset=offset, length=Some(length), kind=SegmentType(c))
}

/**
  * Encapsulates all the information about a segment within a read structure. A segment can either
  * have a definite length, in which case length must be Some(Int), or an indefinite length (can be
  * any length, 0 or more) in which case length must be None.
  */
case class ReadSegment (offset: Int, length: Option[Int], kind: SegmentType) {
  /** Checks some requirements and then calculates the end position for the segment for the given read. */
  private def calculateEnd(bases: String): Int = {
    require(bases.length >= offset, s"Read ends before read segment starts: $this")
    length.foreach(l => require(bases.length >= offset + l, s"Read ends before end of segment: $this"))

    if (hasFixedLength) math.min(offset + fixedLength, bases.length)
    else bases.length
  }

  /** Gets the bases associated with this read segment.  If strict is false then only return
    * the sub-sequence for which we have bases in `bases`, otherwise throw an exception.
    **/
  private[util] def extract(bases: String): SubReadWithoutQuals = {
    val end = calculateEnd(bases)
    SubReadWithoutQuals(bases=bases.substring(offset, end), segment=resized(end))
  }

  /** Gets the bases and qualities associated with this read segment.  If strict is false then only return
    * the sub-sequence for which we have bases in `bases`, otherwise throw an exception.
    **/
  private[util] def extract(bases: String, quals: String): SubReadWithQuals = {
    require(bases.length == quals.length, s"Bases and quals differ in length: $bases $quals")
    val end = calculateEnd(bases)
    SubReadWithQuals(bases=bases.substring(offset, end), quals=quals.substring(offset, end), segment=resized(end))
  }

  private def resized(end: Int): ReadSegment = {
    val newLength = end - offset
    if (length.contains(newLength)) this else this.copy(length=Some(newLength))
  }

  /** Returns true if the read segment has a defined length. */
  def hasFixedLength: Boolean = this.length.nonEmpty

  /** Returns the fixed length if there is one. Throws an exception on segments without fixed lengths! */
  def fixedLength: Int = this.length.getOrElse(throw new IllegalStateException(s"fixedLength called on variable length segment: $this"))

  /** Provides a string representation of this segment (ex. "23T" or "4M"). */
  override def toString: String = s"${length.getOrElse(ReadStructure.AnyLengthChar)}${kind.code}"
}
