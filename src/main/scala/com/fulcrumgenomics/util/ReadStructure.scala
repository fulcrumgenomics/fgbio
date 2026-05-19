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

import scala.collection.immutable

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

  /** Creates a new ReadStructure from a sequence of segments. At most one segment may be the
    * indefinite-length (`+`) segment, and it may appear at any position. */
  def apply(segments: Seq[ReadSegment]): ReadStructure = {
    require(segments.nonEmpty, "Can't create a read structure with zero segments.")
    val indefiniteCount = segments.count(!_.hasFixedLength)
    require(
      indefiniteCount <= 1,
      s"At most one segment may have indefinite length ($AnyLengthChar): ${segments.mkString}"
    )
    segments.find(_.length.exists(_ <= 0)) match {
      case Some(bad) =>
        throw new IllegalArgumentException(s"Read structure contained a zero or negative length segment: $bad in ${segments.mkString}")
      case None =>
        new ReadStructure(segments.toIndexedSeq)
    }
  }

  /** Creates a read structure from its string form (e.g. `"8B+M10T"`), in which each segment is a
    * length (either an integer or `+` for indefinite length) followed by a one-character segment-type
    * code.  Whitespace is ignored. */
  def apply(readStructure: String): ReadStructure = {
    val tidied = readStructure.toUpperCase.filterNot(Character.isWhitespace)
    ReadStructure(segments(rs=tidied))
  }

  /** Creates a sequence of read segments from a string. */
  private def segments(rs: String): Seq[ReadSegment] = {
    var i = 0
    val segs = IndexedSeq.newBuilder[ReadSegment]
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
        case _ =>
          invalid("Read structure missing length information", rs, parsePosition, parsePosition+1)
      }

      // Parse out the operator and make a segment
      if (i == rs.length) {
        invalid("Read structure with invalid segment", rs, parsePosition, i)
      }
      else {
        val code = rs.charAt(i)
        i += 1
        try   { segs += ReadSegment(segLength, SegmentType(code)) }
        catch { case _: Exception => invalid("Read structure segment had unknown type", rs, parsePosition, i) }
      }
    }

    segs.result()
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
  * Describes the structure of a given read.  A read contains one or more read segments. A read segment
  * describes a contiguous stretch of bases of the same type (e.g. template bases) and some length.
  *
  * At most one segment may be the indefinite-length (`+`) segment, meaning "whatever is left of the read".
  * The `+` segment may appear at any position, not only at the end.
  *
  * @param segments the segments composing this read structure.
  */
class ReadStructure private(val segments: IndexedSeq[ReadSegment]) extends immutable.Seq[ReadSegment] {
  import ReadStructure.AnyLengthChar

  /** The index of the indefinite-length segment, if any. */
  private val plusIndex: Option[Int] = segments.indexWhere(!_.hasFixedLength) match {
    case -1 => None
    case i  => Some(i)
  }

  /** Sum of lengths across all fixed-length segments. */
  private val fixedLengthSum: Int = segments.iterator.flatMap(_.length).sum

  /** Sum of fixed-length segment lengths strictly AFTER the `+` segment (zero when there is no `+`
    * or the `+` is trailing). Used to resolve the end of the `+` segment at extract time. */
  private val postPlusLen: Int = plusIndex match {
    case Some(p) => segments.drop(p + 1).iterator.flatMap(_.length).sum
    case None    => 0
  }

  /** Signed per-segment start offset.
    *
    *  - `>= 0` — offset from the start of the read. Used for segments before or at the `+`, and for
    *    every segment when there is no `+`.
    *  - `<  0` — distance from the end of the read, stored as a negative number. Used for segments
    *    strictly after a non-terminal `+`. The actual start position in a read of length `L` is
    *    `L + offset` (i.e. `L - (-offset)`).
    */
  private val offsets: IndexedSeq[Int] = {
    val n   = segments.length
    val out = Array.fill[Int](n)(0)

    // Forward pass up to and including the `+` (or the whole vec if no `+`).
    val forwardEnd = plusIndex.map(_ + 1).getOrElse(n)
    var off        = 0
    var i          = 0
    while (i < forwardEnd) {
      out(i) = off
      off += segments(i).length.getOrElse(0)
      i += 1
    }

    // Backward pass for segments strictly after the `+`, if any.
    plusIndex.foreach { p =>
      var distFromEnd = 0
      var j           = n - 1
      while (j > p) {
        val len = segments(j).length.getOrElse(
          throw new IllegalStateException(s"Post-$AnyLengthChar segment must have a fixed length; got ${segments(j)}")
        )
        distFromEnd += len
        out(j)      = -distFromEnd
        j          -= 1
      }
    }

    out.toIndexedSeq
  }

  /** Returns true if the ReadStructure has a fixed (i.e. non-variable) length. */
  def hasFixedLength: Boolean = plusIndex.isEmpty

  /** Returns the fixed length if there is one. Throws an exception on structures with an indefinite segment! */
  def fixedLength: Int = {
    require(hasFixedLength, s"fixedLength called on variable length structure: $this")
    fixedLengthSum
  }

  /** Length is defined as the number of segments (not bases!) in the read structure. */
  override def length: Int = segments.length

  /** Fetches the segment at the given index. */
  override def apply(idx: Int): ReadSegment = segments(idx)

  /** Provides an iterator over the segments. */
  override def iterator: Iterator[ReadSegment] = segments.iterator

  /** Generates the String format of the ReadStructure that can be re-parsed. */
  override def toString: String = segments.iterator.map(_.toString).mkString

  /** Generates a new ReadStructure that is the same as this one except that the last segment has undefined length.
    * If the structure already has a `+` segment (in any position), it is returned unchanged, since it already
    * handles reads of variable length. */
  def withVariableLastSegment: ReadStructure =
    if (this.hasFixedLength) ReadStructure(segments.dropRight(1) :+ segments.last.copy(length = None)) else this

  /** Validates that a read of the given length can be extracted by this structure.
    * For fixed-length structures we require an exact length match; silent truncation of trailing
    * bases is almost always a bug (wrong structure, stray adapter, etc.), so we require the caller
    * to either supply an exact-length read or convert the structure via [[withVariableLastSegment]]. */
  private def validateReadLength(readLen: Int): Unit = {
    if (plusIndex.isDefined) {
      require(
        readLen >= fixedLengthSum,
        s"Read length $readLen is shorter than required $fixedLengthSum for read structure $this"
      )
    }
    else {
      require(
        readLen == fixedLengthSum,
        s"Read length $readLen does not match fixed read structure length $fixedLengthSum for $this"
      )
    }
  }

  /** Returns the `[start, end)` span within a read of length `readLen` that corresponds to the
    * segment at index `i`. */
  private def spanOf(i: Int, readLen: Int): (Int, Int) = {
    if (plusIndex.contains(i)) {
      // The indefinite-length segment: runs from its stored start offset to just before the post-`+` region.
      (offsets(i), readLen - postPlusLen)
    }
    else {
      val off   = offsets(i)
      val start = if (off >= 0) off else readLen + off
      (start, start + segments(i).length.get)
    }
  }

  /** Extracts the bases corresponding to each read segment.  [[SegmentType.Skip]] segments are omitted
    * from the returned sequence because callers almost always throw those bases away; use the two-argument
    * overload with `includeSkips = true` to keep them. */
  def extract(bases: String): Seq[SubReadWithoutQuals] = extract(bases, includeSkips = false)

  /** Extracts the bases corresponding to each read segment.
    *
    * @param bases the raw read bases
    * @param includeSkips whether to include [[SegmentType.Skip]] segments in the returned sequence
    */
  def extract(bases: String, includeSkips: Boolean): Seq[SubReadWithoutQuals] = {
    validateReadLength(bases.length)
    val readLen = bases.length
    val builder = IndexedSeq.newBuilder[SubReadWithoutQuals]
    var i       = 0
    while (i < segments.length) {
      val seg = segments(i)
      if (includeSkips || seg.kind != SegmentType.Skip) {
        val (start, end) = spanOf(i, readLen)
        builder += SubReadWithoutQuals(bases.substring(start, end), seg)
      }
      i += 1
    }
    builder.result()
  }

  /** Extracts the bases and qualities corresponding to each read segment.
    *
    * @param bases the raw read bases
    * @param quals the raw read qualities; must be the same length as `bases`
    * @param includeSkips whether to include [[SegmentType.Skip]] segments in the returned sequence;
    *                     defaults to `false` because callers almost always throw those bases away.
    */
  def extract(bases: String, quals: String, includeSkips: Boolean = false): Seq[SubReadWithQuals] = {
    require(bases.length == quals.length, s"Bases and quals differ in length: ${bases.length} vs ${quals.length}")
    validateReadLength(bases.length)
    val readLen = bases.length
    val builder = IndexedSeq.newBuilder[SubReadWithQuals]
    var i       = 0
    while (i < segments.length) {
      val seg = segments(i)
      if (includeSkips || seg.kind != SegmentType.Skip) {
        val (start, end) = spanOf(i, readLen)
        builder += SubReadWithQuals(bases.substring(start, end), quals.substring(start, end), seg)
      }
      i += 1
    }
    builder.result()
  }

  /** Returns just the segments of a given kind. */
  def segments(kind: SegmentType): Seq[ReadSegment] = this.segments.filter(_.kind == kind)

  def templateSegments:         Seq[ReadSegment] = segments(SegmentType.Template)
  def sampleBarcodeSegments:    Seq[ReadSegment] = segments(SegmentType.SampleBarcode)
  def molecularBarcodeSegments: Seq[ReadSegment] = segments(SegmentType.MolecularBarcode)
  def cellBarcodeSegments:      Seq[ReadSegment] = segments(SegmentType.CellBarcode)
  def skipSegments:             Seq[ReadSegment] = segments(SegmentType.Skip)
}

/** Sealed class hierarchy for the types of segments that can show up in a read structure. */
sealed abstract class SegmentType(val code: Char) { final override def toString: String = String.valueOf(code) }
object SegmentType {
  case object Template extends SegmentType('T')
  case object SampleBarcode extends SegmentType('B')
  case object MolecularBarcode extends SegmentType('M')
  case object CellBarcode extends SegmentType('C')
  case object Skip extends SegmentType('S')

  /** All the possible types. */
  val values: Seq[SegmentType] = Seq(Template, SampleBarcode, MolecularBarcode, CellBarcode, Skip)

  /** Returns the [[SegmentType]] for the given code/letter. */
  def apply(code: Char): SegmentType = code match {
    case 'T' => Template
    case 'B' => SampleBarcode
    case 'M' => MolecularBarcode
    case 'C' => CellBarcode
    case 'S' => Skip
    case _   => throw new IllegalArgumentException(s"Invalid read segment type: $code")
  }
}

/**
  * Encapsulates all the information about a segment within a read structure. A segment can either
  * have a definite length, in which case length must be Some(Int), or an indefinite length (can be
  * any length, 0 or more) in which case length must be None.  The position of a segment's bases
  * within a read is tracked by the enclosing [[ReadStructure]], not here.
  */
case class ReadSegment(length: Option[Int], kind: SegmentType) {
  /** Returns true if the read segment has a defined length. */
  def hasFixedLength: Boolean = this.length.nonEmpty

  /** Returns the fixed length if there is one. Throws an exception on segments without fixed lengths! */
  def fixedLength: Int = this.length.getOrElse(throw new IllegalStateException(s"fixedLength called on variable length segment: $this"))

  /** Provides a string representation of this segment (ex. "23T" or "4M"). */
  override def toString: String = s"${length.getOrElse(ReadStructure.AnyLengthChar)}${kind.code}"
}

object ReadSegment {
  /** Constructs a [[ReadSegment]] from a length and a segment type character. */
  def apply(length: Int, c: Char): ReadSegment = ReadSegment(length = Some(length), kind = SegmentType(c))
}
