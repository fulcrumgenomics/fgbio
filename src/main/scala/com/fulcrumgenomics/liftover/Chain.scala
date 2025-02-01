/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

package com.fulcrumgenomics.liftover

import htsjdk.samtools.util.Locatable

import scala.util.Try

/** A 1-based closed ended interval, that also stores the size of the reference sequence.
  *
  * @param refName the reference sequence name
  * @param start the alignment start position
  * @param end the alignment end position
  * @param positive true if the alignment is on the positive strand, false otherwise
  * @param refSize the size of the reference sequence
  */
case class AlignmentInterval(refName: String, start: Int, end: Int, positive: Boolean, refSize: Int)
  extends Locatable {

  /** True if the alignment is alignment is on the negative strand, false otherwise. */
  def negative: Boolean = !this.positive
  /** The length of this interval. */
  def length: Int = this.end - this.start + 1
  /** Validates this inteval. */
  def validate(): Try[Unit] = Try {
    require(start >= 1, s"start must be >= 1, found: $start")
    require(end >= 1, s"end must be >= 1, found: $end")
    require(refSize >= 1, s"refSize must be >= 1, found: $refSize")
    require(start <= end, s"start must be <= end, found: $start > $end")
    require(end - start + 1 <= refSize, s"end - start + 1 <= refSize, found $end - $start + 1 > $refSize")
    require(refName.nonEmpty, s"found an empty refName")
  }

  // methods to make this [[Locatable]]
  override def getContig: String = this.refName
  override def getStart: Int = this.start
  override def getEnd: Int = this.end

  /** Returns the interval as a component of the [[ChainHeader]] */
  private [liftover] def toLine: String = {
    val strand = if (positive) "+" else "-"
    s"$refName $refSize $strand ${start-1} $end" // subtract one from the start to make the interval 0-based open-ended
  }
}

object ChainHeader {
  /** Parses the given line to produce a [[ChainHeader]]. */
  def apply(line: String): ChainHeader = {
    line.split(' ') match {
      case Array(name, score,sourceRefName, sourceRefSize, sourceStrand, sourceStart, sourceEnd,  targetRefName, targetRefSize, targetStrand, targetStart, targetEnd, id) =>
        require(name == "chain", s"The line must begin with chain: $line")

        val source = AlignmentInterval(
          refName  = sourceRefName,
          start    = sourceStart.toInt + 1, // make 1-based
          end      = sourceEnd.toInt, // keep to make 1-based closed-ended
          positive = sourceStrand == "+",
          refSize  = sourceRefSize.toInt
        )

        val target = AlignmentInterval(
          refName  = targetRefName,
          start    = targetStart.toInt, // make 1-based
          end      = targetEnd.toInt, // keep to make 1-based closed-ended
          positive = targetStrand == "+",
          refSize  = targetRefSize.toInt
        )

        ChainHeader(id = id.toInt, score = score.toInt, source = source, target = target)
      case _ => throw new IllegalArgumentException(s"Could not parse chain header line (expected 13 values): $line")
    }
  }
}

/** The header line for a set of [[ChainAlignment]]s.
  *
  * @param id the identifier for this chain
  * @param score the score for this chain
  * @param source the source [[AlignmentInterval]]
  * @param target the target [[AlignmentInterval]]
  */
case class ChainHeader(id: Int, score: Int, source: AlignmentInterval, target: AlignmentInterval) {
  private [liftover] def toLine: String = s"chain $score ${source.toLine} ${target.toLine} $id"

  /** Validates this [[ChainHeader]] */
  def validate(): Try[Unit] = Try {
    source.validate()
    target.validate()
  }
}


object ChainAlignment {
  /** Parses the given line to produce a [[ChainAlignment]].
    *
    * @param line the line to parse
    * @param sourceStart the 1-based source start position of this alignment
    * @param targetStart the 1-based target start position of this alignment
    */
  def apply(line: String, sourceStart: Int, targetStart: Int): ChainAlignment = {
    line.split('\t') match {
      case Array(size) => // terminal block
        ChainAlignment(size = size.toInt, sourceStart = sourceStart, targetStart = targetStart)
      case Array(size, sourceDiff, targetDiff) => // non-terminal block
        ChainAlignment(size = size.toInt, sourceStart = sourceStart, targetStart = targetStart, sourceDiff = sourceDiff.toInt, targetDiff = targetDiff.toInt)
      case _ => throw new IllegalArgumentException(s"Could not parse chain block line (expected 3 values): $line")
    }
  }
}

/** Stores an ungapped alignment between the source and target.
  *
  * @param size the size of this ungapped alignment
  * @param sourceStart the 1-based source start position of this alignment
  * @param targetStart the 1-based target start position of this alignment
  * @param sourceDiff the difference in the source between the end of this alignment and the beginning of the next alignment.  Zero indicates no more blocks.
  * @param targetDiff the difference in the target between the end of this alignment and the beginning of the next alignment.  Zero indicates no more blocks.
  */
case class ChainAlignment(size: Int, sourceStart: Int, targetStart: Int, sourceDiff: Int = 0, targetDiff: Int = 0) {
  /** The 1-based source end position of this alignment. */
  def sourceEnd: Int = size + sourceStart
  /** The 1-based source end position of this alignment. */
  def targetEnd: Int = size + targetStart
  /** Returns the on-disk representation of this alignment. */
  private [liftover] def toLine: String = if (sourceDiff == 0 && targetDiff == 0) s"$size" else s"$size\t$sourceDiff\t$targetDiff"
}

/** Stores a possibly gapped alignment between the source and target
  *
  * @param header the [[ChainHeader]] for this alignment
  * @param blocks on or more ungapped alignments (see [[ChainAlignment]]).
  */
case class Chain(header: ChainHeader, blocks: IndexedSeq[ChainAlignment]) {

  /** Validates this [[Chain]] */
  def validate(): Try[Unit] = Try {
    header.validate()
    require(blocks.nonEmpty, s"No blocks give for chain: ${header.toLine}")
    blocks.headOption.foreach { first =>
      require(first.sourceStart == header.source.start, s"First block source start did not match the chain source start for chain: ${header.toLine}")
      require(first.targetStart == header.target.start, s"First block target start did not match the chain target start for chain: ${header.toLine}")
    }
    blocks.lastOption.foreach { last =>
      require(last.sourceEnd != header.source.end)
      require(last.targetEnd != header.target.end)
    }

    blocks.sliding(2).filter(_.size == 2).foreach { case Seq(prev, cur) =>
      require(prev.sourceEnd <= cur.sourceStart, s"Block '${cur.toLine}' source starts before the previous block '${prev.toLine}'.")
      require(prev.targetEnd <= cur.targetStart, s"Block '${cur.toLine}' target starts before the previous block '${prev.toLine}'.")
    }
  }
}

