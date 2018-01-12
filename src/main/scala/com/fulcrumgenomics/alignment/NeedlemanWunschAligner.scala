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
import com.fulcrumgenomics.alignment.Mode.{Global, Glocal}
import com.fulcrumgenomics.alignment.NeedlemanWunschAligner._
import enumeratum.EnumEntry
import htsjdk.samtools.CigarOperator

import scala.collection.{immutable, mutable}
import scala.math.max

/** Trait that entires in Mode will extend. */
sealed trait Mode extends EnumEntry

/** Enum to represent alignment modes supported by the NeedlenameWunschAligner */
object Mode extends FgBioEnum[Mode] {
  /** Alignment mode for global pairwise alignment. */
  case object Global extends Mode
  /** Alignment mode for global alignment of query sequence to local region of target sequence. */
  case object Glocal extends Mode

  override def values: immutable.IndexedSeq[Mode] = findValues
}


object NeedlemanWunschAligner {
  /** Generates a simple scoring function using the match and mismatch scores. */
  def simpleScoringFunction(matchScore: Int, mismatchScore: Int): (Byte, Byte) => Int = {
    (lhs: Byte, rhs: Byte) => if (lhs == rhs) matchScore else mismatchScore
  }

  /** Creates a NW aligner with fixed match and mismatch scores. */
  def apply(matchScore: Int, mismatchScore: Int, gapOpen: Int, gapExtend: Int, mode: Mode = Global): NeedlemanWunschAligner = {
    new NeedlemanWunschAligner(scoringFunction=simpleScoringFunction(matchScore, mismatchScore), gapOpen, gapExtend, mode=mode)
  }

  /** Directions within the trace back matrix. */
  type Direction = Byte
  val Left: Direction     = 1.toByte
  val Up  : Direction     = 2.toByte
  val Diagonal: Direction = 3.toByte
  val Done: Direction     = 4.toByte
}

/**
  * Implementation of the Needleman-Wunsch algorithm for global alignment of two sequences with support for an
  * affine gap penalty.
  *
  * @param scoringFunction a function to score the alignment of a pair of bases
  * @param gapOpen the gap opening penalty
  * @param gapExtend the gap extension penalty
  * @param useEqualsAndX if true use the = and X cigar operators for matches and mismatches,
  *                      else use the M operator for both.
  * @param mode alignment mode to use when generating alignments
  */
class NeedlemanWunschAligner(val scoringFunction: (Byte,Byte) => Int,
                             val gapOpen: Int,
                             val gapExtend: Int,
                             useEqualsAndX: Boolean = true,
                             val mode: Mode = Global) {

  private val (matchOp, mismatchOp) = if (useEqualsAndX) (CigarOperator.EQ, CigarOperator.X) else (CigarOperator.M, CigarOperator.M)

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String): Alignment = align(query.getBytes, target.getBytes)

  /**
    * Align two sequences with the current scoring system.  If the [[Mode]] is `Global` the query and target
    * may be of any length.  If the [[Mode]] is Glocal then the target must be at least as long as the query.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Alignment]] object describing the optimal global alignment of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte]): Alignment = {
    require(query.length <= target.length || this.mode == Mode.Global, "For glocal, query length must be <= target length.")
    val (scoring, trace) = buildMatrices(query, target)
    generateAlignment(query, target, scoring, trace)
  }

  /**
    * Constructs both the scoring and traceback matrices.
    *
    * Matrix is constructed such that:
    *   - Matrix as nRow * nColumn (i.e. coordinates are (row, column))
    *   - Query sequence is down the side (each row represents a query base)
    *   - Target sequence is along the top (each column represents a target base)
    *   - i is iterating over rows/query bases
    *   - j is iterating over columns/target bases
    *   - Up   == Consumes query base  but not target base, i.e. Insertion vs. the target
    *   - Left == Consumes target base but not query base, i.e. Deletion vs. the target
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return a tuple of the scoring matrix and the traceback matrix
    */
  protected def buildMatrices(query: Array[Byte], target: Array[Byte]): (Matrix[Int], Matrix[Direction]) = {
    val scoring = Matrix[Int](query.length + 1, target.length + 1)
    val trace   = Matrix[Direction](query.length + 1, target.length + 1)

    // Top left corner
    scoring(0,0) = 0
    trace(0,0)   = Done

    // Down the left hand side
    forloop(from=1, until=query.length+1) { i =>
      trace(i, 0)   = Up
      scoring(i, 0) = scoring(i-1, 0) + gapExtend + (if (i == 1) gapOpen else 0)
    }

    // Along the top - for global we have to keep going left, for glocal we're done
    mode match {
      case Global =>
        forloop(from=1, until=target.length+1) { j =>
          trace(0, j)   = Left
          scoring(0, j) = scoring(0, j-1) + gapExtend + (if (j == 1) gapOpen else 0)
        }
      case Glocal =>
        forloop(from=1, until=target.length+1) { j =>
          trace(0, j)   = Done
          scoring(0, j) = 0
        }
    }

    // The interior of the matrix
    forloop(from=1, until=query.length+1) { i =>
      forloop(from=1, until=target.length+1) { j =>
        val diagonalScore = scoring(i-1, j-1) + scoringFunction(query(i-1), target(j-1))
        val upScore       = scoring(i-1, j  ) + gapExtend + (if (trace(i-1, j  ) == Up) 0 else gapOpen)
        val leftScore     = scoring(i,   j-1) + gapExtend + (if (trace(i,   j-1) == Left) 0 else gapOpen)

        val maxScore = max(max(diagonalScore, upScore), leftScore)
        scoring(i,j) = maxScore

        // Prefer mismatches first when back-tracing so as to push indels to the left of the alignment
        // The choice to prefer deletions before insertions is largely arbitrary!
        trace(i,j)   = if (maxScore == diagonalScore) Diagonal else if (maxScore == leftScore) Left else Up
      }
    }

    (scoring, trace)
  }

  /**
    * Given the scoring and trace back matrices, construct the [[Alignment]] object for the (or one of the)
    * optimal alignments.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @param scoring the scoring matrix
    * @param trace the trace back matrix
    * @return an [[Alignment]] object representing the alignment
    */
  protected def generateAlignment(query: Array[Byte], target: Array[Byte], scoring: Matrix[Int], trace: Matrix[Direction]): Alignment = {
    var currOperator: CigarOperator = null
    var currLength: Int = 0
    val elems = new mutable.ArrayBuffer[CigarElem]()

    // Always start at the end of the query
    var i = query.length

    // For the target sequence, start at the end for Global, or the highest scoring cell in the last row for Glocal
    var j = this.mode match {
      case Global => target.length
      case Glocal =>
        var (maxScore, maxIndex) = (Int.MinValue, -1)
        forloop(from=0, until=target.length) { i =>
          val score = scoring(query.length-1, i)
          if (score > maxScore) {
            maxScore = score
            maxIndex = i
          }
        }

        maxIndex + 1
    }

    // The score is always the score from the starting cell
    val score = scoring(i,j)

    // For global we have to reach the origin, for glocal we just have to reach the top row
    val done = this.mode match {
      case Global => () => i == 0 && j == 0
      case Glocal => () => i == 0
    }

    while (!done()) {
      val op = trace(i,j) match {
        case Up       =>
          i -= 1
          CigarOperator.INSERTION
        case Left     =>
          j -= 1
          CigarOperator.DELETION
        case Diagonal =>
          i -= 1
          j -= 1
          if (query(i) == target(j)) this.matchOp else this.mismatchOp
      }

      val extending = op == currOperator

      if (extending) {
        currLength += 1
      }
      else {
        if (currLength > 0) elems += CigarElem(currOperator, currLength)
        currOperator = op
        currLength   = 1
      }

      if (done()) elems += CigarElem(currOperator, currLength)
    }

    Alignment(query=query, target=target, queryStart=1, targetStart=j+1, cigar=Cigar(elems.reverse), score=score)
  }
}
