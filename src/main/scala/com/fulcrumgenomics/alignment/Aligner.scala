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
import com.fulcrumgenomics.alignment.Mode.{Global, Glocal, Local}
import com.fulcrumgenomics.alignment.Aligner._
import enumeratum.EnumEntry
import htsjdk.samtools.CigarOperator

import scala.collection.{immutable, mutable}

/** Trait that entries in Mode will extend. */
sealed trait Mode extends EnumEntry

/** Enum to represent alignment modes supported by the Aligner. */
object Mode extends FgBioEnum[Mode] {
  /** Alignment mode for global pairwise alignment. */
  case object Global extends Mode
  /** Alignment mode for global alignment of query sequence to local region of target sequence. */
  case object Glocal extends Mode
  /** Alignment mode for local pairwise alignment. */
  case object Local extends Mode

  override def values: immutable.IndexedSeq[Mode] = findValues
}


object Aligner {
  /** Generates a simple scoring function using the match and mismatch scores. */
  def simpleScoringFunction(matchScore: Int, mismatchScore: Int): (Byte, Byte) => Int = {
    (lhs: Byte, rhs: Byte) => if (lhs == rhs) matchScore else mismatchScore
  }

  /** Creates a NW aligner with fixed match and mismatch scores. */
  def apply(matchScore: Int, mismatchScore: Int, gapOpen: Int, gapExtend: Int, mode: Mode = Global): Aligner = {
    new Aligner(scoringFunction=simpleScoringFunction(matchScore, mismatchScore), gapOpen, gapExtend, mode=mode)
  }

  /** Directions within the trace back matrix. */
  private type Direction = Byte
  private val Left: Direction     = 0.toByte
  private val Up  : Direction     = 1.toByte
  private val Diagonal: Direction = 2.toByte
  private val Done: Direction     = 3.toByte

  // NB: the order of LeftAndDiagonal and UpAndDiagonal matters when breaking ties!
  private val AllDirections: Seq[Direction]   = Seq(Diagonal, Left, Up)
  private val LeftAndUp: Seq[Direction]       = Seq(Left, Up)
  private val LeftAndDiagonal: Seq[Direction] = Seq(Diagonal, Left)
  private val UpAndDiagonal: Seq[Direction]   = Seq(Diagonal, Up)

  /** The minimum score allowed to start an alignment.  This prevents underflow. */
  val MinStartScore: Int = Int.MinValue / 2

  object AlignmentMatrix {
    def apply(direction: Direction, queryLength: Int, targetLength: Int): AlignmentMatrix = {
      AlignmentMatrix(direction=direction, scoring=Matrix[Int](queryLength+1, targetLength+1), trace=Matrix[Direction](queryLength+1, targetLength+1))
    }
  }

  /** A single alignment matrix for a given [[Direction]] storing both the scoring and traceback matrices produce by the aligner. */
  case class AlignmentMatrix(direction: Direction, scoring: Matrix[Int], trace: Matrix[Direction])
}


/**
  * Implementation of an aligner with generic scoring function and affine gap penalty support.
  * Supports multiple alignment [[Mode]]s for global, semi-global and local alignment.
  *
  * A scoring function (`scoringFunction`) is taken to score pair-wise aligned bases. A default
  * implementation is supplied via the companion [[Aligner]] object which uses a fixed match
  * score and mismatch penalty.
  *
  * @param scoringFunction a function to score the alignment of a pair of bases. Parameters
  *                        are the query and target bases as bytes, in that order.
  * @param gapOpen the gap opening penalty, should generally be negative or zero
  * @param gapExtend the gap extension penalty, should generally be negative or zero
  * @param useEqualsAndX if true use the = and X cigar operators for matches and mismatches,
  *                      else use the M operator for both.
  * @param mode alignment mode to use when generating alignments
  */
class Aligner(val scoringFunction: (Byte,Byte) => Int,
              val gapOpen: Int,
              val gapExtend: Int,
              useEqualsAndX: Boolean = true,
              val mode: Mode = Global) {

  private val (matchOp, mismatchOp) = if (useEqualsAndX) (CigarOperator.EQ, CigarOperator.X) else (CigarOperator.M, CigarOperator.M)

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String): Alignment = align(query.getBytes, target.getBytes)

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Alignment]] object describing the optimal global alignment of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte]): Alignment = {
    val matrices = buildMatrices(query, target)
    generateAlignment(query, target, matrices)
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
    * @return an array of alignment matrices, where the indices to the array are the Directions
    */
  protected def buildMatrices(query: Array[Byte], target: Array[Byte]): Array[AlignmentMatrix] = {
    val matrices = AllDirections.sorted.map(dir => AlignmentMatrix(direction=dir, queryLength=query.length, targetLength=target.length)).toArray

    // While we have `matrices` above, it's useful to unpack all the matrices for direct access
    // in the core loop; when we know the exact matrix we need at compile time, it's faster
    val (leftScoreMatrix, leftTraceMatrix, upScoreMatrix, upTraceMatrix, diagScoreMatrix, diagTraceMatrix) = {
      val Seq(l, u, d) = Seq(Left, Up, Diagonal).map(_.toInt).map(matrices.apply)
      (l.scoring, l.trace, u.scoring, u.trace, d.scoring, d.trace)
    }

    // Top left corner - allow all to be zero score but then must be careful to initialize Up and Left to have a gap open
    AllDirections.foreach { direction =>
      matrices(direction).scoring(0, 0) = 0
      matrices(direction).trace(0, 0)   = Done
    }

    // Down the left hand side - for global/glocal we have to keep going up, for local we're done
    mode match {
      case Global | Glocal =>
        forloop(from=1, until=query.length+1) { i =>
          LeftAndDiagonal.foreach { direction =>
            matrices(direction).scoring(i, 0) = MinStartScore
            matrices(direction).trace(i, 0)   = Done
          }
          matrices(Up).scoring(i, 0) = matrices(Up).scoring(i-1, 0) + gapExtend + (if (i == 1) gapOpen else 0)
          matrices(Up).trace(i, 0)   = if (i == 1) Diagonal else Up
        }
      case Local =>
        forloop(from=1, until=query.length+1) { i =>
          forloop(from=1, until=query.length+1) { i =>
            LeftAndUp.foreach { direction =>
              matrices(direction).scoring(i, 0) = MinStartScore
              matrices(direction).trace(i, 0)   = Done
            }
            matrices(Diagonal).scoring(i, 0) = 0
            matrices(Diagonal).trace(i, 0)   = Done
          }
        }
    }

    // Along the top - for global we have to keep going left, for glocal/local we're done
    mode match {
      case Global =>
        forloop(from=1, until=target.length+1) { j =>
          UpAndDiagonal.foreach { direction =>
            matrices(direction).scoring(0, j) = MinStartScore
            matrices(direction).trace(0 , j)   = Done
          }
          matrices(Left).scoring(0, j) = matrices(Left).scoring(0, j-1) + gapExtend + (if (j == 1) gapOpen else 0)
          matrices(Left).trace(0, j)   = if (j == 1) Diagonal else Left
        }
      case Glocal | Local =>
        forloop(from=1, until=target.length+1) { j =>
          LeftAndUp.foreach { direction =>
            matrices(direction).scoring(0, j) = MinStartScore
            matrices(direction).trace(0, j)   = Done
          }
          matrices(Diagonal).scoring(0, j) = 0
          matrices(Diagonal).trace(0, j)   = Done
        }
    }

    // The interior of the matrix
    var i = 1
    while (i <= query.length) {
      val queryBase = query(i-1)
      var j = 1

      while (j <= target.length) {
        { // Diagonal matrix can come from the previous diagonal, up, or left
          val addend = scoringFunction(queryBase, target(j-1))
          val dScore = diagScoreMatrix(i-1, j-1) + addend
          val lScore = leftScoreMatrix(i-1, j-1) + addend
          val uScore = upScoreMatrix(i-1, j-1)   + addend

          if (mode == Local && dScore < 0 && lScore < 0 && uScore < 0) {
            diagScoreMatrix(i, j) = 0
            diagTraceMatrix(i, j) = Done
          }
          else if (dScore >= uScore && dScore >= lScore) {
            diagScoreMatrix(i, j) = dScore
            diagTraceMatrix(i, j) = Diagonal
          }
          else if (lScore >= uScore) {
            diagScoreMatrix(i, j) = lScore
            diagTraceMatrix(i, j) = Left
          }
          else {
            diagScoreMatrix(i, j) = uScore
            diagTraceMatrix(i, j) = Up
          }
        }


        { // Up matrix can come from diagonal or up
          val dScore = gapOpen + diagScoreMatrix(i-1, j)
          val uScore =           upScoreMatrix(i-1, j)
          if (dScore >= uScore) {
            upScoreMatrix(i, j) = dScore + gapExtend
            upTraceMatrix(i, j) = Diagonal
          }
          else {
            upScoreMatrix(i, j) = uScore + gapExtend
            upTraceMatrix(i, j) = Up
          }
        }


        { // Left matrix can come from diagonal or left
          val dScore = gapOpen + diagScoreMatrix(i, j-1)
          val lScore =           leftScoreMatrix(i, j-1)
          if (dScore >= lScore) {
            leftScoreMatrix(i, j) = dScore + gapExtend
            leftTraceMatrix(i, j) = Diagonal
          }
          else {
            leftScoreMatrix(i, j) = lScore + gapExtend
            leftTraceMatrix(i, j) = Left
          }
        }

        j += 1
      }

      i += 1
    }

    matrices
  }

  /**
    * Given the scoring and trace back matrices, construct the [[Alignment]] object for the (or one of the)
    * optimal alignments.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @param matrices the scoring and trace back matrices for the [[Left]], [[Up]], and [[Diagonal]] directions.
    * @return an [[Alignment]] object representing the alignment
    */
  protected def generateAlignment(query: Array[Byte], target: Array[Byte], matrices: Array[AlignmentMatrix]): Alignment = {
    var currOperator: CigarOperator = null
    var currLength: Int = 0
    val elems = new mutable.ArrayBuffer[CigarElem]()

    // For the target sequence, start at the end for Global, or the highest scoring cell in the last row for Glocal, or
    // the highest scoring cell anywhere in the matrix for Local
    var (curI, curJ, curD) = this.mode match {
      case Global =>
        (query.length, target.length, AllDirections.maxBy(d => matrices(d).scoring(query.length, target.length)))
      case Glocal =>
        var (maxScore, maxJ, maxD) = (MinStartScore, -1, Done)
        forloop(from=1, until=target.length+1) { j =>
          val direction = AllDirections.maxBy(d => matrices(d).scoring(query.length, j))
          val score     = matrices(direction).scoring(query.length, j)
          if (score > maxScore) {
            maxScore = score
            maxJ = j
            maxD = direction
          }
        }
        (query.length, maxJ, maxD)
      case Local =>
        var (maxScore, maxI, maxJ, maxD) = (MinStartScore, -1, -1, Done)
        forloop(from=1, until=query.length+1) { i =>
          forloop(from=1, until=target.length+1) { j =>
            val direction = AllDirections.maxBy(d => matrices(d).scoring(i, j))
            val score     = matrices(direction).scoring(i, j)
            if (score > maxScore) {
              maxScore = score
              maxI = i
              maxJ = j
              maxD = direction
            }
          }
        }
        (maxI, maxJ, maxD)
      }

    // The score is always the score from the starting cell
    val score = matrices(curD).scoring(curI, curJ)

    // For global we have to reach the origin, for glocal we just have to reach the top row, and for local we can stop
    // anywhere.  Fortunately, we have initialized the appropriate cells to "Done"
    val done = () => matrices(curD).trace(curI,curJ) == Done
    while (!done()) {
      val nextD = matrices(curD).trace(curI,curJ)
      val op    = curD match {
        case Up       =>
          curI -= 1
          CigarOperator.INSERTION
        case Left     =>
          curJ -= 1
          CigarOperator.DELETION
        case Diagonal =>
          curI -= 1
          curJ -= 1
          if (query(curI) == target(curJ)) this.matchOp else this.mismatchOp
      }
      curD = nextD

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

    Alignment(query=query, target=target, queryStart=curI+1, targetStart=curJ+1, cigar=Cigar(elems.reverse), score=score)
  }
}
