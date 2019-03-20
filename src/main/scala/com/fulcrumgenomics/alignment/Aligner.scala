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
import com.fulcrumgenomics.util.Sequences
import enumeratum.EnumEntry
import htsjdk.samtools.CigarOperator

import scala.collection.{immutable, mutable}
import scala.annotation.switch
import scala.collection.mutable.ArrayBuffer

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
  /** Creates a NW aligner with fixed match and mismatch scores. */
  def apply(matchScore: Int, mismatchScore: Int, gapOpen: Int, gapExtend: Int, mode: Mode = Global): Aligner = {
    val scorer: AlignmentScorer = new AlignmentScorer {
      private val matching      = matchScore
      private val mismatch      = mismatchScore
      private val open          = gapOpen
      private val extend        = gapExtend
      private val openAndExtend = open + extend

      override final def scorePairing(queryBase: Byte, targetBase: Byte): Int = {
        if (queryBase == targetBase) this.matching else this.mismatch
      }
      override final def scoreGap(query: Array[Byte], target: Array[Byte], qOffset: Int, tOffset: Int, inQuery: Boolean, extend: Boolean): Int = {
        if (extend) this.extend else this.openAndExtend
      }
    }

    new Aligner(scorer, mode=mode)
  }

  /** Directions within the trace back matrix. */
  private type Direction = Int
  // NB: types are purposely omitted for match performance.  See: https://stackoverflow.com/questions/16311540/why-cant-scala-optimize-this-match-to-a-switch
  private final val Left     = 0
  private final val Up       = 1
  private final val Diagonal = 2
  private final val Done     = 3

  private val AllDirections: Seq[Direction]   = Seq(Diagonal, Left, Up)

  /** The minimum score allowed to start an alignment.  This prevents underflow. */
  val MinStartScore: Int = Int.MinValue / 2

  object AlignmentMatrix {
    def apply(direction: Direction, queryLength: Int, targetLength: Int): AlignmentMatrix = {
      AlignmentMatrix(
        direction = direction,
        scoring   = Matrix(queryLength+1, targetLength+1),
        trace     = Matrix(queryLength+1, targetLength+1))
    }
  }

  /** A single alignment matrix for a given `Direction` storing both the scoring and traceback matrices produce by the aligner. */
  case class AlignmentMatrix(direction: Direction, scoring: Matrix[Int], trace: Matrix[Direction]) {
    val queryLength: Int  = scoring.x - 1
    val targetLength: Int = scoring.y - 1
  }

  /** Represents a cell within the set of matrices used for alignment. */
  private case class MatrixLocation(queryIndex: Int, targetIndex: Int, direction: Direction)

  /** A trait that specifies how the aligner will ask for scoring information. */
  trait AlignmentScorer {
    /**
      * Provides a score for a pairwise alignment of bases.
      *
      * @param queryBase the base from the query sequence
      * @param targetBase the base from the target sequence
      * @return an integer score
      */
    def scorePairing(queryBase: Byte, targetBase: Byte): Int

    /**
      * Provides a score for a gap open or extension.
      *
      * The position of the gap-base being considered is provided as a single 0-based offset into each of the
      * query sequence array and target sequence array.  On the non-gapped side the offset represents the base
      * opposite the site of the gap insertion.  On the gapped side the offset represents the last base
      * before the gap opening.  For example, in the following alignment:
      *
      * qoffset: 01234567
      * query:   ACGTGCAT
      * aln:     ||||--||
      * target:  ACGT--AT
      * toffset: 0123  45
      *
      * We'd see the following invocations:
      *
      * scoreGap(query, target, qOffset=4, tOffset=3, inQuery=false, extend=false)
      * scoreGap(query, target, qOffset=5, tOffset=3, inQuery=false, extend=true)
      *
      * Offsets of -1, query.length and target.length may be passed to indicate that the gap is occurring
      * before the start of one of the sequences, or after the end of the query or target sequence respectively.
      *
      * @param query: the query sequence as a byte array
      * @param target: the target sequence as a byte array
      * @param qOffset: the offset within the query sequence of the gap-base being scored
      * @param tOffset: the offset within the target sequence of the gap-base being scored
      * @param inQuery: true if the gap is on the query side, false if it's on the target side
      * @param extend: true if this is a gap extension, false if it's a gap open
      * @return an integer score
      */
    def scoreGap(query: Array[Byte], target: Array[Byte], qOffset: Int, tOffset: Int, inQuery: Boolean, extend: Boolean): Int
  }
}

/**
  * Implementation of an aligner with generic scoring function and affine gap penalty support.
  * Supports multiple alignment [[Mode]]s for global, semi-global and local alignment.
  *
  * A scoring function (`scoringFunction`) is taken to score pair-wise aligned bases. A default
  * implementation is supplied via the companion [[Aligner]] object which uses a fixed match
  * score and mismatch penalty.
  *
  * When generating CIGARs for alignments the aligner uses the [[isMatch()]] method to determine
  * whether to treat an aligned pair of bases as a match or mismatch.  The default implementation
  * of this method treats U bases as T bases to allow DNA/RNA alignments.  It also will identify
  * as matches any pair of bases (including IUPAC ambiguity codes) that share at least one base
  * in common.  This behaviour can be modified by overriding the [[isMatch()]] method.
  *
  * @param scorer the AlignmentScorer to use to score base pairings and gaps
  * @param useEqualsAndX if true use the = and X cigar operators for matches and mismatches,
  *                      else use the M operator for both.
  * @param mode alignment mode to use when generating alignments
  */
class Aligner(val scorer: AlignmentScorer,
              useEqualsAndX: Boolean = true,
              val mode: Mode = Global) {

  private val (matchOp, mismatchOp) = if (useEqualsAndX) (CigarOperator.EQ, CigarOperator.X) else (CigarOperator.M, CigarOperator.M)

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String): Alignment = align(query.getBytes, target.getBytes)

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String, minScore: Int): Seq[Alignment] = align(query.getBytes, target.getBytes, minScore)

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Alignment]] object describing the optimal global alignment of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte]): Alignment = {
    val matrices = buildMatrices(query, target)
    val location = findBest(matrices)
    generateAlignment(query, target, matrices, location)
  }

  /** Aligns two sequences and returns all alignments that have score >= `minScore`.
    *
    * Note that "all alignments" does not include every permutation of an alignment. E.g. in the case
    * where there is an indel in a repetitive region the results will not include all possible placements
    * of the indel.  Generally speaking a single alignment is produced from each cell in the alignment
    * matrix selected as meeting the score threshold.  This in turn means a single alignment per _end_
    * position.  Specifically by mode:
    *
    *   - Global alignment will always return either 0 or 1 alignments of the whole query and target
    *   - Glocal alignment will return an alignment per end position on the target sequence at which
    *     the end of the query is aligned and has a score >= minScore
    *   - Local will an alignment per location in the matrix with score >= minScore (this can produce
    *     many sub-alignments so be careful!).
    *
    * @param query the query sequence
    * @param target the target sequence
    * @param minScore the minimum alignment score of alignments to return
    * @return a sequence of 0 or more alignments
    */
  def align(query: Array[Byte], target: Array[Byte], minScore: Int): Seq[Alignment] = {
    val matrices = buildMatrices(query, target)
    val locations = findByScore(matrices, minScore)
    locations.map(l => generateAlignment(query, target, matrices, l))
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
    *   - Up   == Consumes query base but not target base, i.e. Insertion vs. the target
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
      val Seq(l, u, d) = Seq(Left, Up, Diagonal).map(matrices.apply)
      (l.scoring, l.trace, u.scoring, u.trace, d.scoring, d.trace)
    }

    // Top left corner - allow all to be zero score but then must be careful to initialize Up and Left to have a gap open
    AllDirections.foreach { direction =>
      matrices(direction).scoring(0, 0) = 0
      matrices(direction).trace(0, 0)   = Done
    }

    fillLeftmostColumn(query, target, leftScoreMatrix, leftTraceMatrix, upScoreMatrix, upTraceMatrix, diagScoreMatrix, diagTraceMatrix)
    fillTopRow(query, target, leftScoreMatrix, leftTraceMatrix, upScoreMatrix, upTraceMatrix, diagScoreMatrix, diagTraceMatrix)
    fillInterior(query, target, leftScoreMatrix, leftTraceMatrix, upScoreMatrix, upTraceMatrix, diagScoreMatrix, diagTraceMatrix)

    matrices
  }

  /** Fills in the leftmost column of the matrices. */
  private final def fillLeftmostColumn(query: Array[Byte],
                                       target: Array[Byte],
                                       leftScoreMatrix: Matrix[Int],
                                       leftTraceMatrix: Matrix[Direction],
                                       upScoreMatrix: Matrix[Int],
                                       upTraceMatrix: Matrix[Direction],
                                       diagScoreMatrix: Matrix[Int],
                                       diagTraceMatrix: Matrix[Direction]): Unit = {
    mode match {
      case Global | Glocal =>
        forloop(from=1, until=query.length+1) { i =>
          leftScoreMatrix(i, 0) = MinStartScore
          leftTraceMatrix(i, 0) = Done
          diagScoreMatrix(i, 0) = MinStartScore
          diagTraceMatrix(i, 0) = Done
          upScoreMatrix(i, 0)   = upScoreMatrix(i-1, 0) + scorer.scoreGap(query, target, qOffset=i-1, tOffset= -1, inQuery=true, extend= i != 1)
          upTraceMatrix(i, 0)   = if (i == 1) Diagonal else Up
        }
      case Local =>
        forloop(from=1, until=query.length+1) { i =>
          leftScoreMatrix(i, 0) = MinStartScore
          leftTraceMatrix(i, 0) = Done
          upScoreMatrix(i, 0)   = MinStartScore
          upTraceMatrix(i, 0)   = Done
          diagScoreMatrix(i, 0) = 0
          diagTraceMatrix(i, 0) = Done
        }
    }
  }

  /** Fills in the top row of the matrices. */
  private final def fillTopRow(query: Array[Byte],
                               target: Array[Byte],
                               leftScoreMatrix: Matrix[Int],
                               leftTraceMatrix: Matrix[Direction],
                               upScoreMatrix: Matrix[Int],
                               upTraceMatrix: Matrix[Direction],
                               diagScoreMatrix: Matrix[Int],
                               diagTraceMatrix: Matrix[Direction]): Unit = {

    mode match {
      case Global =>
        forloop(from=1, until=target.length+1) { j =>
          upScoreMatrix(0, j)    = MinStartScore
          upTraceMatrix(0 , j)   = Done
          diagScoreMatrix(0, j)  = MinStartScore
          diagTraceMatrix(0 , j) = Done
          leftScoreMatrix(0, j) = leftScoreMatrix(0, j-1) + scorer.scoreGap(query, target, -1, j-1, inQuery=false, extend= j != 1)
          leftTraceMatrix(0, j) = if (j == 1) Diagonal else Left
        }
      case Glocal | Local =>
        forloop(from=1, until=target.length+1) { j =>
          leftScoreMatrix(0, j) = MinStartScore
          leftTraceMatrix(0, j) = Done
          upScoreMatrix(0, j)   = MinStartScore
          upTraceMatrix(0, j)   = Done
          diagScoreMatrix(0, j) = 0
          diagTraceMatrix(0, j) = Done
        }
    }
  }

  /** Fills the interior of the matrix. */
  private final def fillInterior(query: Array[Byte],
                                 target: Array[Byte],
                                 leftScoreMatrix: Matrix[Int],
                                 leftTraceMatrix: Matrix[Direction],
                                 upScoreMatrix: Matrix[Int],
                                 upTraceMatrix: Matrix[Direction],
                                 diagScoreMatrix: Matrix[Int],
                                 diagTraceMatrix: Matrix[Direction]): Unit = {
    var i: Int = 1
    val queryLength = query.length
    val targetLength = target.length
    while (i <= queryLength) {
      val iMinusOne = i - 1
      val queryBase = query(iMinusOne)
      var j = 1

      while (j <= targetLength) {
        val jMinusOne = j - 1

        { // Diagonal matrix can come from the previous diagonal, up, or left
          val addend = scorer.scorePairing(queryBase, target(jMinusOne))
          val dScore = diagScoreMatrix(iMinusOne, jMinusOne) + addend
          val lScore = leftScoreMatrix(iMinusOne, jMinusOne) + addend
          val uScore = upScoreMatrix(iMinusOne, jMinusOne)   + addend
          val maxScore = math.max(math.max(dScore, lScore), uScore)

          if (mode == Local && maxScore < 0) {
            diagScoreMatrix(i, j) = 0
            diagTraceMatrix(i, j) = Done
          }
          else if (dScore == maxScore) {
            diagScoreMatrix(i, j) = dScore
            diagTraceMatrix(i, j) = Diagonal
          }
          else if (lScore == maxScore) {
            diagScoreMatrix(i, j) = lScore
            diagTraceMatrix(i, j) = Left
          }
          else {
            diagScoreMatrix(i, j) = uScore
            diagTraceMatrix(i, j) = Up
          }
        }


        { // Up matrix can come from diagonal or up
          val dScore = diagScoreMatrix(iMinusOne, j) + scorer.scoreGap(query, target, i-1, j-1, inQuery=true, extend=false)
          val uScore = upScoreMatrix(iMinusOne, j)   + scorer.scoreGap(query, target, i-1, j-1, inQuery=true, extend=true)
          if (dScore >= uScore) {
            upScoreMatrix(i, j) = dScore
            upTraceMatrix(i, j) = Diagonal
          }
          else {
            upScoreMatrix(i, j) = uScore
            upTraceMatrix(i, j) = Up
          }
        }


        { // Left matrix can come from diagonal or left
          val dScore = diagScoreMatrix(i, jMinusOne) + scorer.scoreGap(query, target, i-1, j-1, inQuery=false, extend=false)
          val lScore = leftScoreMatrix(i, jMinusOne) + scorer.scoreGap(query, target, i-1, j-1, inQuery=false, extend=true)
          if (dScore >= lScore) {
            leftScoreMatrix(i, j) = dScore
            leftTraceMatrix(i, j) = Diagonal
          }
          else {
            leftScoreMatrix(i, j) = lScore
            leftTraceMatrix(i, j) = Left
          }
        }

        j += 1
      }

      i += 1
    }
  }

  /**
    * Given the scoring and trace back matrices, construct the [[Alignment]] object for the (or one of the)
    * optimal alignments.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @param matrices the scoring and trace back matrices for the `Left`, `Up`, and `Diagonal` directions.
    * @return an [[Alignment]] object representing the alignment
    */
  protected def generateAlignment(query: Array[Byte],
                                  target: Array[Byte],
                                  matrices: Array[AlignmentMatrix],
                                  location: MatrixLocation): Alignment = {
    var currOperator: CigarOperator = null
    var currLength: Direction = 0
    val elems = new mutable.ArrayBuffer[CigarElem]()
    var MatrixLocation(curI, curJ, curD) = location

    // The score is always the score from the starting cell
    val score = matrices(curD).scoring(curI, curJ)

    // For global we have to reach the origin, for glocal we just have to reach the top row, and for local we can stop
    // anywhere.  Fortunately, we have initialized the appropriate cells to "Done"
    var nextD = matrices(curD).trace(curI, curJ)

    while (nextD != Done) {
      val op    = (curD: @switch) match {
        case Diagonal =>
          curI -= 1
          curJ -= 1
          if (isMatch(query(curI), target(curJ))) this.matchOp else this.mismatchOp
        case Up       =>
          curI -= 1
          CigarOperator.INSERTION
        case Left     =>
          curJ -= 1
          CigarOperator.DELETION
      }
      curD = nextD

      if (op == currOperator) {
        currLength += 1
      }
      else {
        if (currLength > 0) elems += CigarElem(currOperator, currLength)
        currOperator = op
        currLength   = 1
      }

      nextD = matrices(curD).trace(curI, curJ)
      if (nextD == Done) elems += CigarElem(currOperator, currLength)
    }

    Alignment(query=query, target=target, queryStart=curI+1, targetStart=curJ+1, cigar=Cigar(elems.reverse), score=score)
  }

  /**
    * Finds the matrix location of the single best alignment from which to backtrack and generate the alignment.
    * When there are multiple equally best alignments the choice of which one to return is arbitrary.
    *
    * @param matrices the alignment matrices to search
    */
  private def findBest(matrices: Array[AlignmentMatrix]): MatrixLocation = {
    // Find location based on alignment mode:
    //   - start at the end of query and target for Global
    //   - the highest scoring cell in the last row for Glocal
    //   - the highest scoring cell anywhere in the matrix for Local
    val qLen = matrices(0).queryLength
    val tLen = matrices(0).targetLength

    this.mode match {
      case Global =>
        MatrixLocation(qLen, tLen, AllDirections.maxBy { d => matrices(d).scoring(qLen, tLen) })
      case Glocal =>
        var (maxScore, maxJ, maxD) = (MinStartScore, -1, Done)
        forloop(from = 1, until = tLen + 1) { j =>
          val direction = AllDirections.maxBy(d => matrices(d).scoring(qLen, j))
          val score = matrices(direction).scoring(qLen, j)
          if (score > maxScore) {
            maxScore = score
            maxJ = j
            maxD = direction
          }
        }
        MatrixLocation(qLen, maxJ, maxD)
      case Local =>
        var (maxScore, maxI, maxJ, maxD) = (MinStartScore, -1, -1, Done)
        forloop(from = 1, until = qLen + 1) { i =>
          forloop(from = 1, until = tLen + 1) { j =>
            val direction = AllDirections.maxBy(d => matrices(d).scoring(i, j))
            val score = matrices(direction).scoring(i, j)
            if (score > maxScore) {
              maxScore = score
              maxI = i
              maxJ = j
              maxD = direction
            }
          }
        }
        MatrixLocation(maxI, maxJ, maxD)
    }
  }

  /**
    * Finds all cells in the matrix from which a mode-specific backtrace could be generated and
    * that have a score >= than the minimum score.
    *
    * @return a [[Seq]] of [[MatrixLocation]] in no particular order.
    */
  private def findByScore(matrices: Array[AlignmentMatrix], minScore: Int): Seq[MatrixLocation] = {
    val qLen = matrices(0).queryLength
    val tLen = matrices(0).targetLength
    val hits = new ArrayBuffer[MatrixLocation]()

    this.mode match {
      case Global =>
        hits += findBest(matrices)
      case Glocal =>
        forloop(from = 1, until = tLen + 1) { j =>
          val direction = AllDirections.maxBy(d => matrices(d).scoring(qLen, j))
          val score = matrices(direction).scoring(qLen, j)
          if (score >= minScore) hits += MatrixLocation(qLen, j, direction)
        }
      case Local =>
        forloop(from = 1, until = qLen + 1) { i =>
          forloop(from = 1, until = tLen + 1) { j =>
            val direction = AllDirections.maxBy(d => matrices(d).scoring(i, j))
            if (matrices(direction).scoring(i, j) >= minScore) hits += MatrixLocation(i, j, direction)
          }
        }
    }

    hits
  }

  /** Returns true if the two bases should be considered a match when generating the alignment from the matrix
    * and false otherwise.
    */
  protected def isMatch(b1: Byte, b2: Byte): Boolean = Sequences.compatible(b1, b2)
}
