/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.identifyprimers

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef.unreachable

import scala.reflect.runtime.universe._

private[identifyprimers] object PrimerMatch {
  val InfoDelimiter: String     = ","
  val LocationName: String      = "Location"
  val UngappedName: String      = "Ungapped"
  val GappedName: String        = "Gapped"
  val NoPrimerMatchName: String = "NoPrimerMatch"

  /** Returns the canonical name of the primer match. */
  def toName[T <: PrimerMatch : TypeTag]: String = typeOf[T] match {
    case t if t =:= typeOf[LocationBasedPrimerMatch]     => LocationName
    case t if t =:= typeOf[UngappedAlignmentPrimerMatch] => UngappedName
    case t if t =:= typeOf[GappedAlignmentPrimerMatch]   => GappedName
    case _ => unreachable(s"Unknown primer match type: ${this.getClass.getSimpleName}.")
  }

  /** Returns the canonical name of the primer match, or `NoPrimerMatchName` otherwise. */
  def toName[T <: PrimerMatch](primerMatch: Option[T]): String = primerMatch.map(_.toName).getOrElse(PrimerMatch.NoPrimerMatchName)
}

/** Base trait for all types of primer matches */
private[identifyprimers] sealed trait PrimerMatch {
  def primer: Primer

  /** Returns a comma-delimited [[String]] of information about a primer match.
    *
    * The following fields are returned:
    * - `pair_id` - the unique identifier for the primer pair
    * - `primer_id` - the unique identifier for the primer in a primer pair
    * - `<ref_name>:<start>-<end>` - the primer chromosome/contig, start (1-based), and end (1-based, inclusive)
    * - `strand` - '+' if the forward primer, '-' otherwise
    * - `read_num` - the read number of the match (1 or 2).
    * - `match_offset` - the offset in the read where the primer matches
    * - `match_length`- the length of the primer match in the read
    * - `match_type` - how the match was found; valid values are 'location', 'gapped', or 'ungapped'
    * - `match_type_info` - additional information based on the `match_type`, containing one or more values
    *
    * The `match_type_info` for each `match_type` is as follow:
    * - location:
    *   - `num_mismatches` - the number of mismatches between the primer and read sequence
    * - ungapped:
    *   - `num_mismatches` - the number of mismatches between the primer and read sequence
    *   - 'next_best` - the number of mismatches in the next best ungapped alignment (or "na" if none was found)
    * - gapped
    *   - `score` - the alignment score
    *   - 'next_best` - the alignment score of then next best alignment
    */
  final def info(rec: SamRecord): String = {
    val baseInfo = Seq(
      primer.pair_id,
      primer.primer_id,
      primer.ref_name + ":" + primer.start + "-" + primer.end,
      if (primer.positive_strand) "+" else "-",
      if (!rec.paired || rec.firstOfPair) 1 else 2,
      0, // TODO: start of the primer match (0-based)
      0, // TODO: length of primer match
      this.toName
    ).map(_.toString)
    (baseInfo ++ this._info(rec)).mkString(PrimerMatch.InfoDelimiter)
  }

  /** All classes should extend this to add any extra match-specific information. */
  protected def _info(rec: SamRecord): Seq[Any]

  def toName: String
}

private[identifyprimers] case class LocationBasedPrimerMatch(primer: Primer, numMismatches: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = Seq(numMismatches)
  def toName: String = PrimerMatch.toName[LocationBasedPrimerMatch]
}

private[identifyprimers] case class UngappedAlignmentPrimerMatch(primer: Primer, numMismatches: Int, nextNumMismatches: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = {
    val nextOrNa = if (nextNumMismatches == Int.MaxValue) "na" else nextNumMismatches
    Seq(numMismatches, nextOrNa)
  }
  def toName: String = PrimerMatch.toName[UngappedAlignmentPrimerMatch]
}

private[identifyprimers] case class GappedAlignmentPrimerMatch(primer: Primer, score: Int, secondBestScore: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = Seq(score, secondBestScore)
  def toName: String = PrimerMatch.toName[GappedAlignmentPrimerMatch]
}
