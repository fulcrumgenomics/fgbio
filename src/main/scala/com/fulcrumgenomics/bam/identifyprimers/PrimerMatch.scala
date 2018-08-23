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
  val NoPrimerMatchName: String = "NoPrimerMatch"

  def toName[T <: PrimerMatch : TypeTag]: String =  typeOf[T] match {
    case t if t =:= typeOf[LocationBasedPrimerMatch]     => "Location"
    case t if t =:= typeOf[UngappedAlignmentPrimerMatch] => "Ungapped"
    case t if t =:= typeOf[GappedAlignmentPrimerMatch]   => "Gapped"
    case _ => unreachable(s"Unknown primer match type: ${this.getClass.getSimpleName}.")
  }

  /** Converts the given [[PrimerMatch]] to a [[String]]. */
  def toName[T <: PrimerMatch](primerMatch: Option[T]): String = primerMatch.map(_.toName).getOrElse(PrimerMatch.NoPrimerMatchName)
}

// TODO: document
private[identifyprimers] sealed trait PrimerMatch {
  def primer: Primer

  final def info(rec: SamRecord): String = {
    val baseInfo = Seq(
      primer.pair_id,
      primer.primer_id,
      primer.ref_name + ":" + primer.start + "-" + primer.end,
      if (primer.positiveStrand) "+" else "-",
      0, // TODO: offset from the 5' end,
      this.toName,
    ).map(_.toString)
    (baseInfo ++ this._info(rec)).mkString(PrimerMatch.InfoDelimiter)
  }

  protected def _info(rec: SamRecord): Seq[Any]

  def toName: String
}

// TODO: document
private[identifyprimers] case class LocationBasedPrimerMatch(primer: Primer, numMismatches: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = Seq(numMismatches)
  def toName: String = PrimerMatch.toName[LocationBasedPrimerMatch]
}

// TODO: document
private[identifyprimers] case class UngappedAlignmentPrimerMatch(primer: Primer, numMismatches: Int, nextNumMismatches: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = {
    val nextOrNa = if (nextNumMismatches == Int.MaxValue) "na" else nextNumMismatches
    Seq(numMismatches, nextOrNa)
  }
  def toName: String = PrimerMatch.toName[UngappedAlignmentPrimerMatch]
}

// TODO: document
private[identifyprimers] case class GappedAlignmentPrimerMatch(primer: Primer, score: Int, secondBestScore: Int) extends PrimerMatch {
  protected def _info(rec: SamRecord): Seq[Any] = Seq(score, secondBestScore)
  def toName: String = PrimerMatch.toName[GappedAlignmentPrimerMatch]
}
