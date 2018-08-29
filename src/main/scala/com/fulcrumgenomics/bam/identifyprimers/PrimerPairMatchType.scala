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

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import enumeratum.EnumEntry

import scala.collection.immutable


private[identifyprimers] sealed trait PrimerPairMatchType extends EnumEntry

/** Identifies how one or two [[PrimerMatch]]es relate. [[PrimerMatch]]es with the same [[Primer.pair_id]] are called
  * "canonical".  For paired end reads, all types are possible, whereas for fragment reads, only
  * [[PrimerPairMatchType.Single]] and [[PrimerPairMatchType.NoMatch]] are possible. */
private[identifyprimers] object PrimerPairMatchType extends FgBioEnum[PrimerPairMatchType] {
  /** Two primer matches that are from the same "canonical" pair, with one match to the forward and one to the reverse. */
  case object Canonical extends PrimerPairMatchType
  /** Two primer matches that are from the same "canonical" pair, but both matches are to the same primer in the pair. */
  case object SelfDimer extends PrimerPairMatchType
  /** Two primers that are not from the same "canonical" pair and are called cross dimers due to their short template/product size. */
  case object CrossDimer extends PrimerPairMatchType
  /** Two primers that are not from the same "canonical" pair but are not [[CrossDimer]]s. */
  case object NonCanonical extends PrimerPairMatchType
  /** Only primer match. */
  case object Single extends PrimerPairMatchType
  /** No primer matches. */
  case object NoMatch extends PrimerPairMatchType

  /** Returns all the types of primer matches. */
  override def values: immutable.IndexedSeq[PrimerPairMatchType] = findValues

  /** Creates a [[PrimerPairMatchType]] based on the [[PrimerMatch]]es for the forward and reverse primer respectively. */
  def apply(forwardFivePrimeMatch: Option[PrimerMatch],
            reverseFivePrimeMatch: Option[PrimerMatch]): PrimerPairMatchType = {
    this.apply(forwardFivePrimeMatch, reverseFivePrimeMatch, Int.MinValue, Int.MaxValue)
  }

  /** Creates a [[PrimerPairMatchType]] based on the [[PrimerMatch]]es for the forward and reverse primer respectively. */
  def apply(forwardFivePrimeMatch: Option[PrimerMatch],
            reverseFivePrimeMatch: Option[PrimerMatch],
            minInsertLength: Int,
            maxInsertLength: Int) : PrimerPairMatchType = {
    (forwardFivePrimeMatch, reverseFivePrimeMatch) match {
      case (Some(fwd), Some(rev)) if fwd.primer.pair_id != rev.primer.pair_id               =>
        // not from the same pair; cross-dimer if the template/product size is too small, non-canonical otherwise
        if (fwd.primer.ref_name != rev.primer.ref_name) CrossDimer
        else {
          val insertLength = this.insertLength(fwd, rev)
          if (insertLength < minInsertLength || maxInsertLength < insertLength ) CrossDimer else NonCanonical
        }
      case (Some(fwd), Some(rev)) if fwd.primer == rev.primer                               => SelfDimer    // same primer pair, but same primer!
      case (Some(fwd), Some(rev)) if fwd.primer.positive_strand == rev.primer.positive_strand => NonCanonical // same primer pair, but same strand!
      case (Some(_),   Some(_))                                                             => Canonical    // same primer pair, different primers
      case (Some(_), None) | (None, Some(_))                                                => Single       // only one primer match
      case _                                                                                => NoMatch      // no primer matches
    }
  }

  /** Returns the insert length for a given primer pair. The primer pair must map to the same chromosome. */
  def insertLength(fwd: PrimerMatch, rev: PrimerMatch): Int = {
    // not from the same pair; cross-dimer if the template/product size is too small, non-canonical otherwise
    require(fwd.primer.ref_name == rev.primer.ref_name, "Primers must be mapped to the same reference.")
    // check the product size using the primer mapping
    val minStart = math.min(fwd.primer.start, rev.primer.end)
    val maxEnd   = math.max(fwd.primer.end, rev.primer.end)
    maxEnd - minStart + 1
  }
}