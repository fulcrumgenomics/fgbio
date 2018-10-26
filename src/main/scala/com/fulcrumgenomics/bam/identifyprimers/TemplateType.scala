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
import com.fulcrumgenomics.bam.api.SamRecord
import enumeratum.EnumEntry

import scala.collection.immutable


private[identifyprimers] sealed trait TemplateType extends EnumEntry

/** The read structure of a template.  Distinguishes between fragment and reads pairs, and mapped and unmapped. */
private[identifyprimers] object TemplateType extends FgBioEnum[TemplateType] {
  /** A read pair where both ends are mapped. */
  case object MappedPair extends TemplateType
  /** A read pair where exactly one end is mapped. */
  case object PartiallyMappedPair extends TemplateType
  /** A read pair where neither end is mapped. */
  case object UnmappedPair extends TemplateType
  /** A mapped fragment. */
  case object MappedFragment extends TemplateType
  /** An unmapped fragment. */
  case object UnmappedFragment extends TemplateType

  override def values: immutable.IndexedSeq[TemplateType] = findValues

  /** Builds a [[TemplateType]] from a fragment read or paired reads. */
  def apply(r1: SamRecord, r2: Option[SamRecord]): TemplateType = {
    r2 match {
      case None     => if (r1.mapped) MappedFragment else UnmappedFragment
      case Some(_r2) =>
        (r1.mapped, _r2.mapped) match {
          case (true, true)   => MappedPair
          case (true, false)  => PartiallyMappedPair
          case (false, true)  => PartiallyMappedPair
          case (false, false) => UnmappedPair
        }
    }
  }
}
