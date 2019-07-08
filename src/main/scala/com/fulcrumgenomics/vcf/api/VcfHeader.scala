/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef._
import enumeratum.EnumEntry
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}

import scala.collection.immutable

sealed trait VcfCount { }

object VcfCount {
  case object OnePerAltAllele  extends VcfCount
  case object OnePerAllele     extends VcfCount
  case object OnePerGenotype   extends VcfCount
  case object Unknown          extends VcfCount
  case class Fixed(count: Int) extends VcfCount
}


sealed trait VcfFieldType extends EnumEntry {
  /** Must parse a String into a single value of the type given. */
  def parse(s: String): Any
}

object VcfFieldType extends FgBioEnum[VcfFieldType] {
  case object Integer   extends VcfFieldType { override def parse(s: String): Any = s.toInt }
  case object Float     extends VcfFieldType { override def parse(s: String): Any = s.toFloat }
  case object String    extends VcfFieldType { override def parse(s: String): Any = s }
  case object Character extends VcfFieldType { override def parse(s: String): Any = s.charAt(0) }
  case object Flag      extends VcfFieldType { override def parse(s: String): Any = "." }

  override def values: immutable.IndexedSeq[VcfFieldType] = findValues
}


sealed trait VcfHeaderEntry {}

case class VcfContigHeader(index: Int,
                           name: String,
                           length: Option[Int] = None,
                           assembly: Option[String] = None
                          ) extends VcfHeaderEntry

case class VcfInfoHeader(id: String,
                         count: VcfCount,
                         kind: VcfFieldType,
                         description: String,
                         source: Option[String] = None,
                         version: Option[String] = None
                        ) extends VcfHeaderEntry {}

case class VcfFormatHeader(id: String,
                           count: VcfCount,
                           kind: VcfFieldType,
                           description: String) extends VcfHeaderEntry {}

case class VcfFilterHeader(id: String, description: String) extends VcfHeaderEntry {}

/** Catch all for header line types we don't care enough to have specific implementations of. */
case class VcfGeneralHeader(headerType: String, id: String, data: Map[String, String])

case class VcfHeader(contigs: IndexedSeq[VcfContigHeader],
                     infos: Seq[VcfInfoHeader],
                     formats: Seq[VcfFormatHeader],
                     filters: Seq[VcfFilterHeader],
                     other: Seq[VcfGeneralHeader],
                     samples: IndexedSeq[String]
                    ) {

  /** The contig lines represented as a SAM sequence dictionary. */
  val dict: SAMSequenceDictionary = {
    val d = new SAMSequenceDictionary()
    contigs.map { c => new SAMSequenceRecord(c.name, c.length.getOrElse(SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH))}
      .foreach(r => d.addSequence(r))
    d
  }

  val info: Map[String, VcfInfoHeader]     = infos.map(i => i.id -> i).toMap
  val format: Map[String, VcfFormatHeader] = formats.map(f => f.id -> f).toMap
  val filter: Map[String, VcfFilterHeader] = filters.map(f => f.id -> f).toMap
}


