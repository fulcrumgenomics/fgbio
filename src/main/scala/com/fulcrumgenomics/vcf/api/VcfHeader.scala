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
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import enumeratum.EnumEntry

import scala.collection.immutable

/** Trait used to represent how many instances of an attribute are to be expected in a VCFs INFO or genotype fields.*/
sealed trait VcfCount { }

object VcfCount {
  /** VcfCount that indicates there should be exactly one value per alternate allele on the Variant. */
  case object OnePerAltAllele  extends VcfCount
  /** VcfCount that indicates there should be exactly one value per allele, including reference, on the Variant. */
  case object OnePerAllele     extends VcfCount
  /** VcfCount that indicates there should be exactly one value per possible genotype.  E.g. for a variant with a
    * ref allele and one alt allele, and a diploid sample, there would be three values. */
  case object OnePerGenotype   extends VcfCount
  /** The attribute may have any number of values. */
  case object Unknown          extends VcfCount
  /** The attribute must have a fixed number of values, indicated by `count`. */
  case class Fixed(count: Int) extends VcfCount
}

/** Trait for the set of acceptable atomic types that can appear in VCF's INFO and genotype fields. */
sealed trait VcfFieldType extends EnumEntry {
  /** Must parse a String into a single value of the type given. */
  def parse(s: String): Any
}

object VcfFieldType extends FgBioEnum[VcfFieldType] {
  /** Integer field type, represented internally as [[Int]]. */
  case object Integer   extends VcfFieldType { override def parse(s: String): Any = if (s == ".") Variant.MissingInt else s.toInt }
  /** Floating point number field type, represented internally as [[Float]]. */
  case object Float     extends VcfFieldType { override def parse(s: String): Any = if (s == ".") Variant.MissingFloat else s.toFloat }
  /** String field type, represented internally as [[String]]. */
  case object String    extends VcfFieldType { override def parse(s: String): Any = s }
  /** Character field type, represented internally as [[Char]]. */
  case object Character extends VcfFieldType { override def parse(s: String): Any = s.charAt(0) }
  /** Character field type, represented internally as [[Char]]. */
  case object Flag      extends VcfFieldType { override def parse(s: String): Any = "." }

  override def values: immutable.IndexedSeq[VcfFieldType] = findValues
}

/** Trait representing an entry/line in a VCF header. */
sealed trait VcfHeaderEntry {}

// TODO add subclasses for ALT and PEDIGREE lines

/** A contig line/entry in the VCF header.
  *
  * @param index the index of the contig in the list of contigs
  * @param name the name of the contig
  * @param length the length of the contig in bases, optional
  * @param assembly the optional name of the assembly from which the contig is drawn
  */
case class VcfContigHeader(index: Int,
                           name: String,
                           length: Option[Int] = None,
                           assembly: Option[String] = None
                          ) extends VcfHeaderEntry


/**
  * An entry in the VCF header describing an INFO field.
  *
  * @param id the ID of the INFO field, i.e. the value that precedes the `=` in the INFO field
  * @param count the expected count of values for the field
  * @param kind the atomic type of values for the field
  * @param description a description of the field, intended for people not for parsing
  * @param source an optional source describing where the field came from
  * @param version an optional version for the INFO field
  */
case class VcfInfoHeader(id: String,
                         count: VcfCount,
                         kind: VcfFieldType,
                         description: String,
                         source: Option[String] = None,
                         version: Option[String] = None
                        ) extends VcfHeaderEntry {}


/**
  * An entry in the VCF header describing a FORMAT field whose values appear in each genotype.
  *
  * @param id the ID of the fields (i.e. the string that appears in the FORMAT field)
  * @param count the expected count of values per sample
  * @param kind the atomic type of each value
  * @param description a human readable description of the field
  */
case class VcfFormatHeader(id: String,
                           count: VcfCount,
                           kind: VcfFieldType,
                           description: String) extends VcfHeaderEntry {}


/**
  * An entry describing a filter used in either the FILTER field or samples' FT field.
  *
  * @param id the name of the filter as it will appear in the FILTER or FT fields
  * @param description a human readable description of the filter
  */
case class VcfFilterHeader(id: String, description: String) extends VcfHeaderEntry {}


/** Catch all for header line types we don't care enough to have specific implementations of. */
case class VcfGeneralHeader(headerType: String, id: String, data: Map[String, String] = Map.empty) extends VcfHeaderEntry


/** The header of a VCF file.
  *
  * @param contigs the ordered list of contig entries in the VCF header (may be empty)
  * @param infos the list of INFO entries in the header
  * @param formats the list of FORMAT entries in the header
  * @param filters the list of FILTER entries in the header
  * @param others the list of all other header entries
  * @param samples the ordered list of samples that appear in the VCF
  */
case class VcfHeader(contigs: IndexedSeq[VcfContigHeader],
                     infos: Seq[VcfInfoHeader],
                     formats: Seq[VcfFormatHeader],
                     filters: Seq[VcfFilterHeader],
                     others: Seq[VcfGeneralHeader],
                     samples: IndexedSeq[String]
                    ) {

  /** A mapping of sample name to the index in samples. */
  private[api] val sampleIndex = samples.zipWithIndex.toMap

  /** The contig lines represented as a SAM sequence dictionary. */
  val dict: SequenceDictionary = {
    val infos = contigs.map { c =>
      SequenceMetadata(name = c.name, length = c.length.getOrElse(0))
    }
    SequenceDictionary(infos:_*)
  }

  /** The INFO entries organized as a map where the key is the ID of the INFO field. */
  val info: Map[String, VcfInfoHeader]     = infos.map(i => i.id -> i).toMap

  /** The FORMAT entries organized as a map where the key is the ID of the FORMAT field. */
  val format: Map[String, VcfFormatHeader] = formats.map(f => f.id -> f).toMap

  /** The FILTER entries organized as a map where the key is the ID of the FILTER. */
  val filter: Map[String, VcfFilterHeader] = filters.map(f => f.id -> f).toMap

  /** Returns an iterator over all entries in the header. */
  def entries: Iterator[VcfHeaderEntry] = contigs.iterator ++ infos.iterator ++ formats.iterator ++ filters.iterator ++ others.iterator

  /** How many lines total are in the header. */
  def size: Int = entries.size
}
