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
  case object OnePerAltAllele        extends VcfCount
  /** VcfCount that indicates there should be exactly one value per allele, including reference, on the Variant. */
  case object OnePerAllele           extends VcfCount
  /** VcfCount that indicates there should be exactly one value per possible genotype.  E.g. for a variant with a
    * ref allele and one alt allele, and a diploid sample, there would be three values. */
  case object OnePerGenotype         extends VcfCount
  /** VcfCount that indicates there should be one value per allele in the sample's GT, e.g. two for a diploid call
    * (`Number=P`, VCF 4.4). */
  case object OnePerGtAllele         extends VcfCount
  /** As [[OnePerAltAllele]], counting only the alternate alleles named by the sample's LAA (`Number=LA`, VCF 4.5). */
  case object OnePerLocalAltAllele   extends VcfCount
  /** As [[OnePerAllele]], counting only the alternate alleles named by the sample's LAA (`Number=LR`, VCF 4.5). */
  case object OnePerLocalAllele      extends VcfCount
  /** As [[OnePerGenotype]], counting only the alternate alleles named by the sample's LAA (`Number=LG`, VCF 4.5). */
  case object OnePerLocalGenotype    extends VcfCount
  /** VcfCount that indicates there should be one value per possible base modification in the alleles of the sample's
    * GT (`Number=M`, VCF 4.5). */
  case object OnePerBaseModification extends VcfCount
  /** The attribute may have any number of values. */
  case object Unknown                extends VcfCount
  /** The attribute must have a fixed number of values, indicated by `count`. */
  case class Fixed(count: Int)       extends VcfCount
}

/** A version of the VCF specification that a VCF header can declare. */
sealed trait VcfVersion { }

object VcfVersion {
  case object V3_2 extends VcfVersion
  case object V3_3 extends VcfVersion
  case object V4_0 extends VcfVersion
  case object V4_1 extends VcfVersion
  case object V4_2 extends VcfVersion
  case object V4_3 extends VcfVersion
  case object V4_4 extends VcfVersion
  case object V4_5 extends VcfVersion
}

/** Trait for the set of acceptable atomic types that can appear in VCF's INFO and genotype fields. */
sealed trait VcfFieldType extends EnumEntry {
  /** Must parse a String into a single value of the type given. */
  def parse(s: String): Any
}

object VcfFieldType extends FgBioEnum[VcfFieldType] {
  /** Integer field type, represented internally as [[scala.Int]]. */
  case object Integer   extends VcfFieldType { override def parse(s: String): Any = if (s == ".") Variant.MissingInt else s.toInt }
  /** Floating point number field type, represented internally as [[scala.Float]].  NaN and infinite values may be
    * written as `NaN`, `Inf` or `Infinity`, in any case and with an optional sign, as the VCF specification allows. */
  case object Float     extends VcfFieldType {
    override def parse(s: String): Any = {
      if (s == ".") Variant.MissingFloat
      else try s.toFloat catch { case ex: NumberFormatException => // Java parses `NaN` and `Infinity`, but not `nan` or `Inf`
        s.toLowerCase(java.util.Locale.ROOT) match {
          case "nan" | "+nan" | "-nan"                   => java.lang.Float.NaN
          case "inf" | "+inf" | "infinity" | "+infinity" => java.lang.Float.POSITIVE_INFINITY
          case "-inf" | "-infinity"                      => java.lang.Float.NEGATIVE_INFINITY
          case _                                         => throw ex
        }
      }
    }
  }
  /** String field type, represented internally as [[java.lang.String]]. */
  case object String    extends VcfFieldType { override def parse(s: String): Any = s }
  /** Character field type, represented internally as [[scala.Char]]. */
  case object Character extends VcfFieldType { override def parse(s: String): Any = s.charAt(0) }
  /** Character field type, represented internally as [[scala.Char]]. */
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
  * @param md5 the optional MD5 checksum of the contig's sequence
  * @param species the optional species from which the contig is drawn
  * @param url the optional URL of the contig's sequence: the `URL` key of the contig line, and the `uri` of a
  *            [[com.fulcrumgenomics.fasta.SequenceMetadata]]
  */
case class VcfContigHeader(index: Int,
                           name: String,
                           length: Option[Int] = None,
                           assembly: Option[String] = None,
                           md5: Option[String] = None,
                           species: Option[String] = None,
                           url: Option[String] = None
                          ) extends VcfHeaderEntry {
  /** Returns the contig as a sequence in a sequence dictionary, keeping its index, length, assembly, MD5, species and
    * URL.  A contig of unknown length is given
    * [[com.fulcrumgenomics.fasta.SequenceMetadata$.UnknownSequenceLength]]. */
  def toSequenceMetadata: SequenceMetadata = SequenceMetadata(
    name     = name,
    length   = length.getOrElse(SequenceMetadata.UnknownSequenceLength),
    index    = index,
    assembly = assembly,
    md5      = md5,
    species  = species,
    uri      = url
  )
}

object VcfContigHeader {
  /** Builds a contig header line from a sequence in a sequence dictionary, keeping its length, assembly, MD5, species
    * and URL.  A sequence of unknown length is written without a length.
    */
  def fromSequenceMetadata(metadata: SequenceMetadata): VcfContigHeader = VcfContigHeader(
    index    = metadata.index,
    name     = metadata.name,
    length   = if (metadata.length == SequenceMetadata.UnknownSequenceLength) None else Some(metadata.length),
    assembly = metadata.assembly,
    md5      = metadata.md5,
    species  = metadata.species,
    url      = metadata.uri
  )

  /** Builds the contig header lines for the sequences in a sequence dictionary (see [[fromSequenceMetadata]]). */
  def fromSequenceDictionary(dict: SequenceDictionary): IndexedSeq[VcfContigHeader] = {
    dict.iterator.map(fromSequenceMetadata).toIndexedSeq
  }
}


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
  * @param version the VCF version the header declares, if any.  A VCF is written at this version, or at
  *                [[VcfVersion.V4_2]] if it is older, and a version too old to express all of the header's lines is
  *                refused when the VCF is written.  A header that declares no version is written at the oldest
  *                version that can express all of its lines, and at least [[VcfVersion.V4_2]].
  */
case class VcfHeader(contigs: IndexedSeq[VcfContigHeader],
                     infos: Seq[VcfInfoHeader],
                     formats: Seq[VcfFormatHeader],
                     filters: Seq[VcfFilterHeader],
                     others: Seq[VcfGeneralHeader],
                     samples: IndexedSeq[String],
                     version: Option[VcfVersion] = None
                    ) {
  /** A mapping of sample name to the index in samples. */
  private[api] val sampleIndex = samples.zipWithIndex.toMap

  /** The contig lines represented as a SAM sequence dictionary. */
  val dict: SequenceDictionary = SequenceDictionary(contigs.map(_.toSequenceMetadata):_*)

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

  /**
    * Returns a copy of this header with the given INFO, FORMAT and FILTER entries added.  Each entry replaces any entry
    * of the same kind with the same ID, and is placed after this header's other entries of its kind.
    *
    * @param entries the INFO, FORMAT and FILTER entries to add
    * @throws java.lang.IllegalArgumentException if any entry is a contig or other entry, which this method does not add,
    *                                            or if two entries of the same kind have the same ID
    */
  def withEntries(entries: Iterable[VcfHeaderEntry]): VcfHeader = addEntries(entries, replace=true)

  /**
    * Returns a copy of this header with the given INFO, FORMAT and FILTER entries added, except for those whose ID is
    * already used by an entry of the same kind, which are kept as they are.  New entries are placed after this header's
    * other entries of their kind.
    *
    * @param entries the INFO, FORMAT and FILTER entries to add
    * @throws java.lang.IllegalArgumentException if any entry is a contig or other entry, which this method does not add,
    *                                            or if two entries of the same kind have the same ID
    */
  def withEntriesIfAbsent(entries: Iterable[VcfHeaderEntry]): VcfHeader = addEntries(entries, replace=false)

  /** Adds the entries for [[withEntries]] and [[withEntriesIfAbsent]], replacing existing entries if `replace` is true. */
  private def addEntries(entries: Iterable[VcfHeaderEntry], replace: Boolean): VcfHeader = {
    entries.foreach {
      case _: VcfInfoHeader | _: VcfFormatHeader | _: VcfFilterHeader => ()
      case entry => throw new IllegalArgumentException(s"Only INFO, FORMAT and FILTER entries can be added to a header, found: $entry")
    }

    def merge[E](existing: Seq[E], added: Iterable[E], id: E => String): Seq[E] = {
      added.groupBy(id).collectFirst { case (dup, es) if es.size > 1 => dup }.foreach { dup =>
        throw new IllegalArgumentException(s"More than one entry of the same kind with ID $dup: ${added.filter(id(_) == dup).mkString(", ")}")
      }
      if (replace) {
        val addedIds = added.iterator.map(id).toSet
        existing.filterNot(e => addedIds.contains(id(e))) ++ added
      }
      else {
        val existingIds = existing.iterator.map(id).toSet
        existing ++ added.filterNot(e => existingIds.contains(id(e)))
      }
    }

    this.copy(
      infos   = merge[VcfInfoHeader](infos, entries.collect { case e: VcfInfoHeader => e }, _.id),
      formats = merge[VcfFormatHeader](formats, entries.collect { case e: VcfFormatHeader => e }, _.id),
      filters = merge[VcfFilterHeader](filters, entries.collect { case e: VcfFilterHeader => e }, _.id)
    )
  }
}
