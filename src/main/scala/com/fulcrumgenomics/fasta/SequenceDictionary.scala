/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.fasta

import java.io.StringWriter

import com.fulcrumgenomics.FgBioDef
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceMetadata.Keys
import com.fulcrumgenomics.util.Io
import enumeratum.EnumEntry
import htsjdk.samtools.util.BufferedLineReader
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceDictionaryCodec, SAMSequenceRecord, SAMTextHeaderCodec}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object SequenceMetadata {
  /** Keys for standard attributes in [[SequenceMetadata]] */
  object Keys {
    val Aliases               : String = "AN"
    val AlternateLocus        : String = "AH"
    val Assembly              : String = SAMSequenceRecord.ASSEMBLY_TAG
    val Description           : String = SAMSequenceRecord.DESCRIPTION_TAG
    private[fasta] val Length : String = SAMSequenceRecord.SEQUENCE_LENGTH_TAG
    val Md5                   : String = SAMSequenceRecord.MD5_TAG
    private[fasta] val Name   : String = SAMSequenceRecord.SEQUENCE_NAME_TAG
    val Species               : String = SAMSequenceRecord.SPECIES_TAG
    val Topology              : String = "TP"
    val Uri                   : String = SAMSequenceRecord.URI_TAG

    val values: Set[String] = Set(
      Aliases, AlternateLocus, Assembly, Description, Length, Md5, Name, Species, Topology, Uri
    )
  }

  /** Builds a new [[SequenceMetadata]]
    *
    * @param name name of the reference sequence (contig)
    * @param length length of the contig
    * @param index index of the contig in a set of contigs (ex. genome or assembly)
    * @param aliases list of reference sequence name aliases
    * @param alternateLocus genomic range in the genome or assembly for which this sequence is an alternate
    * @param assembly genome assembly identifier
    * @param description description of the reference sequence
    * @param md5 MD5 checksum
    * @param species species name
    * @param topology molecule topology
    * @param uri URI of the sequence
    * @param customAttributes any custom attributes (non-standard SAM attributes)
    */
  def apply(name: String,
            length: Int                            = SequenceMetadata.UnknownSequenceLength,
            index: Int                             = SequenceMetadata.UnknownSequenceIndex,
            aliases: Seq[String]                   = Seq.empty,
            alternateLocus: Option[AlternateLocus] = None,
            assembly: Option[String]               = None,
            description: Option[String]            = None,
            md5: Option[String]                    = None,
            species: Option[String]                = None,
            topology: Option[Topology]             = None,
            uri: Option[String]                    = None,
            customAttributes: Map[String, String]        = Map.empty
           ): SequenceMetadata = {
    Keys.values.find(customAttributes.contains).foreach { key =>
      throw new IllegalArgumentException(s"Attributes contains a standard key: $key")
    }
    // Build the attributes
    val buffer: ListBuffer[(String, String)] = ListBuffer()
    if (aliases.nonEmpty)       buffer += ((Keys.Aliases,        aliases.mkString(",")))
    alternateLocus.foreach(l => buffer += ((Keys.AlternateLocus, l.toString)) )
    assembly.foreach      (a => buffer += ((Keys.Assembly,       a)) )
    description.foreach   (d => buffer += ((Keys.Description,    d)) )
    md5.foreach           (m => buffer += ((Keys.Md5,            m)) )
    species.foreach       (s => buffer += ((Keys.Species,        s)) )
    topology.foreach      (t => buffer += ((Keys.Topology,       t.name)) )
    uri.foreach           (u => buffer += ((Keys.Uri,            u)) )

    SequenceMetadata(
      name       = name,
      length     = length,
      index      = index,
      attributes = buffer.toMap ++ customAttributes
    )
  }

  /** Generates a new SequenceMetadata from a SAMSequenceRecord. */
  def apply(rec: SAMSequenceRecord): SequenceMetadata = {
    val attributes: Map[String, String] = rec.getAttributes.map { entry => entry.getKey -> entry.getValue }.toMap

    SequenceMetadata(
      name       = rec.getSequenceName,
      length     = rec.getSequenceLength,
      index      = rec.getSequenceIndex,
      attributes = attributes
    )
  }

  /** The value for sequences of unknown length */
  val UnknownSequenceLength: Int = SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH


  /** The value for sequences of unknown length */
  val UnknownSequenceIndex: Int = SAMSequenceRecord.UNAVAILABLE_SEQUENCE_INDEX

  /** Stores information about the coordinates of the alternate locus  */
  type AlternateLocus = GenomicRange

  object AlternateLocus {
    @inline def apply(refName: String, start: Int, end: Int): AlternateLocus = {
      GenomicRange(refName=refName, start=start, end=end)
    }
  }
}

/** Stores information about a single Sequence (ex. chromosome, contig)
  *
  * @param name the primary name of the sequence
  * @param length the length of the sequence, or zero if unknown
  * @param index the index in the sequence dictionary
  * @param attributes attributes of this sequence
  */
case class SequenceMetadata private[fasta]
(name: String,
 length: Int,
 index: Int,
 attributes: Map[String, String]) {
  import SequenceMetadata.AlternateLocus

  allNames.foreach { name => SAMSequenceRecord.validateSequenceName(name) }
  require(length >= 0, s"Length must be >= 0 for '$name'")
  require(!attributes.contains(Keys.Name), f"`${Keys.Name}` should not given in the list of attributes")
  require(!attributes.contains(Keys.Length), s"`${Keys.Length}` should not given in the list of attributes")

  @inline final def apply(key: String): String = this.attributes(key)
  @inline final def get(key: String): Option[String] = this.attributes.get(key)
  @inline final def getOrElse(key: String, default: String): String = this.attributes.getOrElse(key, default)
  @inline final def contains(key: String): Boolean = this.attributes.contains(key)
  lazy val aliases: Seq[String] = this.get(Keys.Aliases).map(_.split(',').toSeq).getOrElse(Seq.empty[String])
  /** All names, including aliases */
  @inline final def allNames: Seq[String] = name +: aliases
  @inline final def isAlternate: Boolean = this.alternate.isDefined
  lazy val alternate: Option[AlternateLocus] = {
    this.get(Keys.AlternateLocus).flatMap {
      case "*"   => None
      case range =>
        val locus = FgBioDef.parseRange(range) match {
          case GenomicRange("=", start, end) => AlternateLocus(refName=this.name, start=start.toInt, end=end.toInt)
          case _locus                        => _locus
        }
        Some(locus)
    }
  }
  @inline final def md5: Option[String] = this.get(Keys.Md5)
  @inline final def md5Int: Option[BigInt] = this.md5.map(BigInt(_, 16))
  @inline final def assembly: Option[String] = this.get(Keys.Assembly)
  @inline final def uri: Option[String] = this.get(Keys.Uri)
  @inline final def species: Option[String] = this.get(Keys.Species)
  @inline final def description: Option[String] = this.get(Keys.Description)
  @inline final def topology: Option[Topology] = {
    this.get(Keys.Topology).flatMap(tp => Topology.values.find(_.name == tp))
  }

  /** Returns true if the the sequences share a common reference name (including aliases), have the same length, and
    * the same MD5 if both have MD5s. */
  def sameAs(that: SequenceMetadata): Boolean = {
    if (this.length != that.length) false
    else if (this.name != that.name && !this.allNames.exists(that.allNames.contains)) false
    else {
      (this.md5Int, that.md5Int) match {
        case (Some(thisMd5Int), Some(thatMd5Int)) => thisMd5Int == thatMd5Int
        case _                                    => true
      }
    }
  }

  /** Generates a copy of this record as a SAMSequenceRecord. */
  def toSam: SAMSequenceRecord = {
    val rec = new SAMSequenceRecord(name, length)
    rec.setSequenceIndex(index)
    attributes.foreach { case (key, value) => rec.setAttribute(key, value) }
    rec
  }

  override def toString: FilenameSuffix = {
    import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceRecord
    this.asSam.getSAMString
  }
}

object SequenceDictionary {
  /** Builds a sequence dictionary from one or more sequence metadatas.  This will set the sequence index. */
  def apply(info: SequenceMetadata*): SequenceDictionary = {
    val infos = info.zipWithIndex.map { case (info, index) =>
      if (info.index != index) info.copy(index=index) else info
    }.toIndexedSeq
    SequenceDictionary(infos=infos)
  }

  /** Builds a sequence dictionary from the given path */
  def apply(path: PathToSequenceDictionary): SequenceDictionary = {
    import Converters.FromSAMSequenceDictionary
    val codec  = new SAMTextHeaderCodec()
    val reader = new BufferedLineReader(Io.toInputStream(path))
    val dict   = codec.decode(reader, path.toString).getSequenceDictionary
    reader.close()
    dict.fromSam
  }

  /** Generates a new SequenceDictionary from a SAMSequenceDictionary */
  def apply(dict: SAMSequenceDictionary): SequenceDictionary = {
    apply(dict.getSequences.toSeq.map(SequenceMetadata.apply):_*)
  }

  /** Extracts a [[SequenceDictionary]] from SAM/BAM/CRAM/FASTA/DICT files. */
  def extract(path: FilePath): SequenceDictionary = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    SAMSequenceDictionaryExtractor.extractDictionary(path).fromSam
  }
}

/** Contains an ordered collection of sequences. */
case class SequenceDictionary(infos: IndexedSeq[SequenceMetadata]) extends Iterable[SequenceMetadata] {

  private val mapping: Map[String, SequenceMetadata] = {
    // validates that the same name is not present twice across allNames (so includes aliases)
    val names: mutable.Set[String] = mutable.HashSet[String]()
    infos.iterator.zipWithIndex.flatMap { case (info, index) =>
      require(info.index == index, s"Infos must be given with index set correctly. See ${index}th with name: $names")
      info.allNames.map { name =>
        require(!names.contains(name), f"Found duplicate sequence name: $name")
        names.add(name)
        name -> info
      }
    }.toMap
  }

  def apply(name: String): SequenceMetadata = this.mapping(name)
  def get(name: String): Option[SequenceMetadata] = this.mapping.get(name)
  def contains(name: String): Boolean = this.mapping.contains(name)
  def apply(index: Int): SequenceMetadata = this.infos(index)
  def get(index: Int): Option[SequenceMetadata] = this.infos.lift(index)
  override def iterator: Iterator[SequenceMetadata] = this.infos.iterator
  def length: Int = this.infos.length

  /** Returns true if the the sequences share a common reference name (including aliases), have the same length, and
    * the same MD5 if both have MD5s. */
  def sameAs(that: SequenceDictionary): Boolean = {
    this.length == that.length && this.zip(that).forall { case (thisInfo, thatInfo) => thisInfo.sameAs(thatInfo) }
  }

  /** Returns a copy of this sequence dictionary as a SAMSequenceDictionary */
  def toSam: SAMSequenceDictionary = {
    val recs = infos.iterator.map(_.toSam).toJavaList
    new SAMSequenceDictionary(recs)
  }

  /** Writes the sequence dictionary to the given path */
  def write(path: FilePath): Unit = {
    val writer = Io.toWriter(path)
    this.write(writer=writer)
    writer.close()
  }

  /** Writes the sequence dictionary to the given writer */
  def write(writer: java.io.Writer): Unit = {
    import Converters.ToSAMSequenceRecord
    val codec  = new SAMSequenceDictionaryCodec(writer)
    codec.encodeHeaderLine(false)
    this.foreach { metadata: SequenceMetadata =>
      codec.encodeSequenceRecord(metadata.asSam)
    }
  }

  override def toString(): String = {
    val writer = new StringWriter()
    this.write(writer=writer)
    writer.toString
  }
}


/** Contains useful converters to and from HTSJDK objects. */
object Converters {

  /**
    * Converter from a [[SequenceMetadata]] to a [[SAMSequenceRecord]]
    * @deprecated use [[SequenceMetadata.toSam]] instead.
    */
  @Deprecated
  implicit class ToSAMSequenceRecord(info: SequenceMetadata) {
    def asSam: SAMSequenceRecord = info.toSam
  }

  /** Converter from a [[SAMSequenceRecord]] to a [[SequenceMetadata]] */
  implicit class FromSAMSequenceRecord(rec: SAMSequenceRecord) {
    def fromSam: SequenceMetadata = SequenceMetadata(rec)
    def toScala: SequenceMetadata = SequenceMetadata(rec)
  }

  /** Converter from a [[SequenceDictionary]] to a [[SAMSequenceDictionary]]
    * @deprecated use [[SequenceDictionary.toSam]] instead.
    */
  @Deprecated
  implicit class ToSAMSequenceDictionary(dict: SequenceDictionary) {
    def asSam: SAMSequenceDictionary = dict.toSam
  }

  /** Converter from a [[SAMSequenceDictionary]] to a [[SequenceDictionary]] */
  implicit class FromSAMSequenceDictionary(dict: SAMSequenceDictionary) {
    require(dict != null, "The reference provided does not have a sequence dictionary (.dict)")
    def fromSam: SequenceDictionary = SequenceDictionary(dict)
    def toScala: SequenceDictionary = SequenceDictionary(dict)
  }
}

/** The base trait for all topologies. */
sealed trait Topology extends EnumEntry {
  def name: String = this.entryName.toLowerCase
}

/** An enumeration representing the various reference sequence topologies. */
object Topology extends FgBioEnum[Topology] {
  def values: scala.collection.immutable.IndexedSeq[Topology] = findValues
  /** The sequence is linear. */
  case object Linear extends Topology
  /** The sequence is circular. */
  case object Circular extends Topology
}