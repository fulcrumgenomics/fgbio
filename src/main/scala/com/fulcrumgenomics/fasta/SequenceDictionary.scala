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

import com.fulcrumgenomics.FgBioDef
import com.fulcrumgenomics.FgBioDef._
import enumeratum.EnumEntry
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}
import scala.collection.mutable
import SequenceMetadata.Keys

object SequenceMetadata {
  /** Keys for standard attributes in [[SequenceMetadata]] */
  object Keys {
    val Aliases              : String = "AN"
    val AlternateLocus       : String = "AH"
    val Assembly             : String = SAMSequenceRecord.ASSEMBLY_TAG
    val Description          : String = SAMSequenceRecord.DESCRIPTION_TAG
    private[fasta] val Length: String = SAMSequenceRecord.SEQUENCE_LENGTH_TAG
    val Md5                  : String = SAMSequenceRecord.MD5_TAG
    private[fasta] val Name  : String = SAMSequenceRecord.SEQUENCE_NAME_TAG
    val Species              : String = SAMSequenceRecord.SPECIES_TAG
    val Topology             : String = "TP"
    val Uri                  : String = SAMSequenceRecord.URI_TAG
  }
}

/** Stores information about a single Sequence (ex. chromosome, contig)
  *
  * @param name the primary name of the sequence
  * @param length the length of the sequence, or zero if unknown
  * @param index the index in the sequence dictionary
  * @param attributes attributes of this sequence
  */
case class SequenceMetadata(name: String,
                            length: Int,
                            index: Int = 0,
                            attributes: Map[String, String] = Map.empty) {
  allNames.foreach { name => SAMSequenceRecord.validateSequenceName(name) }
  require(length >= 0, s"Length must be >= 0 for '$name'")
  require(attributes.keys.forall(_ != Keys.Name),
    f"`${Keys.Name}` should not given in the list of attributes")
  require(attributes.keys.forall(_ != Keys.Length),
    s"`${Keys.Length}` should not given in the list of attributes")

  @inline final def apply(key: String): String = this.attributes(key)
  @inline final def get(key: String): Option[String] = this.attributes.get(key)
  @inline final def contains(key: String): Boolean = this.attributes.contains(key)
  @inline final def aliases: Seq[String] = this.get(Keys.Aliases).map(_.split(',').toSeq).getOrElse(Seq.empty[String])
  /** All names, including aliases */
  @inline final def allNames: Seq[String] = name +: aliases
  @inline final def isAlternate: Boolean = this.alternate.isDefined
  lazy val alternate: Option[AlternateLocus] = {
    this.get(Keys.AlternateLocus).flatMap {
      case "*"   => None
      case range =>
        val (refName: String, start: Int, end: Int) = FgBioDef.parseRange(range)
        val locus: AlternateLocus = {
          if (refName == "=") AlternateLocus(refName=this.name, start=start.toInt, end=end.toInt)
          else AlternateLocus(refName=refName, start=start.toInt, end=end.toInt)
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
    val md5Match = (this.md5Int, that.md5Int) match {
      case (Some(thisMd5Int), Some(thatMd5Int)) => thisMd5Int == thatMd5Int
      case _                                    => true
    }

    this.length == that.length &&
      (this.name == that.name || this.allNames.exists(that.allNames.contains)) &&
      md5Match
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
  override def iterator: Iterator[SequenceMetadata] = this.infos.iterator
  def length: Int = this.infos.length
}


/** Contains useful converters to and from HTSJDK objects. */
object Converters {

  /** Converter from a [[SequenceMetadata]] to a [[SAMSequenceRecord]] */
  implicit class ToSAMSequenceRecord(info: SequenceMetadata) {
    def asSam: SAMSequenceRecord = {
      val rec = new SAMSequenceRecord(info.name, info.length)
      rec.setSequenceIndex(info.index)
      info.attributes.foreach { case (key, value) => rec.setAttribute(key, value.toString) }
      rec
    }
  }

  /** Converter from a [[SAMSequenceRecord]] to a [[SequenceMetadata]] */
  implicit class FromSAMSequenceRecord(rec: SAMSequenceRecord) {
    def fromSam: SequenceMetadata = {
      val attributes: Map[String, String] = rec.getAttributes.map { entry =>
        entry.getKey -> entry.getValue
      }.toMap
      SequenceMetadata(
        name       = rec.getSequenceName,
        length     = rec.getSequenceLength,
        index      = rec.getSequenceIndex,
        attributes = attributes
      )
    }
  }

  /** Converter from a [[SequenceDictionary]] to a [[SAMSequenceDictionary]] */
  implicit class ToSAMSequenceDictionary(infos: SequenceDictionary) {
    def asSam: SAMSequenceDictionary = {
      val recs = infos.iterator.zipWithIndex.map { case (info, index) =>
        info.copy(index=index).asSam
      }.toJavaList
      new SAMSequenceDictionary(recs)
    }
  }

  /** Converter from a [[SAMSequenceDictionary]] to a [[SequenceDictionary]] */
  implicit class FromSAMSequenceDictionary(dict: SAMSequenceDictionary) {
    def fromSam: SequenceDictionary = SequenceDictionary(dict.getSequences.map(_.fromSam).toIndexedSeq)
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

/** Stores information about the coordinates of the alternate locus  */
case class AlternateLocus(refName: String, start: Int, end: Int)