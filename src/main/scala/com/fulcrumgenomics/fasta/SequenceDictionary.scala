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


/** Stores information about a single Sequence (ex. chromosome, contig)
  *
  * @param name the primary name of the sequence
  * @param length the length of the sequence, or zero if unknown
  * @param attributes attributes of this sequence
  */
case class SequenceMetadata(name: String,
                            length: Int,
                            attributes: Map[String, String] = Map.empty) {
  allNames.foreach { name => SAMSequenceRecord.validateSequenceName(name) }
  require(length >= 0, s"Length must be >= 0 for '$name'")
  require(attributes.keys.forall(_ != SAMSequenceRecord.SEQUENCE_NAME_TAG),
    f"`${SAMSequenceRecord.SEQUENCE_NAME_TAG}` should not given in the list of attributes")
  require(attributes.keys.forall(_ != SAMSequenceRecord.SEQUENCE_LENGTH_TAG),
    s"`${SAMSequenceRecord.SEQUENCE_LENGTH_TAG}` should not given in the list of attributes")

  @inline final def apply(key: String): String = this.attributes(key)
  @inline final def string(key: String): String = this.attributes(key).toString
  @inline final def get(key: String): Option[String] = this.attributes.get(key)
  @inline final def getString(key: String): Option[String] = this.attributes.get(key).map(_.toString)
  @inline final def contains(key: String): Boolean = this.attributes.contains(key)

  @inline final def aliases: Seq[String] = this.get("AN").getOrElse("").split(',')
  /** All names, including aliases */
  @inline final def allNames: Seq[String] = name +: aliases

  @inline final def isAlternate: Boolean = this.alternate.isDefined
  lazy val alternate: Option[AlternateLocus] = {
    this.get("AH").flatMap {
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

  @inline final def md5: Option[String] = this.get(SAMSequenceRecord.MD5_TAG)
  @inline final def md5Int: Option[BigInt] = this.md5.map(BigInt(_, 16))

  @inline final def assembly: Option[String] = this.get(SAMSequenceRecord.ASSEMBLY_TAG)
  @inline final def uri: Option[String] = this.get(SAMSequenceRecord.URI_TAG)
  @inline final def species: Option[String] = this.get(SAMSequenceRecord.SPECIES_TAG)
  //@inline final def description: Option[String] = this.get(SAMSequenceRecord.DESCRIPTION_TAG) // FIXME
  @inline final def topology: Option[Topology] = {
    this.get("TP").flatMap(tp => Topology.values.find(_.name == tp))
  }

  /** Returns true if the the sequences share a common reference name (including aliases), have the same length, and
    * the same MD5 if both have MD5s. */
  def same(that: SequenceMetadata): Boolean = {
    val md5Match = (this.md5Int, that.md5Int) match {
      case (Some(thisMd5Int), Some(thatMd5Int)) => thisMd5Int == thatMd5Int
      case _                                    => true
    }

    this.length == that.length &&
      (this.name == that.name || this.allNames.exists(that.allNames.contains)) &&
      md5Match
  }

}

/** Contains an ordered collection of sequences. */
case class SequenceDictionary(infos: IndexedSeq[SequenceMetadata]) extends Iterable[SequenceMetadata] {

  private val mapping: Map[String, SequenceMetadata] = {
    // validates that the same name is not present twice across allNames (so includes aliases)
    val names: mutable.Set[String] = mutable.HashSet[String]()
    infos.flatMap { info =>
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
}


/** Contains useful converters to and from HTSJDK objects. */
object Converters {

  implicit class ToSequenceRecord(info: SequenceMetadata) {
    def asSam(index: Option[Int] = None): SAMSequenceRecord = {
      val rec = new SAMSequenceRecord(info.name, info.length)
      index.foreach(rec.setSequenceIndex)
      info.attributes.foreach { case (key, value) => rec.setAttribute(key, value.toString) }
      rec
    }
  }

  implicit class FromSequenceRecord(rec: SAMSequenceRecord) {
    def fromSam(): SequenceMetadata = {
      val attributes: Map[String, String] = rec.getAttributes.map { entry =>
        entry.getKey -> entry.getValue
      }.toMap
      SequenceMetadata(
        name       = rec.getSequenceName,
        length     = rec.getSequenceLength,
        attributes = attributes
      )
    }
  }

  implicit class ToSequenceDictionary(infos: SequenceDictionary) {
    def asSam(): SAMSequenceDictionary = {
      val recs = infos.iterator.zipWithIndex.map { case (info, index) =>
        info.asSam(index=Some(index))
      }.toJavaList
      new SAMSequenceDictionary(recs)
    }
  }

  implicit class FromSequenceDictionary(dict: SAMSequenceDictionary) {
    def fromSam(): SequenceDictionary = SequenceDictionary(dict.getSequences.map(_.fromSam()).toIndexedSeq)
  }
}

/** The base trait for all topologies. */
sealed trait Topology extends EnumEntry {
  def name: String = this.name.toLowerCase
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