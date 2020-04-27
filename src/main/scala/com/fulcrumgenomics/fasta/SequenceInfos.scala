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

import com.fulcrumgenomics.FgBioDef._
import enumeratum.EnumEntry
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}


/** Stores information about a single Sequence (ex. chromosome, contig)
  *
  * @param name the primary name of the sequence
  * @param length the length of the sequence, or zero if unknown
  * @param aliases the list of aliases for this sequence
  * @param attributes attributes of this sequence
  */
case class SequenceInfo(name: String,
                        length: Int,
                        aliases: Seq[String] = Seq.empty,
                        attributes: Map[String, Any] = Map.empty
                       ) {

  allNames.foreach { name => SAMSequenceRecord.validateSequenceName(name) }
  require(length > 0, s"Length must be > 0 for '$name'")

  /** All names, including aliases */
  @inline final def allNames: Seq[String] = name +: aliases
  @inline final def apply[T](key: String): T = this.attributes(key).asInstanceOf[T]
  @inline final def get[T](key: String): Option[T] = this.attributes.get(key).map(_.asInstanceOf[T])
  @inline final def contains(key: String): Boolean = this.attributes.contains(key)

  @inline final def isAlternate: Boolean = this.contains("AH")
  def alternate: Option[AlternateLocus] = {
    this.apply[String]("AH") match {
      case "*" => None
      case string =>
        val Array(refName, rest) = string.split(':')
        val Array(start, end) = rest.split('-')
        val locus = refName match {
          case "=" => AlternateLocus(refName=this.name, start=start.toInt, end=end.toInt)
          case _   => AlternateLocus(refName=refName, start=start.toInt, end=end.toInt)
        }
        Some(locus)
    }
  }

  @inline final def m5: String = this.apply(SAMSequenceRecord.MD5_TAG)
  @inline final def m5Int: BigInt = BigInt(this.m5, 16)

  @inline final def assembly: String = this.apply(SAMSequenceRecord.ASSEMBLY_TAG)
  @inline final def uri: String = this.apply(SAMSequenceRecord.URI_TAG)
  @inline final def species: String = this.apply(SAMSequenceRecord.SPECIES_TAG)
  @inline final def description: String = this.apply(SAMSequenceRecord.DESCRIPTION_TAG)
  @inline final def topology: Topology = Topology.values.find(_.name == this.apply("TP")).get

  /** Returns true if the the sequences share a common reference name (including aliases), have the same length, and
    * the same MD5 if both have MD5s. */
  def same(that: SequenceInfo): Boolean = {
    val thisM5: Option[BigInt] = if (this.contains(SAMSequenceRecord.MD5_TAG)) Some(this.m5Int) else None
    val thatM5: Option[BigInt] = if (that.contains(SAMSequenceRecord.MD5_TAG)) Some(that.m5Int) else None

    this.allNames.exists(that.allNames.contains) &&
      this.length == that.length &&
      thisM5 == thatM5
  }

}

/** Contains an ordered collection of sequences. */
class SequenceInfos(infos: IndexedSeq[SequenceInfo]) extends Iterable[SequenceInfo] {
  // Ensure that all names, even aliases, are unique.
  {
    val duplicateNames: String = infos.flatMap(_.allNames)
      .groupBy(identity)
      .filter(_._2.length > 1)
      .keys
      .mkString(", ")
    require(duplicateNames.isEmpty, f"Found duplicate names: $duplicateNames")
  }

  private val mapping: Map[String, SequenceInfo] = infos.flatMap { info =>
    info.allNames.map { name => name -> info }
  }.toMap

  def apply(name: String): SequenceInfo = this.mapping(name)
  def get(name: String): Option[SequenceInfo] = this.mapping.get(name)
  def contains(name: String): Boolean = this.mapping.contains(name)
  override def iterator: Iterator[SequenceInfo] = this.infos.iterator
}


/** Contains useful converters to and from HTSJDK objects. */
object Converters {

  implicit class ToSequenceRecord(info: SequenceInfo) {
    def asSam(index: Option[Int] = None): SAMSequenceRecord = {
      val rec = new SAMSequenceRecord(info.name, info.length)
      index.foreach(rec.setSequenceIndex)
      info.attributes.foreach { case (key, value) => rec.setAttribute(key, value.toString) }
      rec
    }
  }

  implicit class FromSequenceRecord(rec: SAMSequenceRecord) {
    def fromSam(): SequenceInfo = {
      val attributes: Map[String, Any] = rec.getAttributes.map { entry =>
        entry.getKey -> entry.getValue
      }.toMap
      val aliases = attributes.getOrElse[String]("AN", "").split(',')
      SequenceInfo(
        name       = rec.getSequenceName,
        length     = rec.getSequenceLength,
        aliases    = aliases,
        attributes = attributes
      )
    }
  }

  implicit class ToSequenceDictionary(infos: SequenceInfos) {
    def asSam(): SAMSequenceDictionary = {
      val recs = infos.iterator.zipWithIndex.map { case (info, index) =>
        info.asSam(index=Some(index))
      }.toJavaList
      new SAMSequenceDictionary(recs)
    }
  }

  implicit class FromSequenceDictionary(dict: SAMSequenceDictionary) {
    def fromSam(): SequenceInfos = new SequenceInfos(dict.getSequences.map(_.fromSam()).toIndexedSeq)
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