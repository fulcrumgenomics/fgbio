/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.bam.{Bams, AlignmentInfo, Template}
import com.fulcrumgenomics.umi.ConsensusTags
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.Murmur3
import htsjdk.samtools.{SAMFileHeader, SAMUtils}
import org.apache.commons.math3.genetics.RandomKey


/** Trait for specifying BAM orderings. */
sealed trait SamOrder extends Product {
  type A <: Ordered[A]

  /** Sets the appropriate fields in the SAM header for the order. */
  def applyTo(header: SAMFileHeader): SAMFileHeader = {
    header.setSortOrder(sortOrder)
    header.setGroupOrder(groupOrder)
    header.setAttribute("SS", subSort.map(ss => sortOrder.name() + ":" + ss).orNull)
    header
  }

  /** Returns the order's name. */
  def name: String = productPrefix

  /** Returns the sort order that should be set in the SAM header. */
  def sortOrder : SortOrder

  /** Returns the group order that should be set in the SAM header. None will cause GO to be unset. */
  def groupOrder: GroupOrder

  /** The subsort string to be placed *after* the sort order in the SS tag. */
  def subSort: Option[String]

  /** Function to generate the sort key for sorting records into this ordering. */
  def sortkey: (SamRecord) => A
}

object SamOrder {
  /** The set of all possible values. */
  val values: Seq[SamOrder] = Seq(Coordinate, Queryname, Random, RandomQuery, TemplateCoordinate, Unsorted, Unknown)

  /**
    * Returns the SamOrder for the given name. Throws an exception if the name doesn't match a known ordering.
    * Performs case insensitive matching so that we can use capitalized class names, but have backwards
    * compatibility to SAM spec sort orders that are all lower-case.
    * */
  def apply(name: String): SamOrder = values.find(_.name.equalsIgnoreCase(name)).getOrElse {
    throw new NoSuchElementException("No SamOrder: " + name)
  }

  /** If the header represents a known sort order returns that SortOrder otherwise None. */
  def apply(header: SAMFileHeader): Option[SamOrder] = {
    // TODO: Update this to allow partial matches?
    val so = header.getSortOrder
    val go = Option(header.getGroupOrder).getOrElse(GroupOrder.none)

    // TODO: Update this so that if the SS prefix is not SO, discard SS?
    val ss = Option(header.getAttribute("SS")).map(s => s.substring(s.indexOf(':') + 1))
    values.find(sam => sam.sortOrder == so && sam.groupOrder == go && sam.subSort == ss)
  }

  /** Ordering object for coordinate order per the SAM spec - reads without refIndex/coordinate are sorted to the end. */
  case object Coordinate extends SamOrder {
    override type A = CoordinateKey
    override val sortOrder:  SortOrder      = SortOrder.coordinate
    override val groupOrder: GroupOrder     = GroupOrder.none
    override val subSort:    Option[String] = None
    override val sortkey: SamRecord => A = rec => {
      val ref = rec.refIndex
      val idx = if (ref < 0) Int.MaxValue else ref
      CoordinateKey(idx, rec.start, rec.asSam.getFlags)
    }
  }

  /** Ordered key object for Coordinate order. */
  case class CoordinateKey(refIndex: Int, pos: Int, flags: Int) extends Ordered[CoordinateKey] {
    override def compare(that: CoordinateKey): Int = {
      var retval = this.refIndex - that.refIndex
      if (retval == 0) retval = this.pos - that.pos
      if (retval == 0) retval = this.flags - that.flags
      retval
    }
  }

  /** Ordering object for queryname order per the SAM spec. */
  case object Queryname extends SamOrder {
    override type A = QuerynameKey
    override val sortOrder:  SortOrder      = SortOrder.queryname
    override val groupOrder: GroupOrder     = GroupOrder.none
    override val subSort:    Option[String] = None
    override val sortkey: (SamRecord => A)  = r => QuerynameKey(r.name, r.flags)
  }

  /** Ordered key object for Queryname order. */
  case class QuerynameKey(name: String, flags: Int) extends Ordered[QuerynameKey] {
    override def compare(that: QuerynameKey): Int = {
      var retval = this.name.compareTo(that.name)
      if (retval == 0) retval = this.flags - that.flags
      retval
    }
  }

  /** Ordering object for generating a random order over all reads. */
  case object Random extends SamOrder {
    override type A = RandomKey
    private  val hasher = new Murmur3(42)
    override val sortOrder:  SortOrder      = SortOrder.unsorted
    override val groupOrder: GroupOrder     = GroupOrder.none
    override val subSort:    Option[String] = Some("random")
    override val sortkey: (SamRecord => A)  = rec => RandomKey(hasher.hashUnencodedChars(rec.id + rec.basesString), rec.flags)
  }

  /** Key object used by Random sort. */
  final case class RandomKey(hash: Int, flags: Int) extends Ordered[RandomKey] {
    override def compare(that: RandomKey): Int = {
      var retval = Integer.compare(this.hash, that.hash)
      if (retval == 0) retval = Integer.compare(this.flags, that.flags)
      retval
    }
  }

  /** Ordering object for generating a random order with queryname grouping. */
  case object RandomQuery extends SamOrder {
    val HashSeed: Int = 42
    override type A = RandomQueryKey
    private  val hasher = new Murmur3(HashSeed)
    override val sortOrder:  SortOrder      = SortOrder.unsorted
    override val groupOrder: GroupOrder     = GroupOrder.query
    override val subSort:    Option[String] = Some("random-query")
    override val sortkey: SamRecord => A    = rec => RandomQueryKey(hasher.hashUnencodedChars(rec.name), rec.name, rec.flags)
  }

  /** Key object used by RandomQuery sort. */
  final case class RandomQueryKey(hash: Int, name: String, flags: Int) extends Ordered[RandomQueryKey] {
    override def compare(that: RandomQueryKey): Int = {
      var retval = Integer.compare(this.hash, that.hash)
      if (retval == 0) retval = this.name.compareTo(that.name)
      if (retval == 0) retval = Integer.compare(this.flags, that.flags)
      retval
    }
  }

  /**
    * The sort order used by GroupReadByUmi. Sorts reads by the earlier unclipped 5' coordinate of the read
    * pair, the higher unclipped 5' coordinate of the read pair, library, the molecular identifier (see
    * [[com.fulcrumgenomics.umi.ConsensusTags.MolecularId]]), read name, and if R1 has the lower coordinates of the pair.
    */
  case object TemplateCoordinate extends SamOrder {
    override type A = TemplateCoordinateKey
    override val sortOrder:  SortOrder      = SortOrder.unsorted
    override val groupOrder: GroupOrder     = GroupOrder.query
    override val subSort:    Option[String] = Some("template-coordinate")
    override val sortkey: SamRecord => A = rec => {
      // For non-secondary/non-supplementary alignments, or unpaired or mate unmapped reads, use the info in the record.
      // Otherwise, use the info in the pa/pm tags.
      val primary = if (!rec.secondary && !rec.supplementary) AlignmentInfo(rec) else {
        val readPrimaryTag = rec.get[String](Template.ReadPrimarySamTag).getOrElse {
          throw new IllegalStateException(
            s"Missing '${Template.ReadPrimarySamTag}' tag required for TemplateCoordinate; try using ZipperBams or SetMateInformation"
          )
        }
        AlignmentInfo(readPrimaryTag, rec.header)
      }
      val mate = if ((!rec.secondary && !rec.supplementary) || rec.unpaired || rec.mateUnmapped) AlignmentInfo(rec, mate=true) else {
        val matePrimaryTag = rec.get[String](Template.MatePrimarySamTag).getOrElse {
          throw new IllegalStateException(
            s"Missing '${Template.MatePrimarySamTag}' tag required for TemplateCoordinate; try using ZipperBams or SetMateInformation"
          )
        }
        AlignmentInfo(matePrimaryTag, rec.header)
      }
      val readPos = if (primary.unmapped) Int.MaxValue else if (primary.negativeStrand) primary.unclippedEnd else primary.unclippedStart
      val matePos = if (mate.unmapped) Int.MaxValue else if (mate.negativeStrand) mate.unclippedEnd else mate.unclippedStart
      val lib = Option(rec.readGroup).flatMap(rg => Option(rg.getLibrary)).getOrElse("Unknown")
      val mid = rec.get[String](ConsensusTags.MolecularId).map { m =>
        val index: Int = m.lastIndexOf('/')
        if (index >= 0) m.substring(0, index) else m
      }.getOrElse("")

      if (primary.refIndex < mate.refIndex || (primary.refIndex == mate.refIndex && readPos < matePos) ||
           (primary.refIndex == mate.refIndex && readPos == matePos && primary.positiveStrand)) {
        TemplateCoordinateKey(primary.refIndex, mate.refIndex, readPos, matePos, primary.negativeStrand, mate.negativeStrand, lib, mid, rec.name, false)
      }
      else {
        TemplateCoordinateKey(mate.refIndex, primary.refIndex, matePos, readPos, mate.negativeStrand, primary.negativeStrand, lib, mid, rec.name, true)
      }
    }
  }

  /** Sorting key used by the [[TemplateCoordinate]] sort. */
  case class TemplateCoordinateKey(refIndex1: Int,
                                   refIndex2: Int,
                                   pos1: Int,
                                   pos2: Int,
                                   neg1: Boolean,
                                   neg2: Boolean,
                                   library: String,
                                   mid: String,
                                   name : String,
                                   isUpperOfPair: Boolean) extends Ordered[TemplateCoordinateKey] {
    override def compare(that: TemplateCoordinateKey): Int = {
      var retval = Integer.compare(this.refIndex1, that.refIndex1)
      if (retval == 0) retval = Integer.compare(this.refIndex2, that.refIndex2)
      if (retval == 0) retval = Integer.compare(this.pos1, that.pos1)
      if (retval == 0) retval = Integer.compare(this.pos2, that.pos2)
      if (retval == 0) retval = this.neg1.compare(that.neg1)
      if (retval == 0) retval = this.neg2.compare(that.neg2)
      if (retval == 0) retval = this.mid.length.compareTo(that.mid.length)
      if (retval == 0) retval = this.mid.compare(that.mid)
      if (retval == 0) retval = this.name.compareTo(that.name)
      if (retval == 0) retval = this.library.compareTo(that.library)
      if (retval == 0) retval = this.isUpperOfPair.compareTo(that.isUpperOfPair)
      retval
    }
  }

  /** Ordering for the official "unsorted" ordering. */
  case object Unsorted extends SamOrder {
    override type A = Nothing
    override val sortOrder:  SortOrder      = SortOrder.unsorted
    override val groupOrder: GroupOrder     = GroupOrder.none
    override val subSort:    Option[String] = None
    override def sortkey: SamRecord => A = throw new UnsupportedOperationException("Sorting not supported for Unsorted order.")
  }

  /** Ordering for the official "unknown" ordering. */
  case object Unknown extends SamOrder {
    override type A = Nothing
    override val sortOrder:  SortOrder      = SortOrder.unsorted
    override val groupOrder: GroupOrder     = GroupOrder.none
    override val subSort:    Option[String] = None
    override def sortkey: SamRecord => A = throw new UnsupportedOperationException("Sorting not supported for Unknown order.")
  }
}

