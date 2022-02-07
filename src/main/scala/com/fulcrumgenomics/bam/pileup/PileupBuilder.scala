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

package com.fulcrumgenomics.bam.pileup

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.bam.pileup.PileupBuilder.BamAccessPattern.{RandomAccess, Streaming}
import com.fulcrumgenomics.bam.pileup.PileupBuilder.PileupParameters
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import enumeratum.EnumEntry
import htsjdk.samtools.CigarOperator.{INSERTION => Insertion}

import java.io.Closeable
import scala.collection.mutable.ArrayBuffer

/** Companion object for [[PileupBuilder]]. */
object PileupBuilder {

  /** Sensible defaults for various implementations of [[PileupBuilder]]. */
  trait PileupParameters {

    /** The minimum mapping quality a record must have to contribute to a pileup by default. */
    val minMapQ: Int = 20

    /** The minimum base quality that a base must have to contribute to a pileup by default. */
    val minBaseQ: Int = 20

    /** Only allow paired records with both records mapped to contribute to a pileup by default. */
    val mappedPairsOnly: Boolean = true

    /** Allow records flagged as duplicates to contribute to a pileup by default. */
    val includeDuplicates: Boolean = false

    /** Allow records flagged as secondary alignments to contribute to a pileup by default. */
    val includeSecondaryAlignments: Boolean = false

    /** Allow records flagged as supplementary alignments to contribute to a pileup by default. */
    val includeSupplementalAlignments: Boolean = false

    /** Exclude any record of an FR pair where the site requested is outside the insert by default. */
    val includeMapPositionsOutsideFrInsert: Boolean = false
  }

  /** A singleton version of all sensible default pileup parameters. */
  object PileupDefaults extends PileupParameters

  /** A trait that all enumerations of BAM access pattern must extend. */
  sealed trait BamAccessPattern extends EnumEntry

  /** The various types of BAM access patterns. */
  object BamAccessPattern extends FgBioEnum[BamAccessPattern] {
    override def values: IndexedSeq[BamAccessPattern] = findValues

    /** The type of BAM access pattern when queries are completed using random access. */
    case object RandomAccess extends BamAccessPattern

    /** The type of BAM access pattern when queries are completed using full BAM streaming. */
    case object Streaming extends BamAccessPattern
  }

  /** Build a [[PileupBuilder]] from a [[SamSource]]. */
  def apply(
    source: SamSource,
    accessPattern: BamAccessPattern             = RandomAccess,
    minMapQ: Int                                = PileupDefaults.minMapQ,
    minBaseQ: Int                               = PileupDefaults.minBaseQ,
    mappedPairsOnly: Boolean                    = PileupDefaults.mappedPairsOnly,
    includeDuplicates: Boolean                  = PileupDefaults.includeDuplicates,
    includeSecondaryAlignments: Boolean         = PileupDefaults.includeSecondaryAlignments,
    includeSupplementalAlignments: Boolean      = PileupDefaults.includeSupplementalAlignments,
    includeMapPositionsOutsideFrInsert: Boolean = PileupDefaults.includeMapPositionsOutsideFrInsert,
  ): PileupBuilder with Closeable = accessPattern match {
    case RandomAccess => RandomAccessPileupBuilder(
      source                             = source,
      minMapQ                            = minMapQ,
      minBaseQ                           = minBaseQ,
      mappedPairsOnly                    = mappedPairsOnly,
      includeDuplicates                  = includeDuplicates,
      includeSecondaryAlignments         = includeSecondaryAlignments,
      includeSupplementalAlignments      = includeSupplementalAlignments,
      includeMapPositionsOutsideFrInsert = includeMapPositionsOutsideFrInsert,
    )
    case Streaming => StreamingPileupBuilder(
      source                             = source,
      minMapQ                            = minMapQ,
      minBaseQ                           = minBaseQ,
      mappedPairsOnly                    = mappedPairsOnly,
      includeDuplicates                  = includeDuplicates,
      includeSecondaryAlignments         = includeSecondaryAlignments,
      includeSupplementalAlignments      = includeSupplementalAlignments,
      includeMapPositionsOutsideFrInsert = includeMapPositionsOutsideFrInsert,
    )
  }

  /** Returns true if <rec> is in a mapped FR pair but the position <pos> is outside the insert coordinates of <rec>.
    * Returns false if <rec> is in a mapped FR pair and the position <pos> is inside the insert coordinates of <rec> or
    * <rec> is not in a mapped FR pair.
    */
  private def positionIsOutsideFrInsert(rec: SamRecord, refIndex: Int, pos: Int): Boolean = {
    rec.isFrPair && {
      val (start, end) = Bams.insertCoordinates(rec)
      rec.refIndex == refIndex && pos >= start && pos <= end
    }
  }

  /** Returns true if the first non-clipping operator is an insertion. */
  private def startsWithInsertion(rec: SamRecord): Boolean = {
    rec.cigar.find(elem => !elem.operator.isClipping).exists(_.operator == Insertion)
  }
}

/** A trait that all pileup builders must extends. */
trait PileupBuilder extends PileupParameters {

  /** The sequence dictionary associated with the records we will pileup. */
  val dict: SequenceDictionary

  /** Pileup records at this position. */
  def pileup(refName: String, pos: Int): Pileup[PileupEntry]

  /** Quickly check the SAM record to see if all the simple static per-read filters accept the read. */
  @inline private final def quickAcceptRecord(rec: SamRecord): Boolean = {
    var compare = true
    if (compare && minMapQ > 0) compare = rec.mapq >= minMapQ
    if (compare && mappedPairsOnly) compare = rec.paired && rec.mapped && rec.mateMapped
    if (compare && !includeDuplicates) compare = !rec.duplicate
    if (compare && !includeSecondaryAlignments) compare = !rec.secondary
    if (compare && !includeSupplementalAlignments) compare = !rec.supplementary
    compare
  }

  /** Quickly check the pileup entry to see if all the simple static per-base filters accept the entry. */
  @inline private final def quickAcceptEntry(entry: PileupEntry): Boolean = {
    entry match {
      case p: BaseEntry => p.qual >= minBaseQ
      case _            => true
    }
  }

  /** Custom SAM record filters to further filter down pileups made from this builder. */
  protected val recFilters: ArrayBuffer[SamRecord => Boolean] = ArrayBuffer.empty[SamRecord => Boolean]

  /** Custom pileup entry filters to further filter down pileups made from this builder. */
  protected val entryFilters: ArrayBuffer[PileupEntry => Boolean] = ArrayBuffer.empty[PileupEntry => Boolean]

  /** Adds a filter to the set of SAM record filters; filters should return true to retain records and false to discard. */
  def withReadFilter(fn: SamRecord => Boolean): this.type = { recFilters.addOne(fn); this }

  /** Adds a filter to the set of pileup entry filters; filters should return true to retain entries and false to discard. */
  def withEntryFilter(fn: PileupEntry => Boolean): this.type = { entryFilters.addOne(fn); this }

  /** Checks to see if all SAM record filters accept the record. */
  protected def acceptRecord(rec: SamRecord): Boolean = quickAcceptRecord(rec) && recFilters.forall(fn => fn(rec))

  /** Checks to see if all pileup entry filters accept the pileup entry. */
  protected def acceptEntry(entry: PileupEntry): Boolean = quickAcceptEntry(entry) && entryFilters.forall(fn => fn(entry))

  /** Constructs a pileup, at the single requested position, from the records returned by the iterator.
    *
    * @param recs the collection of coordinate-sorted SAM records that may or may not overlap the position
    * @param refName the name of the reference sequence on which the position resides
    * @param pos the 1-based position on the reference sequence at which to construct the pileup
    */
  protected def build(recs: IterableOnce[SamRecord], refName: String, pos: Int): Pileup[PileupEntry] = {
    val refIndex = dict(refName).index.ensuring(_ >= 0, s"Unknown reference sequence name: $refName")
    val pile     = IndexedSeq.newBuilder[PileupEntry]
    if (recs.knownSize >= 0) pile.sizeHint(pile.knownSize)

    @inline def testAndAdd(entry: PileupEntry): Unit = if (this.acceptEntry(entry)) pile.addOne(entry)

    recs.iterator.bufferBetter
      .dropWhile(rec => rec.refIndex < refIndex)
      .takeWhile(rec => rec.refIndex == refIndex && rec.start <= pos + 1)
      .filter { rec =>
        var compare: Boolean = !rec.unmapped
        if (compare) compare = this.acceptRecord(rec)
        if (compare) compare = rec.end >= pos
        if (compare) compare = rec.start <= pos || PileupBuilder.startsWithInsertion(rec)
        if (compare) compare = if (!includeMapPositionsOutsideFrInsert && rec.isFrPair) {
          PileupBuilder.positionIsOutsideFrInsert(rec, refIndex = refIndex, pos = pos)
        } else { compare }
        compare
      }
      .foreach { rec =>
        if (rec.start == pos + 1) { // This site must be an insertion in the read that is at the start of the read.
          testAndAdd(InsertionEntry(rec, 0))
        } else {
          val offset = rec.readPosAtRefPos(pos, returnLastBaseIfDeleted = false)
          if (offset == 0) { // This site must be deleted within the read.
            val deletionPosition = rec.readPosAtRefPos(pos, returnLastBaseIfDeleted = true)
            testAndAdd(DeletionEntry(rec, deletionPosition - 1))
          } else { // This site must be a matched site within the read.
            testAndAdd(BaseEntry(rec, offset - 1))
            // Also check to see if the subsequent base represents an insertion.
            if (offset < rec.length - 1 && rec.refPosAtReadPos(offset + 1) == 0) testAndAdd(InsertionEntry(rec, offset))
          }
        }
      }

    Pileup(refName = refName, refIndex = refIndex, pos = pos, pile = pile.result())
  }
}
