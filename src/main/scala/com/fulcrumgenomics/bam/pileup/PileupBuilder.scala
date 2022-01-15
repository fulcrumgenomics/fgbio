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

import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.CigarOperator.{INSERTION => Insertion}

import scala.collection.mutable.ArrayBuffer

/** Class that provides methods to build and filter pileups.
  *
  * @param dict the sequence dictionary associated with the reference sequences.
  * @param minBaseQ the minimum base quality that a base must have to contribute to a pileup.
  * @param minMapQ the minimum mapping quality a read must have to contribute to a pileup.
  * @param mappedPairsOnly if true, only allow read-pairs with both records mapped to contribute to a pileup.
  * @param includeDuplicates if true, allow records flagged as duplicates to contribute to a pileup.
  * @param includeSecondaryAlignments if true, allow records flagged as secondary alignments to contribute to a pileup.
  * @param includeSupplementalAlignments if true, allow records flagged as supplementary alignments to contribute to a pileup.
  * @param excludeMapPositionOutsideFrInsert if true, exclude any record of an FR pair where the site requested is outside the insert.
  */
class PileupBuilder(
  val dict: SequenceDictionary,
  val minMapQ: Int                               = 20,
  val minBaseQ: Int                              = 20,
  val mappedPairsOnly: Boolean                   = true,
  val includeDuplicates: Boolean                 = false,
  val includeSecondaryAlignments: Boolean        = false,
  val includeSupplementalAlignments: Boolean     = false,
  val excludeMapPositionOutsideFrInsert: Boolean = true,
) {

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

  /** Quickly check the pileup entry to see if all the simple static per-base filters accept the read. */
  @inline private final def quickAcceptEntry(entry: PileupEntry): Boolean = {
    entry match {
      case p: BaseEntry => p.qual >= minBaseQ
      case _            => true
    }
  }

  /** Custom SAM record filters to further filter down pileups made from this builder. */
  private val recFilters = ArrayBuffer.empty[SamRecord => Boolean]

  /** Custom pileup entry filters to further filter down pileups made from this builder. */
  private val entryFilters = ArrayBuffer.empty[PileupEntry => Boolean]

  /** Adds a filter to the set of SAM record filters. */
  def withReadFilter(fn: SamRecord => Boolean): PileupBuilder = yieldAndThen(this) { recFilters += fn }

  /** Adds a filter to the set of pileup entry filters. */
  def withEntryFilter(fn: PileupEntry => Boolean): PileupBuilder = yieldAndThen(this) { entryFilters += fn }

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
  def build(recs: IterableOnce[SamRecord], refName: String, pos: Int): Pileup[PileupEntry] = {
    val refIndex = dict(refName).index.ensuring(_ >= 0, s"Unknown reference sequence name: $refName")
    val pile     = IndexedSeq.newBuilder[PileupEntry]
    if (recs.knownSize >= 0) pile.sizeHint(pile.knownSize)

    @inline def testAndAdd(entry: PileupEntry): Unit = if (this.acceptEntry(entry)) pile.addOne(entry)

    recs.iterator.bufferBetter
      .dropWhile(rec => rec.refIndex < refIndex)
      .takeWhile(rec => rec.refIndex == refIndex && rec.start <= pos + 1)
      .filter { rec =>
        lazy val startsWithInsertion = PileupBuilder.startsWithInsertion(rec)
        var compare: Boolean = !rec.unmapped
        if (compare) compare = this.acceptRecord(rec)
        if (compare) compare = rec.end >= pos
        if (compare) compare = rec.start <= pos || startsWithInsertion
        if (compare) compare = if (excludeMapPositionOutsideFrInsert && rec.isFrPair) {
          PileupBuilder.positionInsideFrInsert(rec, refName = refName, pos = pos).contains(true)
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
          } else { // This site must be a matched site within the read. Also check and add read-end insertions.
            testAndAdd(BaseEntry(rec, offset - 1))
            if (offset < rec.length - 1 && rec.refPosAtReadPos(offset + 1) == 0) testAndAdd(InsertionEntry(rec, offset))
          }
        }
      }

    Pileup(refName = refName, refIndex = refIndex, pos = pos, pile = pile.result())
  }
}

/** Companion object for [[PileupBuilder]]. */
object PileupBuilder {

  /** Returns true if the position <refName>:<pos> is inside the insert of <rec>, if <rec> is in a mapped FR pair. */
  private def positionInsideFrInsert(rec: SamRecord, refName: String, pos: Int): Option[Boolean] = {
    Option.when(rec.refName == refName && rec.isFrPair) {
      val (start, end) = Bams.insertCoordinates(rec)
      pos >= start && pos <= end
    }
  }

  /** Returns true if the read is mapped and the first non-clipping operator is an insertion. */
  private def startsWithInsertion(rec: SamRecord): Boolean = {
    rec.mapped &&
      rec.cigar.iterator.bufferBetter.dropWhile(_.operator.isClipping).headOption.exists(_.operator == Insertion)
  }
}
