/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.bam.api.SamOrder.Coordinate
import com.fulcrumgenomics.bam.api.{SamIterator, SamOrder, SamRecord, SamSource}
import com.fulcrumgenomics.bam.pileup.ByCoordinatePileupBuilder.DefaultInitialCacheSize
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.coord.LocatableUtil.GenomicOrdering
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util.{Interval, Locatable}

import java.io.Closeable
import scala.collection.mutable

/** A class to lazily pileup SAM records for coordinate-maintaining or coordinate-advancing loci.
  *
  * @param records a by-name iterator of SAM records that is assumed to be coordinate-sorted.
  * @param dict the sequence dictionary associated with the reference sequences.
  * @param minBaseQ the minimum base quality that a base must have to contribute to a pileup.
  * @param minMapQ the minimum mapping quality a read must have to contribute to a pileup.
  * @param mappedPairsOnly if true, only allow read-pairs with both reads mapped to contribute to a pileup.
  * @param includeDuplicates if true, allow reads flagged as duplicates to contribute to a pileup.
  * @param includeSecondaryAlignments if true, allow reads flagged as secondary alignments to contribute to a pileup.
  * @param includeSupplementalAlignments if true, allow reads flagged as supplementary alignments to contribute to a pileup.
  * @param excludeMapPositionOutsideFrInsert if true, exclude any record of an FR pair where the site requested is outside the insert.
  * @param initialCacheSize the initial size for the internal SAM record cache, set this to your expected pileup depth.
  */
class ByCoordinatePileupBuilder private(
  records: => SamIterator,
  override val dict: SequenceDictionary,
  override val minMapQ: Int                               = 20,
  override val minBaseQ: Int                              = 20,
  override val mappedPairsOnly: Boolean                   = true,
  override val includeDuplicates: Boolean                 = false,
  override val includeSecondaryAlignments: Boolean        = false,
  override val includeSupplementalAlignments: Boolean     = false,
  override val excludeMapPositionOutsideFrInsert: Boolean = true,
  val initialCacheSize: Int                               = DefaultInitialCacheSize,
) extends PileupBuilder(
  dict                              = dict,
  minMapQ                           = minMapQ,
  minBaseQ                          = minBaseQ,
  mappedPairsOnly                   = mappedPairsOnly,
  includeDuplicates                 = includeDuplicates,
  includeSecondaryAlignments        = includeSecondaryAlignments,
  includeSupplementalAlignments     = includeSupplementalAlignments,
  excludeMapPositionOutsideFrInsert = excludeMapPositionOutsideFrInsert,
) with Closeable {
  import com.fulcrumgenomics.bam.pileup.ByCoordinatePileupBuilder.LocatablePileup

  /** Whether this builder is able to pileup more records from the input iterator of SAM records or not. */
  private var closed: Boolean = false

  /** The last pileup we built over the input SAM records and is cached to save time if the locus is re-requested. */
  private var lastPileup: Option[Pileup[PileupEntry]] = None

  /** Advance this builder to the next requested locus and add all possibly overlapping records to the cache. */
  @inline private def seek(refIndex: Int, pos: Int): Unit = {
    // Drop records up until the next record stands a chance of overlapping the requested locus. Then, take records up
    // until the next record stands a chance of having a start greater than the requested locus plus one. All records in
    // this query have a chance of overlapping the locus so we must then filter down to only those that have an end
    // location greater than or equal to the requested locus. Finally, these records will be added to the cache so that
    // they can be evaluated for pileup at the requested locus. We add one to <pos> for the case when a record starts
    // with an insertion.
    val maybeOverlapping = underlying
      .dropWhile(rec => rec.refIndex < refIndex || (rec.refIndex == refIndex && rec.end < pos))
      .takeWhile(rec => rec.refIndex == refIndex && rec.start <= pos + 1)
      .filter(rec => rec.end >= pos)
    cache.addAll(maybeOverlapping)
  }

  /** A genomic ordering for any locatable that utilizes the sequence dictionary corresponding to the input records. */
  private lazy val byCoordinate: Ordering[Locatable] = GenomicOrdering(dict)

  /** Records that we've accumulated that could overlap another coordinate-advancing call to <advanceTo>. */
  private lazy val cache: mutable.ArrayBuffer[SamRecord] = new mutable.ArrayBuffer[SamRecord](initialCacheSize)

  /** The underlying buffered stream of input SAM records which is lazily summoned. */
  private lazy val underlying: Iterator[SamRecord] = records.filter(_.mapped).bufferBetter

  /** Efficiently advance to the next coordinate-maintaining or coordinate-advancing locus and build a pileup there. */
  def advanceTo(refName: String, pos: Int): Pileup[PileupEntry] = {
    require(!closed, "Cannot advance to a new locus if the pileup builder was closed!")
    val currentLocus = new Interval(refName, pos, pos)
    val refIndex     = dict(refName).index.ensuring(_ >= 0, s"Reference name not in sequence dictionary: $refName")

    // If there is a last pileup defined and it was a locus prior to the requested locus then purge the cache of any
    // records that will not overlap the requested locus, advance the iterator to accumulate maybe-overlapping records
    // back into the cache, and finally pileup records from the requested locus using the cache. If there is a last
    // pileup defined and the requested locus is at the same locus as the last pileup, then return the cached pileup.
    // If there is a last pileup defined and the requested locus is prior to the last pileup then an exception will be
    // raised since it means we are not coordinate-maintaining or coordinate-advancing. When there is no last pileup
    // defined, advance the iterator to accumulate maybe-overlapping records into a fresh cache and then pileup records
    // from the requested locus using the cache.
    val pileup = lastPileup match {
      case Some(last) if byCoordinate.lt(last, currentLocus) =>
        lastPileup = None // Set to `None` now so we can drop the object reference ASAP for garbage collection.
        cache.filterInPlace(rec => (rec.refIndex == refIndex && rec.end >= pos) || rec.refIndex > refIndex)
        this.seek(refIndex = refIndex, pos = pos)
        this.build(cache, refName = refName, pos = pos)
      case Some(last) if byCoordinate.equiv(last, currentLocus) => last
      case Some(last) =>
        throw new IllegalArgumentException(
          s"Queried locus $refName:$pos is behind the previous pileup and should be " +
          s"greater than or equal to: ${last.refName}:${last.pos}"
        )
      case None =>
        this.seek(refIndex = refIndex, pos = pos)
        this.build(cache, refName = refName, pos = pos)
    }

    lastPileup = Some(pileup)
    pileup
  }

  /** Close this pileup builder and clear internal state so that all object references can be garbage collected.q */
  override def close(): Unit = { closed = true; records.safelyClose(); lastPileup = None; cache.clear() }
}

/** Companion object for [[ByCoordinatePileupBuilder]]. */
object ByCoordinatePileupBuilder {

  /** The default initial cache size for pre-allocating an array for a pileup of reads. Set to 300x coverage by default. */
  val DefaultInitialCacheSize: Int = 300

  /** Helper class to ensure pileups are locatable and can be used in coordinate comparison and ordering. */
  private implicit class LocatablePileup(pileup: Pileup[PileupEntry]) extends Locatable {
    override def getContig: String = pileup.refName
    override def getStart: Int     = pileup.pos
    override def getEnd: Int       = getStart
  }

  /** Build a [[ByCoordinatePileupBuilder]] from a coordinate-sorted [[SamSource]]. */
  def apply(
    source: SamSource,
    minMapQ: Int                               = 20,
    minBaseQ: Int                              = 20,
    mappedPairsOnly: Boolean                   = true,
    includeDuplicates: Boolean                 = false,
    includeSecondaryAlignments: Boolean        = false,
    includeSupplementalAlignments: Boolean     = false,
    excludeMapPositionOutsideFrInsert: Boolean = true,
    initialCacheSize: Int                      = DefaultInitialCacheSize,
  ): ByCoordinatePileupBuilder = {
    require(SamOrder(source.header).contains(Coordinate), "SAM source must be coordinate sorted!")
    new ByCoordinatePileupBuilder(
      records                           = source.iterator,
      dict                              = source.dict,
      minMapQ                           = minMapQ,
      minBaseQ                          = minBaseQ,
      mappedPairsOnly                   = mappedPairsOnly,
      includeDuplicates                 = includeDuplicates,
      includeSecondaryAlignments        = includeSecondaryAlignments,
      includeSupplementalAlignments     = includeSupplementalAlignments,
      excludeMapPositionOutsideFrInsert = excludeMapPositionOutsideFrInsert,
      initialCacheSize                  = initialCacheSize
    )
  }
}
