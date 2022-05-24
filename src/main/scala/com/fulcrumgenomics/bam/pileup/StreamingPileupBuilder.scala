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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource}
import com.fulcrumgenomics.bam.pileup.PileupBuilder._
import com.fulcrumgenomics.bam.pileup.StreamingPileupBuilder.DefaultInitialCacheSize
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.coord.LocatableOrdering
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util.{Interval, Locatable}

import java.io.Closeable
import scala.collection.mutable
import scala.language.reflectiveCalls

/** Companion object for [[StreamingPileupBuilder]]. */
object StreamingPileupBuilder {

  /** The default initial cache size for pre-allocating an array for a pileup of reads. Set to 300x coverage by default. */
  val DefaultInitialCacheSize: Int = 300

  /** Build a streaming pileup builder from a coordinate sorted SAM source. */
  def apply(
    source: SamSource,
    minMapQ: Int                                = PileupDefaults.minMapQ,
    minBaseQ: Int                               = PileupDefaults.minBaseQ,
    mappedPairsOnly: Boolean                    = PileupDefaults.mappedPairsOnly,
    includeDuplicates: Boolean                  = PileupDefaults.includeDuplicates,
    includeSecondaryAlignments: Boolean         = PileupDefaults.includeSecondaryAlignments,
    includeSupplementalAlignments: Boolean      = PileupDefaults.includeSupplementalAlignments,
    includeMapPositionsOutsideFrInsert: Boolean = PileupDefaults.includeMapPositionsOutsideFrInsert,
    initialCacheSize: Int                       = DefaultInitialCacheSize,
  ): StreamingPileupBuilder = {
    require(SamOrder(source.header).contains(Coordinate), "SAM source must be coordinate sorted!")
    lazy val iterator = source.iterator
    new StreamingPileupBuilder(
      records                            = iterator,
      dict                               = source.dict,
      minMapQ                            = minMapQ,
      minBaseQ                           = minBaseQ,
      mappedPairsOnly                    = mappedPairsOnly,
      includeDuplicates                  = includeDuplicates,
      includeSecondaryAlignments         = includeSecondaryAlignments,
      includeSupplementalAlignments      = includeSupplementalAlignments,
      includeMapPositionsOutsideFrInsert = includeMapPositionsOutsideFrInsert,
      initialCacheSize                   = initialCacheSize,
      source                             = Some(iterator)
    )
  }

  /** Helper class to ensure pileups are locatable and can be used in coordinate comparison and ordering. */
  private implicit class LocatablePileup(pileup: Pileup[PileupEntry]) extends Locatable {
    override def getContig: String = pileup.refName
    override def getStart: Int     = pileup.pos
    override def getEnd: Int       = pileup.pos
  }
}

/** A pileup builder that builds coordinate-maintaining or advancing pileups using lazy end-to-end BAM streaming.
  *
  * For every call to <pileup>, this class will advance <records> until SAM records are found that overlap the requested
  * locus. These records will be loaded into an internal cache and then a pileup will be built from them. The pileup
  * will be cached in the event a coordinate-maintaining call to <pileup> is made in the future. If another
  * coordinate-advancing call is made to <pileup> then any previous SAM records are evicted from the internal cache that
  * do not overlap the requested locus and <records> will be advanced to find all additional SAM records that overlap
  * the requested locus.
  *
  * @param records a by-name iterator of SAM records that is assumed to be coordinate-sorted.
  * @param dict the sequence dictionary associated with the records we will pileup.
  * @param minBaseQ the minimum base quality that a base must have to contribute to a pileup.
  * @param minMapQ the minimum mapping quality a record must have to contribute to a pileup.
  * @param mappedPairsOnly if true, only allow paired records with both records mapped to contribute to a pileup.
  * @param includeDuplicates if true, allow records flagged as duplicates to contribute to a pileup.
  * @param includeSecondaryAlignments if true, allow records flagged as secondary alignments to contribute to a pileup.
  * @param includeSupplementalAlignments if true, allow records flagged as supplementary alignments to contribute to a pileup.
  * @param includeMapPositionsOutsideFrInsert if true, include any record of an FR pair where the site requested is outside the insert.
  * @param initialCacheSize the initial size for the internal SAM record cache, set this to your expected pileup depth.
  * @param source an optional source to additionally close upon closing this pileup builder.
  */
class StreamingPileupBuilder private(
  records: => Iterator[SamRecord],
  override val dict: SequenceDictionary,
  override val minMapQ: Int                                = PileupDefaults.minMapQ,
  override val minBaseQ: Int                               = PileupDefaults.minBaseQ,
  override val mappedPairsOnly: Boolean                    = PileupDefaults.mappedPairsOnly,
  override val includeDuplicates: Boolean                  = PileupDefaults.includeDuplicates,
  override val includeSecondaryAlignments: Boolean         = PileupDefaults.includeSecondaryAlignments,
  override val includeSupplementalAlignments: Boolean      = PileupDefaults.includeSupplementalAlignments,
  override val includeMapPositionsOutsideFrInsert: Boolean = PileupDefaults.includeMapPositionsOutsideFrInsert,
  initialCacheSize: Int                                    = DefaultInitialCacheSize,
  source: => Option[{ def close(): Unit }]                 = None,
) extends PileupBuilder with Closeable {
  import com.fulcrumgenomics.bam.pileup.StreamingPileupBuilder.LocatablePileup

  /** Records that we've accumulated that could overlap another coordinate-advancing call to <advanceTo>. */
  private val cache: mutable.ArrayBuffer[SamRecord] = new mutable.ArrayBuffer[SamRecord](initialCacheSize)

  /** Whether this builder is able to pileup more records from the input iterator of SAM records or not. */
  private var closed: Boolean = false

  /** A genomic ordering for any locatable that utilizes the sequence dictionary corresponding to the input records. */
  private val byCoordinate: Ordering[Locatable] = LocatableOrdering(dict)

  /** The last pileup we built over the input SAM records and is cached to save time if the locus is re-requested. */
  private var lastPileup: Option[Pileup[PileupEntry]] = None

  /** The underlying buffered stream of input SAM records which is lazily summoned. */
  private lazy val underlying: Iterator[SamRecord] = records.filter(_.mapped).bufferBetter

  /** Advance this builder to the next requested locus and add all possibly overlapping records to the cache. */
  @inline private def seekAndFillCache(refIndex: Int, pos: Int): Unit = {
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

  /** Efficiently advance to the next coordinate-maintaining or coordinate-advancing locus and build a pileup there. */
  def pileup(refName: String, pos: Int): Pileup[PileupEntry] = {
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
        if (last.refIndex != refIndex) cache.clear() // If we've moved to a new reference sequence index, simply clear.
        else cache.filterInPlace(rec => (rec.refIndex == refIndex && rec.end >= pos) || rec.refIndex > refIndex)
        this.seekAndFillCache(refIndex = refIndex, pos = pos)
        this.build(cache, refName = refName, pos = pos)
      case Some(last) if byCoordinate.equiv(last, currentLocus) => last
      case Some(last) =>
        throw new IllegalArgumentException(
          s"Queried locus $refName:$pos is behind the previous pileup and should be " +
          s"greater than or equal to: ${last.refName}:${last.pos}"
        )
      case None =>
        this.seekAndFillCache(refIndex = refIndex, pos = pos)
        this.build(cache, refName = refName, pos = pos)
    }

    lastPileup = Some(pileup)
    pileup
  }

  /** Close this pileup builder and clear internal state so that all object references can be garbage collected. */
  override def close(): Unit = { closed = true; source.foreach(_.close()); cache.clear(); lastPileup = None }
}
