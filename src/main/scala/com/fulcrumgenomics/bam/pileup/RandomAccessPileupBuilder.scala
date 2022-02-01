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

import com.fulcrumgenomics.bam.api.QueryType.Overlapping
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.bam.pileup.PileupBuilder._
import com.fulcrumgenomics.fasta.SequenceDictionary

import java.io.Closeable

/** Companion object for [[RandomAccessPileupBuilder]]. */
object RandomAccessPileupBuilder {

  /** Build a random access pileup builder from an indexed SAM source. */
  def apply(
    source: SamSource,
    minMapQ: Int                                = DefaultMinMapQ,
    minBaseQ: Int                               = DefaultMinBaseQ,
    mappedPairsOnly: Boolean                    = DefaultMappedPairsOnly,
    includeDuplicates: Boolean                  = DefaultIncludeDuplicates,
    includeSecondaryAlignments: Boolean         = DefaultIncludeSecondaryAlignments,
    includeSupplementalAlignments: Boolean      = DefaultIncludeSupplementalAlignments,
    includeMapPositionsOutsideFrInsert: Boolean = DefaultIncludeMapPositionsOutsideFrInsert,
  ): RandomAccessPileupBuilder = {
    require(source.indexed, "SAM source must be indexed for random access!")
    new RandomAccessPileupBuilder(
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
}

/** A pileup builder that builds pileups using index-based BAM random access.
  *
  * @param source the indexed [[SamSource]] of records we will pileup.
  * @param minBaseQ the minimum base quality that a base must have to contribute to a pileup.
  * @param minMapQ the minimum mapping quality a record must have to contribute to a pileup.
  * @param mappedPairsOnly if true, only allow paired records with both records mapped to contribute to a pileup.
  * @param includeDuplicates if true, allow records flagged as duplicates to contribute to a pileup.
  * @param includeSecondaryAlignments if true, allow records flagged as secondary alignments to contribute to a pileup.
  * @param includeSupplementalAlignments if true, allow records flagged as supplementary alignments to contribute to a pileup.
  * @param includeMapPositionsOutsideFrInsert if true, include any record of an FR pair where the site requested is outside the insert.
  */
class RandomAccessPileupBuilder private(
  source: SamSource,
  override val minMapQ: Int                                = DefaultMinMapQ,
  override val minBaseQ: Int                               = DefaultMinBaseQ,
  override val mappedPairsOnly: Boolean                    = DefaultMappedPairsOnly,
  override val includeDuplicates: Boolean                  = DefaultIncludeDuplicates,
  override val includeSecondaryAlignments: Boolean         = DefaultIncludeSecondaryAlignments,
  override val includeSupplementalAlignments: Boolean      = DefaultIncludeSupplementalAlignments,
  override val includeMapPositionsOutsideFrInsert: Boolean = DefaultIncludeMapPositionsOutsideFrInsert,
) extends PileupBuilder with Closeable {

  /** The sequence dictionary associated with the records we will pileup. */
  override val dict: SequenceDictionary = source.dict

  /** Pileup records at this position. */
  override def pileup(refName: String, pos: Int): Pileup[PileupEntry] = {
    val end         = Math.min(pos + 1, dict(refName).length) // Add 1bp to "end" in case a read has an insertion there.
    val overlapping = source.query(refName, start = pos, end = end, Overlapping)
    this.build(overlapping, refName = refName, pos = pos)
  }

  /** Close the wrapped SAM source. */
  override def close(): Unit = source.close()
}
