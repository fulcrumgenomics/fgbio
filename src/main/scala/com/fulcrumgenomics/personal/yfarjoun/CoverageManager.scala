/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.personal.yfarjoun.CoverageManager.{convertSamToInterval, readFilterForCoverage}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList, Locatable, OverlapDetector}

import scala.jdk.CollectionConverters.{IteratorHasAsJava, IteratorHasAsScala}

object CoverageManager {

  /**
    * a Filter for reads determining if they count towards coverage.
    *
    * This filter examins the reads flags and mapping quality. It will pass (i.e. return true)
    * primary and supplementary (non-secondary) non-duplicate reads that are mapped (with mapping quality >= input minMapQ,)
    * and pass PF.
    *
    * @param record  Input record to filter
    * @param minMapQ minimum mapping quality for passing reads
    * @return true if record is mapped with mapQ>=minMapQ, non-duplicate, pf, and non-secondary.
    */
  def readFilterForCoverage(record: SamRecord, minMapQ: PhredScore): Boolean = {
    !record.secondary &&
      record.mapped &&
      !record.duplicate &&
      record.pf &&
      record.mapq >= minMapQ
  }

  /**
    * Method that takes in a SamRecord and returns an Interval with the same position.
    * requires that the read is mapped
    *
    * @param samRecord SamRecord to convert to an Interval
    * @return Interval covering the same position as the alignment of the read
    */
  def convertSamToInterval(samRecord: SamRecord): Interval = {
    assert(samRecord.mapped)
    new Interval(samRecord.refName, samRecord.start, samRecord.end)
  }

  /**
    * Takes an iterator of Intervals and reduces it to a *non-overlapping* Seq of the same.
    *
    * Result will cover the same original territory, and are guaranteed to not have overlapping Intervals.
    */
  def squish(ls: IterableOnce[Interval]): IterableOnce[Interval] = {
    new IntervalList.IntervalMergerIterator(ls.iterator.asJava, true, false, false)
  }.asScala
}

/**
  * Class for managing coverage over an interval list.
  *
  * Holds an IntervalList of interest, a target coverage value and a minimum mapping quality for reads to be considered
  * as providing coverage.
  *
  */
class CoverageManager(val intervals: IndexedSeq[Interval], val minMapQ: PhredScore, val coverageTarget: Int) {

  private val coverage: OverlapDetector[LocusTrack] = new OverlapDetector[LocusTrack](0, 0)

  intervals.foreach(i => coverage.addLhs(new LocusTrack(i), i))

  private def getCoverageTracks(locus: Locatable): Iterator[LocusTrack] = {
    coverage.getOverlaps(locus).map(lt => lt.sliceToLocus(locus))
  }

  /**
    * Determines if a Template contains a read that provides coverage over a region that needs it.
    */
  private[yfarjoun] def needsCoverage(template: Template): Boolean = {
    template.allReads.filter(readFilterForCoverage(_, minMapQ))
      .exists(r => needsCoverage(r.asSam))
  }

  /**
    * Determines if a locus contains overlaps a region that needs coverage.
    */
  def needsCoverage(locus: Locatable): Boolean = {
    getCoverageTracks(locus).exists(lt => lt.track.exists(_ < coverageTarget))
  }

  /**
    * Increments the coverage map over an interval.
    * Will not increment beyond Short.MaxValue to avoid overflow
    */
  private[yfarjoun] def incrementCoverage(interval: Interval): Unit = {
    getCoverageTracks(interval).foreach(locusTrack =>
      locusTrack.track.indices
        .foreach(j => locusTrack.track(j) = Math.min(Short.MaxValue, locusTrack.track(j) + 1).toShort))
  }

  /**
    * Increments the coverage map over the range of the relevant reads of the template.
    * Will not increment beyond Short.MaxValue to avoid overflow
    */
  private[yfarjoun] def incrementCoverage(template: Template): Unit = {
    CoverageManager
      .squish(template.allReads.filter(readFilterForCoverage(_, minMapQ)).map(convertSamToInterval))
      .iterator.foreach(read => incrementCoverage(read))
  }

  /** for testing */
  private[yfarjoun] def resetCoverage(): Unit = {
    coverage.getAll.foreach(locusTrack =>
      locusTrack.track.indices
        .foreach(j => locusTrack.track.update(j, 0.toShort)))
  }

  /**
    * Filters, an Iterator[Template] while recording the resulting coverage
    * **with side-effects** as follows:
    *
    * 1. Only templates that have reads that count towards coverage, AND are needed in order to hit the target coverage
    * pass the filters.
    * 2. coverage over the overlapping region is incremented (only by 1 if overlapping a read which counts towards coverage)
    * 3. templates are returned in their entirety, or not at all. no partial Templates are returned.
    * 4. Any template which is returned contains at least one record whose coverage was used towards the target.
    * 5. returns a Seq rather than an Iterator so that the side-effects will happen for sure.
    *
    * @param templateIterator input Template Iterator
    * @return filtered and processed templates.
    */
  def processTemplates(templateIterator: Iterator[Template]): Seq[Template] = {
    templateIterator
      // Only consider templates that contain reads which cover regions
      // that need coverage, i.e. that at least one base that they cover is less
      // than the target coverage
      .filter(t => needsCoverage(t))
      //increment coverages and emit reads from template
      .map(t => {
        incrementCoverage(t)
        t
      }).toSeq
  }
}

