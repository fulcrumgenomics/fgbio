package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.personal.yfarjoun.CoverageManager.{convertSamToInterval, readFilterForCoverage}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}

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
    * will return an empty Interval on contig "" if read is unmapped
    *
    * @param samRecord SamRecord to convert to an Interval
    * @return Interval covering the same position as the alignment of the read
    */
  def convertSamToInterval(samRecord: SamRecord): Interval = {
    if (samRecord.mapped) {
      new Interval(samRecord.refName, samRecord.start, samRecord.end)
    } else {
      new Interval("", 1, 0)
    }
  }

  /**
    * Takes an iterator of Intervals and reduces it to a *non-overlapping* Seq of the same.
    *
    * Result will cover the same original territory, and are guarranteed to not have overlapping Intervals.
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
class CoverageManager(val intervals: IntervalList, val minMapQ:PhredScore, val coverageTarget:Int) {

  private val coverage: OverlapDetector[LocusTrack] = new OverlapDetector[LocusTrack](0, 0)

  IntervalList.getUniqueIntervals(intervals, false)
    .foreach(i => coverage.addLhs(new LocusTrack(i), i))



  private def getCoverageTracks(locus: Interval): Iterator[LocusTrack] = {
    coverage.getOverlaps(locus).map(lt => lt.sliceToLocus(locus))
  }

  /**
    * provides the minimum value in the coverage map in the region that is both
    * a region of interest, and covered by one of coverage-worthy reads in the template.
    *
    * returns Short.MaxValue if template doesn't overlap with requested intervals
    */
  private[yfarjoun] def getMinTemplateCoverage(template: Template): Short = {
    template.allReads.filter(readFilterForCoverage(_, minMapQ)).map(read => {
      getMinCoverage(convertSamToInterval(read))
    }).minOption.getOrElse(Short.MaxValue)
  }

  /**
    * returns the minimum value in the region of interest and in the locus provided
    *
    * returns Short.MaxValue if template doesn't overlap with requested intervals (so that it will
    * register as "not needed" to hit coverage no matter what the target is.
    */
  def getMinCoverage(locus: Interval): Short = {
    getCoverageTracks(locus)
      .map(lt => lt.track.min)
      .minOption.
      getOrElse(Short.MaxValue)
  }


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

  /** Mostly for testing */
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
      .filter(t => getMinTemplateCoverage(t) < coverageTarget)
      //increment coverages and emit reads from template
      .map(t => {
        incrementCoverage(t)
        t
      }).toSeq
  }
}

