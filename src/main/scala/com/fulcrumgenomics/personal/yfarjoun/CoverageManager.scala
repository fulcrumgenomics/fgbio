package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.personal.yfarjoun.CoverageManager.{convertSamToInterval, readFilterForCoverage}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}

import scala.jdk.CollectionConverters.{IteratorHasAsJava, IteratorHasAsScala}

object CoverageManager {


  def readFilterForCoverage(record: SamRecord, minMapQ: PhredScore): Boolean = {
    !record.secondary && // I'm assuming here that secondary reads to not count towards coverage, but supplementary do
      record.mapped &&
      !record.duplicate &&
      record.pf &&
      record.mapq >= minMapQ
  }

  def convertSamToInterval(samRecord: SamRecord): Interval = {
    new Interval(samRecord.refName, samRecord.start, samRecord.end)
  }

  /**
    * Takes an iterator of Intervals and reduces it to a non-overlapping Seq of the same.
    *
    * Result will not have overlapping Intervals, but will cover the same original territory
    */
  def squish(ls: IterableOnce[Interval]): IterableOnce[Interval] = {
    new IntervalList.IntervalMergerIterator(ls.iterator.asJava, true, false, false)
  }.asScala

}

/**
  *
  */
class CoverageManager(val intervals: IntervalList) {

  private val coverage: OverlapDetector[LocusTrack] = new OverlapDetector[LocusTrack](0, 0)

  IntervalList.getUniqueIntervals(intervals, false)
    .foreach(i => coverage.addLhs(new LocusTrack(i), i))

  def getCoverageTracks(locus: Interval): Iterator[LocusTrack] = {
    coverage.getOverlaps(locus).map(lt => lt.sliceToLocus(locus))
  }

  /**
    * provides the minimum value in the coverage map in the region that is both
    * a region of interest, and covered by one of coverage-worthy reads in the template.
    *
    * returns Short.MaxValue if template doesn't overlap with requested intervals
    */
  def getMinCoverage(template: Template, minMapQ: PhredScore): Short = {
    template.allReads.filter(readFilterForCoverage(_, minMapQ)).map(read => {
      getMinCoverage(convertSamToInterval(read))
    }).minOption.getOrElse(Short.MaxValue)
  }

  /**
    * returns the minimum value in the region of interest and in the locus
    *
    * returns Short.MaxValue if template doesn't overlap with requested intervals
    */
  def getMinCoverage(locus: Interval): Short = {
    getCoverageTracks(locus)
      .map(lt => lt.track.min)
      .minOption.
      getOrElse(Short.MaxValue)
  }

  def incrementCoverage(interval: Interval): Unit = {
    getCoverageTracks(interval).foreach(locusTrack =>
      locusTrack.track.indices
        .foreach(j => locusTrack.track(j) = Math.min(Short.MaxValue, locusTrack.track(j) + 1).toShort))
  }

  /**
    * Increments the coverage map over the range of the relevant reads of the template.
    * Will not increment beyond Short.MaxValue to avoid overflow
    */
  def incrementCoverage(template: Template, minMapQ: PhredScore): Unit = {
    CoverageManager
      .squish(template.allReads.filter(readFilterForCoverage(_, minMapQ)).map(convertSamToInterval))
      .iterator.foreach(read => incrementCoverage(read))
  }

  /** mostly for testing, but might come in handy */
  def resetCoverage(): Unit = {
    coverage.getAll.foreach(locusTrack =>
      locusTrack.track.indices
        .foreach(j => locusTrack.track.update(j, 0.toShort)))
  }

  def processTemplates(templateIterator: Iterator[Template], minMapQ: PhredScore, coverageTarget: Short): Iterator[SamRecord] = {
    templateIterator
      // only consider templates that contain reads which cover regions
      // that need coverage, i.e. that at least one base that they cover is less
      // than the target coverage
      .filter(t => getMinCoverage(t, minMapQ) < coverageTarget)
      //increment coverages and emit reads from template
      .flatMap(t => {
        incrementCoverage(t, minMapQ)
        //output
        t.allReads.iterator
      })
  }

}

