package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.personal.yfarjoun.CoverageManager.{coverageReadsFromTemplate, readFilterForCoverage}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList, OverlapDetector}

import scala.jdk.CollectionConverters.{IteratorHasAsJava, IteratorHasAsScala}

object CoverageManager {

  def minMapQ(template: Template): PhredScore = {
    template.allReads.filter(readFilterForCoverage)
      .map(_.mapq)
      .minOption
      .getOrElse(0)
      .toByte
  }

  def readFilterForCoverage(record: SamRecord): Boolean = {
    !record.secondary && // I'm assuming here that secondary reads to not count towards coverage, but supplementary do
      record.mapped &&
      !record.duplicate &&
      record.pf
  }

  /**
    * Takes an iterator of Intervals and reduces it to a non-overlapping Seq of the same.
    *
    * Result will not have overlapping Intervals, but will cover the same original territory
    */
  def squish(ls: IterableOnce[Interval]): IterableOnce[Interval] = {
    new IntervalList.IntervalMergerIterator(ls.iterator.asJava,true,false,false)
  }.asScala

  /**
    * Gets the minimum value of two input Shorts. Used for left-folding
    */
  def min(s1: Short, s2: Short): Short = {
    s1 min s2
  }

  def coverageReadsFromTemplate(t: Template): Iterator[SamRecord] = {
    t.allReads.filter(readFilterForCoverage)
  }
}

/**
  *
  */
class CoverageManager(val intervals: IntervalList) {

  private val coverage: OverlapDetector[LocusTrack] = new OverlapDetector[LocusTrack](0, 0)

  IntervalList.getUniqueIntervals(intervals, false)
    .foreach(i => coverage.addLhs(new LocusTrack(i), i))

  /**
    * For an input template, returns True if one of its coverage-worthy reads overlap the region of interest
    *
    * @param template input Template
    * @return True only if one of the template's coverage-worthy reads coincide with the region of interest
    */
  def primariesOverlapTargets(template: Template): Boolean = {
    coverageReadsFromTemplate(template).exists(r => overlapsAny(r))
  }

  def overlapsAny(read: SamRecord): Boolean = {
    coverage.overlapsAny(convertSamToInterval(read))
  }

  def overlapsAny(Interval: Interval): Boolean = {
    coverage.overlapsAny(Interval)
  }

  def getCoverageTracks(locus: Interval): Iterator[LocusTrack] = {
    coverage.getOverlaps(locus).map(lt => lt.sliceToLocus(locus))
  }

  /**
    * provides the minimum value in the coverage map in the region that is both
    * a region of interest, and covered by one of coverage-worthy reads in the template.
    *
    * returns Short.MaxValue if template doesn't overlap with requested intervals
    */
  def getMinCoverage(template: Template): Short = {
    template.allReads.filter(readFilterForCoverage).map(read => {
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

  def convertSamToInterval(samRecord: SamRecord): Interval = {
    new Interval(samRecord.refName, samRecord.start, samRecord.end)
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
  def incrementCoverage(template: Template): Unit = {
    CoverageManager
      .squish(template.allReads.filter(readFilterForCoverage).map(convertSamToInterval))
      .iterator.foreach(read => incrementCoverage(read))
  }

  /** mostly for testing, but might come in handy */
  def resetCoverage(): Unit = {
    coverage.getAll.foreach(locusTrack =>
      locusTrack.track.indices
        .foreach(j => locusTrack.track.update(j, 0.toShort)))
  }
}

