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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.personal.yfarjoun.NormalizeCoverageOptions._
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.{IntervalListSource, Io}
import htsjdk.samtools.util.{Interval, Locatable, OverlapDetector}

import scala.jdk.CollectionConverters._
import scala.collection.mutable

object NormalizeCoverageOptions {
  /** Various default values for the Coverage Normalizer. */
  val DefaultMinMapQ: PhredScore = 30.toByte
}

@clp(description =
  """
Normalizes coverage of an input bam to be at a given coverage for every base in the reference/input interval list.
  """,
  group = ClpGroups.SamOrBam)
class NormalizeCoverage(
                         @arg(flag = 'i', doc = "The input SAM or BAM file.") val input: PathToBam,
                         @arg(flag = 'o', doc = "Output BAM file to write.") val output: PathToBam,
                         @arg(flag = 'l', doc = "The interval list over which to normalize coverage. " +
                           "If provided, templates whose primary reads do not overlap with this region will be dropped entirely.") val targetIntervals: Option[PathToIntervals] = None,
                         @arg(name = "mq"  , doc = "MAPQ lower bound for considering read-template.") val minMQ: PhredScore = DefaultMinMapQ,
                         @arg(flag = 'c', doc = "Target coverage.") val coverageTarget: Int,
                         @arg(flag = 's', doc = "Order in which to emit the records.") val sortOrder: SamOrder = SamOrder.Coordinate,
                         @arg(flag = 'm', doc = "Max records in RAM.") val maxRecordsInRam: Int = SamWriter.DefaultMaxRecordsInRam

                       ) extends FgBioTool with LazyLogging {

  var overlapDetector: OverlapDetector[Interval] = _
  val coverage: mutable.Map[String, mutable.IndexedSeq[Short]] = mutable.Map()

  def initialize(bam: PathToBam, intervalsMaybe: Option[PathToIntervals]): Unit = {
    val in = SamSource(bam)

    //QUESTIOM: do I need to close and 'f'?
    overlapDetector = intervalsMaybe.map(f => IntervalListSource(f))
      .map(i => OverlapDetector.create[Interval](i.toIntervalList.getIntervals))
      .getOrElse(makeOverlapDetectorFromDictionary(in.dict))
    in.dict.foreach(seq => coverage.update(seq.name, mutable.IndexedSeq[Short](seq.length.toShort)))
    in.close()
  }

  def union(i1: Locatable, i2: Locatable): Locatable = {
    new Interval(i1.getContig, Math.min(i1.getStart, i2.getStart), Math.max(i1.getEnd, i2.getEnd))
  }

  def readToSubsettedLocus(r: Locatable): Option[Locatable] = {
    val overlappingIntervals: Set[Interval] = overlapDetector.getOverlaps(r).asScala.toSet
    val unionOverlapMaybe = overlappingIntervals.reduceOption(union)
    unionOverlapMaybe.map(unionOverlap => new Interval(
      r.getContig,
      Math.max(r.getStart, unionOverlap.getStart),
      Math.min(r.getEnd, unionOverlap.getEnd)))
  }

  def templateToCoverageTuple(template: Template): Option[Locatable] = {
    val locus: Locatable = template
      //take primary reads
      .primaryReads
      //make sure they are good
      .filter(readFilterForCoverage)

      //replace each read with its start/end tuple
      .map[Locatable](r => r.asSam)

      //find the union of all the reads
      .reduce(union)

    //Subset the read to the overlapping intervals
    readToSubsettedLocus(locus)
  }

  def primariesOverlapTargets(template: Template): Boolean = {
    template.primaryReads.exists(r => overlapDetector.overlapsAny(r.asSam))
  }

  //returns Short.MaxValue if template doesn't overlap with requested intervals
  def minCoverage(template: Template): Short = {
    val cov: IndexedSeq[Short] = coverage(template.r1.get.refName).toIndexedSeq
    val covLocusMaybe: Option[Locatable] = templateToCoverageTuple(template)

    val retVal:Short=covLocusMaybe.map[Short](covLocus => cov.slice(covLocus.getStart - 1, covLocus.getEnd).minOption[Short].getOrElse(0)).getOrElse(Short.MaxValue)
    retVal
  }

  def incrementCoverage(template: Template): Unit = {
    val cov: mutable.IndexedSeq[Short] = coverage(template.r1.get.refName)
    val covLocusMaybe: Option[Locatable] = templateToCoverageTuple(template)

    covLocusMaybe.foreach(covLocus =>
      Range(covLocus.getStart - 1, covLocus.getEnd)
        .map(i => i.toShort)
        .foreach(j => cov.update(j, (cov(j) + 1).toShort)))
  }

  def readFilterForCoverage(record: SamRecord): Boolean = {
    record.mapped &&
      record.properlyPaired &&
      !record.duplicate &&
      record.pf
  }

  def minMapQ(template: Template): PhredScore = {
    template.primaryReads.filter(readFilterForCoverage)
      .map(_.mapq)
      .minOption
      .getOrElse(0)
      .toByte
  }

  def makeOverlapDetectorFromDictionary(dict: SequenceDictionary): OverlapDetector[Interval] = {
    OverlapDetector.create[Interval](dict.iterator.map(seq => new Interval(seq.name, 1, seq.length)).toList.asJava)
  }

  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    targetIntervals.foreach(Io.assertReadable)
    assert(coverageTarget > 0, s"Target must be >0 found $coverageTarget.")

    initialize(input, targetIntervals)

    val in = SamSource(input)
    val out = SamWriter(output, in.header.clone(), sort = Some(sortOrder), maxRecordsInRam = maxRecordsInRam)

    val templateIterator = Bams.templateSortedIterator(in, maxRecordsInRam)

    templateIterator
      //see if primary reads overlap region of interest
      .filter(t => primariesOverlapTargets(t))
      //check that all the reads in the template pass minMQ
      .filter(t => minMapQ(t) >= minMQ)
      // only include template whose primary reads are all on the same reference
      .filter(t => t.primaryReads.map(r => r.refName).distinct.length == 1)
      //only consider reads that cover regions that need coverage, i.e. that
      //at least one base that they cover is less than the target
      .filter(t => minCoverage(t) < coverageTarget)
      //
      .foreach(t => {
        //increment coverages and emit template
        incrementCoverage(t)
        //write
        t.allReads.foreach(out.write)
      })
  }
}
