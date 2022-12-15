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
import com.fulcrumgenomics.personal.yfarjoun.NormalizeCoverage._
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.{IntervalListSource, Io}
import htsjdk.samtools.util.{Interval, Locatable, OverlapDetector}
import htsjdk.tribble.SimpleFeature

import java.util
import scala.util.Using
import scala.jdk.CollectionConverters._
import scala.collection.mutable

object NormalizeCoverage {

  def getIntervalsFromDictionary(dict: SequenceDictionary): util.List[Interval] = {
    dict.iterator.map(seq => new Interval(seq.name, 1, seq.length)).toList.asJava
  }

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

  def subsetReadToLocus(rs: Seq[Locatable], overlappingIntervals: Seq[Locatable]): Seq[Locatable] = {
    val unionOverlap: Seq[Locatable] = squish(overlappingIntervals)
    rs.flatMap(r => {
      unionOverlap.filter(l => l.overlaps(r))
        .map(l => intersect(l, r).getOrElse(Set.empty[Locatable]))
    })
  }

  def union(i1: Locatable, i2: Locatable): Set[Locatable] = {
    if (i1.withinDistanceOf(i2, 1))
      Set(new SimpleFeature(i1.getContig, Math.min(i1.getStart, i2.getStart), Math.max(i1.getEnd, i2.getEnd)))
    else
      Set(i1, i2)
  }

  // calculates the locus covered by two loci (assuming they are on the same reference)
  def union(i1: Locatable, i2: Set[Locatable]): Set[Locatable] = {
    i2.flatMap(i => union(i1, i))
  }

  def intersect(i1: Locatable, i2: Locatable): Option[Locatable] = {
    if (!i1.contigsMatch(i2))
      None
    else
      Some(new SimpleFeature(
        i1.getContig,
        Math.max(i1.getStart, i2.getStart),
        Math.min(i1.getEnd, i2.getEnd)))
  }

//  //result will not have overlapping locatables, but will cover the same original territory
//  def squish(ls: IterableOnce[SamRecord]): Seq[Locatable] = {
//    squish(ls.iterator.map(x=>x.asSam))
//  }
  //result will not have overlapping locatables, but will cover the same original territory
  def squish(ls: IterableOnce[Locatable]): Seq[Locatable] = {
    ls.iterator.foldLeft(Seq.empty[Locatable])(union)
  }

  def min(s1: Short, s2: Short): Short = {
    s1 min s2
  }


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
                         @arg(name = "mq", doc = "MAPQ lower bound for considering read-template.") val minMQ: PhredScore = DefaultMinMapQ,
                         @arg(flag = 'c', doc = "Target coverage.") val coverageTarget: Short,
                         @arg(flag = 's', doc = "Order in which to emit the records.") val sortOrder: SamOrder = SamOrder.Coordinate,
                         @arg(flag = 'm', doc = "Max records in RAM.") val maxRecordsInRam: Int = SamWriter.DefaultMaxRecordsInRam

                       ) extends FgBioTool with LazyLogging {


  val overlapDetector: OverlapDetector[Interval] = new OverlapDetector[Interval](0, 0)

  //keeps track of the coverage already added to output
  val coverage: mutable.Map[String, mutable.IndexedSeq[Short]] = mutable.Map()

  //sets up the overlap detector and the coverage tracker
  def initialize(bam: PathToBam, intervalsMaybe: Option[PathToIntervals]): Unit = {
    Using.Manager { use =>
      val in = use(SamSource(bam))
      val intervals: util.List[Interval] = intervalsMaybe.map(f => use(IntervalListSource(f)))
        .map(i => i.toIntervalList.getIntervals)
        .getOrElse(getIntervalsFromDictionary(in.dict))
      overlapDetector.addAll(intervals, intervals)

      in.dict.foreach(seq => {
        coverage(seq.name) = mutable.IndexedSeq.fill(seq.length)(0.toShort)
      }
      )
    }
  }

  def readToSubsettedLocus(rs: Seq[Locatable]): Seq[Locatable] = {
    val squished = squish(rs)
    val overlapping: Seq[Locatable] = rs.flatMap(r => overlapDetector.getOverlaps(r).asScala)
    subsetReadToLocus(squished, overlapping)
  }

  def templateToCoverageLocatable(template: Template): Seq[Locatable] = {
    val locus: Seq[Locatable] = template

      //take all reads
      .allReads

      //make sure they are good
      .filter(readFilterForCoverage)

      //replace each read with its start/end tuple
      .map[Locatable](r => r.asSam).toSeq

    //Subset the read to the overlapping intervals
    readToSubsettedLocus(locus)
  }


  def primariesOverlapTargets(template: Template): Boolean = {
    template.allReads.filter(readFilterForCoverage).exists(r => overlapDetector.overlapsAny(r.asSam))
  }

  //returns Short.MaxValue if template doesn't overlap with requested intervals
  def minCoverage(template: Template): Short = {
    //FIX
    template.allReads.filter(readFilterForCoverage).map(read => {
      val cov: IndexedSeq[Short] = coverage(read.refName).toIndexedSeq
      val covLocusSet: Seq[Locatable] = templateToCoverageLocatable(template)
      val covPerTemplate = covLocusSet.flatMap(covLocus => cov.slice(covLocus.getStart - 1, covLocus.getEnd).minOption)
      covPerTemplate.reduce(min)
    }).reduceOption(min)
      .getOrElse(Short.MaxValue)
  }

  // increments the coverage over the range of the template.
  // will not increment beyond Short.MaxValue
  def incrementCoverage(template: Template): Unit = {
    squish(template.allReads.filter(readFilterForCoverage).map(x=>x.asSam)).foreach(read=>{
    val cov: mutable.IndexedSeq[Short] = coverage(read.getContig)
    val covLocusSet: Seq[Locatable] = templateToCoverageLocatable(template)

    covLocusSet.foreach(covLocus =>
      Range(covLocus.getStart - 1, covLocus.getEnd)
        .foreach(j => cov.update(j, Math.min(Short.MaxValue, cov(j) + 1).toShort)))
  })}

  override def execute(): Unit = {
    checkArguments()

    initialize(input, targetIntervals)

    Using.Manager { use =>
      val in = use(SamSource(input))
      val out = use(SamWriter(output, in.header.clone(), sort = Some(sortOrder), maxRecordsInRam = maxRecordsInRam))

      val templateIterator: Iterator[Template] = Bams.templateSortedIterator(in, maxRecordsInRam)

      templatesToReads(templateIterator).foreach(out.write)
    }
  }

  private def checkArguments(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    targetIntervals.foreach(Io.assertReadable)
    assert(coverageTarget > 0, s"Target must be >0 found $coverageTarget.")
  }

  def templatesToReads(templateIterator: Iterator[Template]): Iterator[SamRecord] = {
    templateIterator
      //see if primary reads overlap region of interest
      .filter(t => primariesOverlapTargets(t))
      //check that all the reads in the template pass minMQ
      .filter(t => minMapQ(t) >= minMQ)
      //only consider reads that cover regions that need coverage, i.e. that
      //at least one base that they cover is less than the target
      .filter(t => minCoverage(t) < coverageTarget)
      //
      .flatMap(t => {
        //increment coverages and emit template
        incrementCoverage(t)
        //output
        t.allReads.iterator
      })
  }
}
