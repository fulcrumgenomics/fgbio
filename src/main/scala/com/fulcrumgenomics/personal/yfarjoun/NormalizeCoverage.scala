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
import com.fulcrumgenomics.personal.yfarjoun.NormalizeCoverageOptions._
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.{IntervalListSource, Io}
import htsjdk.samtools.util.{Interval, OverlapDetector}

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
                         @arg(flag = 'm', doc = "MAPQ lower bound for considering read-template.") val minMQ: PhredScore = DefaultMinMapQ,
                         @arg(flag = 'c', doc = "Target coverage.") val coverageTarget: Int,
                         @arg(flag = 's', doc = "Order in which to emit the records.") val sortOrder: SamOrder = SamOrder.Coordinate,
                         @arg(flag = 'm', doc = "Max records in RAM.") val maxRecordsInRam: Int = SamWriter.DefaultMaxRecordsInRam

                       ) extends FgBioTool with LazyLogging {


  def templateToCoverageTuple(template: Template): (Int, Int) = {
    template
      //take primary reads
      .primaryReads
      //make sure they are good
      .filter(readFilterForCoverage)

      //replace each read with its start/end tuple
      .map(r => (r.start, r.end)) reduce ((_, _) => {
      case ((r1_s: Int, r1_e: Int), (r2_s: Int, r2_e: Int)) =>
        (Seq(r1_s, r2_s).min, Seq(r1_e, r2_e).max)
    })
  }

  def primariesOverlapTargets(template: Template, overlapDetector: Option[OverlapDetector[Interval]]): Boolean = {
    if (overlapDetector.isEmpty) {
      true
    } else {
      template.primaryReads
        .map(r => overlapDetector.get.overlapsAny(r.asSam))
        //if any read overlaps target
        .reduce(_ || _)
    }
  }

  def minCoverage(template: Template, cov: Map[String, mutable.IndexedSeq[Short]]): Short = {
    val coverage: IndexedSeq[Short] = cov.getOrElse(template.r1.get.refName, Nil).toIndexedSeq
    templateToCoverageTuple(template) match {
      case (start: Int, end: Int) =>
        coverage.slice(start - 1, end).minOption.getOrElse(0)
    }
  }

  def incrementCoverage(template: Template, cov: Map[String, mutable.IndexedSeq[Short]]): Unit = {
    val coverage: mutable.IndexedSeq[Short] = cov.getOrElse(template.r1.get.refName, Nil)
    val test: Short = coverage(0)
    coverage.update(1, test)

    templateToCoverageTuple(template) match {
      case (start: Int, end: Int) =>
        Range(start - 1, end)
          .map(i => i.toShort)
          .foreach(j => coverage.update(j) = (coverage(j) + 1).toShort)
    }
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

  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    targetIntervals.foreach(Io.assertReadable)
    assert(coverageTarget > 0, s"Target must be >0 found $coverageTarget.")

    val in = SamSource(input)
    val out = SamWriter(output, in.header.clone(), sort = Some(sortOrder), maxRecordsInRam = maxRecordsInRam)

    val overlapDetector: Option[OverlapDetector[Interval]] = targetIntervals.map(f => IntervalListSource(f))
      .map(i => OverlapDetector.create[Interval](i.toIntervalList.getIntervals))

    val templateIterator = Bams.templateSortedIterator(in, maxRecordsInRam)
    val coverages = in.dict.map(seq => seq.name -> mutable.IndexedSeq[Short](seq.length.toShort)).toMap

    templateIterator
      //see if primary reads overlap region of interest
      .filter(t => primariesOverlapTargets(t, overlapDetector))
      //check that all the reads in the template pass minMQ
      .filter(t => minMapQ(t) >= minMQ)
      // only include template whose primary reads are all on the same reference
      .filter(t => t.primaryReads.map(r => r.refName).distinct.length == 1)
      //only consider reads that cover regions that need coverage, i.e. that
      //at least one base that they cover is less than the target
      .filter(t => minCoverage(t, coverages) < coverageTarget)
      //
      .foreach(t => {
        //increment coverages and emit template
        incrementCoverage(t, coverages)
        //write
        t.allReads.foreach(out.write)
      })
  }
}
