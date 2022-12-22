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
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.{Interval, IntervalList}

import scala.jdk.CollectionConverters._

object NormalizeCoverage {
  /** Various default values for the Coverage Normalizer. */
  val DefaultMinMapQ: PhredScore = 30.toByte


  def getIntervalsFromDictionary(dict: SequenceDictionary): List[Interval] = {
    dict.iterator.map(seq => new Interval(seq.name, 1, seq.length)).toList
  }
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


  //sets up the overlap detector and the coverage tracker
  def createCoverageManager(bam: PathToBam, intervalsMaybe: Option[PathToIntervals]): CoverageManager = {
    val intervals: IntervalList = intervalsMaybe match {
      case None =>
        val dict = SequenceDictionary.extract(bam)
        val il: IntervalList = new IntervalList(dict.toSam)
        il.addall(getIntervalsFromDictionary(dict).asJava)
        il
      case Some(intervals) => IntervalList.fromPath(intervals)
    }
    new CoverageManager(intervals)
  }


  def templateMinMapQ(t: Template): PhredScore = {
    CoverageManager.coverageReadsFromTemplate(t).map(r=>r.mapq).min.toByte
  }

  def filterTemplatesToCoverage(templateIterator: Iterator[Template], coverageManager: CoverageManager, minMapQ: PhredScore): Iterator[SamRecord] = {
    templateIterator
      //see if primary reads overlap region of interest
      .filter(t => coverageManager.primariesOverlapTargets(t))
      //check that all the reads in the template pass minMQ
      .filter(t => templateMinMapQ(t) >= minMapQ)
      //only consider reads that cover regions that need coverage, i.e. that
      //at least one base that they cover is less than the target coverage
      .filter(t => coverageManager.getMinCoverage(t) < coverageTarget)
      //increment coverages and emit reads from template
      .flatMap(t => {
        coverageManager.incrementCoverage(t)
        //output
        t.allReads.iterator
      })
  }

  override def execute(): Unit = {
    checkArguments()

    val coverageManager = createCoverageManager(input, targetIntervals)


    val in = SamSource(input)
    val header:SAMFileHeader=in.header.clone()
    header.setGroupOrder(GroupOrder.query)
    header.setSortOrder(SortOrder.unsorted)

    val out = SamWriter(output, header, sort = Some(sortOrder), maxRecordsInRam = maxRecordsInRam)

    val templateIterator: Iterator[Template] = Bams.queryGroupedRandomIterator(in.iterator,
      in.header,
      maxRecordsInRam)

    filterTemplatesToCoverage(templateIterator,coverageManager,minMQ).foreach(out.write)

    out.close()
    in.close()
  }


  private def checkArguments(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    targetIntervals.foreach(Io.assertReadable)
    assert(coverageTarget > 0, s"Target must be >0. Found $coverageTarget.")
  }
}
