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
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.personal.yfarjoun.NormalizeCoverage._
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList, SequenceUtil}

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
    |Normalizes coverage of an input bam to be at a given coverage for every base in the reference/input interval list.
    |
    |Only non-secondary, non-duplicate, mapped (with mapping quality >= mq argument), pf reads are considered. A set of reads
    |from the same queryname will increment the depth by 1 if it covers the requested region. Templates (the complete collection
    |of reads with the same queryname) will be emitted to the output file if they were found to be needed to provide coverage
    |for at least one base in the region of interest (or the entire reference specified in the input bam, if no such interval
    |is provided.)
    |
    |The order at which the templates are considered is random (based on the query name) but deterministic.
    |
    |Reads are assumed to provide coverage from the "alignment start" to the "alignment end" positions. No special
    |consideration is provided for deletions or gaps in the alignment.
""",
  group = ClpGroups.SamOrBam)
class NormalizeCoverage(
                         @arg(flag = 'i', doc = "The input SAM or BAM file.") val input: PathToBam,
                         @arg(flag = 'o', doc = "Output BAM file to write the templates that are needed to bring " +
                           "the coverage to the target value.") val output: PathToBam,
                         @arg(flag = 'l', doc = "The interval list over which to normalize coverage. If not provided program will " +
                           "assume that coverage over entire reference is required.") val targetIntervals: Option[PathToIntervals] = None,
                         @arg(name = "mq", doc = "Mapping Quality lower bound for considering a read.") val minMQ: PhredScore = DefaultMinMapQ,
                         @arg(flag = 'c', doc = "Target coverage. ") val coverageTarget: Short,
                         @arg(flag = 's', doc = "Order in which to emit the records.") val sortOrder: SamOrder = SamOrder.Coordinate,
                         @arg(flag = 'm', doc = "Max records in RAM.") val maxRecordsInRam: Int = SamWriter.DefaultMaxRecordsInRam

                       ) extends FgBioTool with LazyLogging {


  //sets up the overlap detector and the coverage tracker
  private def createCoverageManager(bam: PathToBam, intervalsMaybe: Option[PathToIntervals]): CoverageManager = {
    val dict = SequenceDictionary.extract(bam)
    val intervals = intervalsMaybe match {
      // if no intervals are provided, make one from the embedded dictionary
      case None =>
        val il = new IntervalList(dict.toSam)
        il.addall(getIntervalsFromDictionary(dict).asJava)
        il
      case Some(intervals) =>
        val il = IntervalList.fromPath(intervals)
        // make sure that input il and bam have compatible dictionaries.
        SequenceUtil.assertSequenceDictionariesEqual(dict.toSam, il.getHeader.getSequenceDictionary,true)
        il
    }
    new CoverageManager(intervals, minMapQ = minMQ, coverageTarget = coverageTarget)
  }

  override def execute(): Unit = {
    checkArguments()

    val coverageManager = createCoverageManager(input, targetIntervals)

    val in = SamSource(input)

    // create an output writer with the correct output order
    val out = SamWriter(output,
      in.header.clone(),
      sort = Some(sortOrder),
      maxRecordsInRam = maxRecordsInRam)

    // create a random-ordered template iterator
    val templateIterator = Bams.templateRandomIterator(in.iterator,
      in.header,
      maxRecordsInRam)

    // find the templates required for coverage and add them to the output writer
    coverageManager
      .processTemplates(templateIterator)
      .foreach(t => t.allReads.foreach(out.write))

    // close streams
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
