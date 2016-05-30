/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.wgs

import com.fulcrumgenomics.cmdline.{FgBioTool, ClpGroups}
import com.fulcrumgenomics.util.ProgressLogger
import dagr.commons.CommonsDef._
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import htsjdk.samtools.util.{IntervalList, Interval, SamLocusIterator, CloserUtil}
import htsjdk.samtools.{SAMFileHeader, SamReaderFactory, SamReader}

import scala.collection.JavaConversions._
import scala.collection.mutable.ListBuffer

/**
  * Creates a list of intervals for WGS variant calling.
  */
@clp(description =
  """
    |Creates a list of intervals for WGS variant calling.
    |
    |The output intervals will only contain sites that have both a minimum and maximum coverage as set by the user.
    |This tool is useful for when sites with too much coverage cause too much slow down in variant calling.
  """,
  group = ClpGroups.Utilities)
class CreateWgsCallableIntervalList
( @arg(flag="i", doc="Input BAM file.") val input: PathToBam,
  @arg(flag="o", doc="Output interval list") val output: PathToIntervals,
  @arg(flag="m", doc="Minimum coverage.") val minCoverage: Int = 1,
  @arg(flag="M", doc="Maximum coverage.") val maxCoverage: Int = 1024
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val progress: ProgressLogger = new ProgressLogger(logger, noun="positions", unit=10*1000*1000)
    val in: SamReader = SamReaderFactory.make.open(input.toFile)
    val locusIterator = new SamLocusIterator(in)
    locusIterator.setMaxReadsToAccumulatePerLocus(maxCoverage+1) // +1 so we know when we have hit the maximum
    if (0 < minCoverage) locusIterator.setEmitUncoveredLoci(false)

    val converter = new CoveredLociToIntervals(header=in.getFileHeader, minCoverage=minCoverage, maxCoverage=maxCoverage)
    locusIterator.iterator().foreach { locus =>
      converter.add(locus)
      progress.record(locus.getSequenceName, locus.getPosition)
    }
    converter.intervalList.write(output.toFile)

    CloserUtil.close(locusIterator)
  }
}

private class CoveredLociToIntervals(val header: SAMFileHeader, val minCoverage: Int, val maxCoverage: Int) {

  private var counter: Int = 0
  private var refName: String = ""
  private var refIndex: Int = -1
  private var start: Int = -1
  private var end: Int = -1
  private val intervals = ListBuffer[Interval]()

  def add(locus: LocusInfo): Boolean = {
    val coverage = locus.getRecordAndPositions.size()
    if (coverage < minCoverage || maxCoverage < coverage) return false

    if (locus.getSequenceIndex == refIndex && locus.getPosition == end + 1) {
      end = end + 1
    }
    else {
      if (0 <= refIndex) intervals.append(new Interval(refName, start, end, false, s"$counter"))

      refName  = locus.getSequenceName
      refIndex = locus.getSequenceIndex
      start    = locus.getPosition
      end      = locus.getPosition
      counter  = counter + 1
    }

    true
  }

  def intervalList: IntervalList = {
    if (0 <= refIndex) intervals.append(new Interval(refName, start, end, false, s"$counter"))

    val outHeader = new SAMFileHeader()
    outHeader.setSequenceDictionary(header.getSequenceDictionary)

    val intervalList = new IntervalList(outHeader)
    intervalList.addall(intervals.toList)
    intervalList
  }
}
