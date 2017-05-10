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

package com.fulcrumgenomics.bam

import java.nio.file.Path
import java.text.DecimalFormat

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.ProgressLogger
import com.fulcrumgenomics.commons.CommonsDef.{PathToBam, PathToIntervals}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import htsjdk.samtools._
import htsjdk.samtools.util.{CloserUtil, IntervalList}

/**
  * Program which takes in a BAM file and filters out all reads for templates that match one or more
  * criteria.  Designed to be used to filter out reads that might confuse variant callers and lead
  * to false positive variant calls.
  *
  * @author Tim Fennell
  */
@clp(description = 
  """
     |Filters reads out of a BAM file. Removes reads that may not be useful in downstream processing, in order
     |to reduce the size of the file. By default will remove unmapped reads, reads with MAPQ=0, reads
     |marked as secondary alignments, reads marked as duplicates, and if a set of Intervals are provided,
     |reads that do not overlap any of the intervals.
     |
     |NOTE: this will usually produce a BAM file in which some mate-pairs are orphaned (i.e. read 1 or
     |read 2 is included, but not both), but does not update any flag fields.
  """,
  group = ClpGroups.SamOrBam)
class FilterBam
( @arg(flag='i', doc="Input BAM file.")                                           val input: PathToBam,
  @arg(flag='o', doc="Output BAM file.")                                          val output: PathToBam,
  @arg(flag='l', doc="Optionally remove reads not overlapping intervals.")        val intervals: Option[PathToIntervals] = None,
  @arg(flag='D', doc="If true remove all reads that are marked as duplicates.")   val removeDuplicates: Boolean = true,
  @arg(flag='U', doc="Remove all unmapped reads.")                                val removeUnmappedReads: Boolean = true,
  @arg(flag='M', doc="Remove all mapped reads with MAPQ lower than this number.") val minMapQ: Int = 1,
  @arg(flag='S', doc="Remove all reads marked as secondary alignments.")          val removeSecondaryAlignments: Boolean = true
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertReadable)

  override def execute(): Unit = {
    val progress = new ProgressLogger(logger, verb="written", unit=5e6.toInt)
    val in       = SamReaderFactory.make.open(input)
    val iterator = buildInputIterator(in, intervals)
    val out      = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, null)
    val kept = iterator.count { rec => {
        val throwOut = (removeDuplicates && rec.getDuplicateReadFlag) ||
          (removeUnmappedReads && rec.getReadUnmappedFlag) ||
          (!rec.getReadUnmappedFlag && rec.getMappingQuality < minMapQ) ||
          (removeSecondaryAlignments && !rec.getReadUnmappedFlag && rec.getNotPrimaryAlignmentFlag)

        if (throwOut) {
          false
        }
        else {
          out.addAlignment(rec)
          progress.record(rec)
          true
        }
      }
    }
    logger.info("Kept " + new DecimalFormat("#,##0").format(kept) + " records.")
    out.close()
    iterator.safelyClose()
  }
  /**
    * If intervalListFile is null return an interator over all the input, otherwise returns an
    * iterator over only those reads that overlap one or more of the intervals in the file.
    */
  protected def buildInputIterator(in: SamReader, intervalListFile: Option[Path]): SAMRecordIterator = {
    intervalListFile match {
      case None => in.iterator()
      case Some(file) =>
        val intervals = IntervalList.fromFile(file.toFile).uniqued
        val dict      = intervals.getHeader.getSequenceDictionary
        val qs        = intervals.getIntervals.map(interval =>
          new QueryInterval(dict.getSequenceIndex(interval.getContig), interval.getStart, interval.getEnd)).toArray
        in.queryOverlapping(qs)
    }
  }
}
