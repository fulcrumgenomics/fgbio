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

import com.fulcrumgenomics.bam.ConsensusCallerOptions._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{PhredValue, ProgressLogger}
import dagr.commons.CommonsDef.PathToBam
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import dagr.sopt.cmdline.ValidationException
import htsjdk.samtools._
import htsjdk.samtools.util.CloserUtil

import scala.collection.JavaConverters._


@clp(description =
  """
    |Calls consensus sequences from reads with the same UMI identifier.
    |
    |This tool assumes that reads with the same identifier are grouped together (consecutive in the file).
  """,
  group = ClpGroups.SamOrBam)
class CallConsensusFromUmis
(  @arg(flag="i", doc="The input SAM or BAM file.") val input: PathToBam,
   @arg(flag="o", doc="Output SAM or BAM file to write consensus reads.") val output: PathToBam,
   @arg(flag="r", doc="Output SAM or BAM file to write reads not used.") val rejects: PathToBam,
   @arg(flag="t", doc="The SAM attribute with the unique UMI identifier.") val attribute: String = DefaultAttribute,
   @arg(flag="p", doc="The Prefix all consensus read names") val readNamePrefix: Option[String] = None,
   @arg(flag="R", doc="The Read Group for all the consensus reads.") val readGroup: String = "A",
   @arg(flag="1", doc="The Phred-scaled error rate for an error prior to the UMIs being integrated.") val errorRatePreUmi: PhredValue = DefaultErrorRatePreUmi,
   @arg(flag="2", doc="The Phred-scaled error rate for an error post the UMIs have been integrated.") val errorRatePostUmi: PhredValue = DefaultErrorRatePostUmi,
   @arg(flag="q", doc="Cap the maximum base quality in the input.") val maxBaseQuality: PhredValue = DefaultMaxBaseQuality,
   @arg(flag="s", doc="Subtract this base quality from the input base qualities.") val baseQualityShift: PhredValue = DefaultBaseQualityShift,
   @arg(flag="N", doc="Mask (make 'N') consensus bases with quality less than this threshold.") val minConsensusBaseQuality: PhredValue = DefaultMinConsensusBaseQuality,
   @arg(flag="M", doc="The minimum number of reads to produce a consensus base.") val minReads: Int = DefaultMinReads,
   @arg(flag="Q", doc="The minimum mean base quality across a consensus base to output.") val minMeanConsensusBaseQuality: PhredValue = DefaultMinMeanConsensusBaseQuality
  ) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Seq(output, rejects).foreach(Io.assertCanWriteFile(_))
  if (attribute.length != 2) throw new ValidationException("attribute must be of length 2")
  if (errorRatePreUmi < 0)   throw new ValidationException("Error rate pre UMI must be >= 0")
  if (errorRatePostUmi < 0)  throw new ValidationException("Error rate post UMI must be >= 0")

  /** TODO
    * - what about when a read has an indel causing the rest of the bases to mismatch?
    * - what about if a read is soft-clipped and the rest of the bases mismatch (ex. 30M120S) [second pass]?
    *
    * TODO: hard clipping?
    *
    */

  /** Main method that does the work of reading input files, creating the consensus reads, and writing the output file. */
  override def execute(): Unit = {
    val in  = SamReaderFactory.make().open(input.toFile)
    val out = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, null)
    val rej = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, rejects.toFile, null)
    // TODO: metrics...

    val options = new ConsensusCallerOptions(
      attribute                   = attribute,
      errorRatePreUmi             = errorRatePreUmi,
      errorRatePostUmi            = errorRatePostUmi,
      maxBaseQuality              = maxBaseQuality,
      baseQualityShift            = baseQualityShift,
      minConsensusBaseQuality     = minConsensusBaseQuality,
      minReads                    = minReads,
      minMeanConsensusBaseQuality = minMeanConsensusBaseQuality
    )

    val progress = new ProgressLogger(logger)
    val consensusCaller = new ConsensusCaller(
      input          = in.iterator().asScala,
      header         = in.getFileHeader,
      readNamePrefix = readNamePrefix,
      readGroup      = readGroup,
      options        = options,
      rejects        = Some(rej),
      progress       = Some(progress)
    )

    consensusCaller.foreach { rec => out.addAlignment(rec) }

    CloserUtil.close(in)
    out.close()
    rej.close()
    logger.info(s"Processed ${progress.getCount} records.")
  }
}

// TODO
class ConsensusCallingMetricsCollector {
}
