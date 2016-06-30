/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.umi.ConsensusCallerOptions._
import com.fulcrumgenomics.util.ProgressLogger
import dagr.commons.CommonsDef.PathToBam
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import dagr.sopt.cmdline.ValidationException
import htsjdk.samtools._
import htsjdk.samtools.util.CloserUtil
import com.fulcrumgenomics.FgBioDef._

import scala.collection.JavaConverters._


@clp(description =
  """
    |Calls consensus sequences from reads with the same unique molecular tag.
    |
    |Reads with the same unique molecular tag are examined base-by-base to assess the likelihood of each base in the
    |source molecule.  The likelihood model is as follows:
    |   1. First, the base qualities are adjusted. The base qualities are assumed to represent the probability of a
    |      sequencing error (i.e. the sequencer observed the wrong base present on the cluster/flowcell/well).  A fixed
    |      value is subtracted from the phred-scaled base qualities (ex. Q30 with a shift of 10 becomes Q20).  Next, the
    |      base qualities are capped to a maximum phred-scaled value.  Finally the base quality are converted to a
    |      probability to incorporate a probability representing the chance of an error from the time the unique
    |      molecular tags were integrated to just prior to sequencing.  The resulting probability is the error rate of
    |      all processes from right after integrating the molecular tag through to the end of sequencing.
    |   2. Next, a consensus sequence is called for all reads with the same unique molecular tag base-by-base.  For a
    |      given base position in the reads, the likelihoods that an A, C, G, or T is the base for the underlying
    |      source molecule respectively are computed by multiplying the likelihood of each read observing the base
    |      position being considered.  The probability of error (from 1.) is used when the observed base does not match
    |      the hypothesized base for the underlying source molecule, while one minus that probability is used otherwise.
    |      The computed likelihoods are normalized by dividing them by the sum of all four likelihoods to produce a
    |      posterior probability, namely the probability that the source molecule was an A, C, G, or T from just after
    |      integrating molecular tag through to sequencing, given the observations.  The base with the maximum posterior
    |      probability as the consensus call, and the posterior probability is used as its raw base quality.
    |   3. Finally, the consensus raw base quality is modified by incorporating the probability of an error prior to
    |      integrating the unique molecular tags.  Therefore, the probability used for the final consensus base
    |      quality is the posterior probability of the source molecule having the consensus base given the observed
    |      reads with the same molecular tag, all the way from sample extraction and through sample and library
    |      preparation, through preparing the library for sequencing (ex. amplification, target selection), and finally,
    |      through sequencing.
    |
    |This tool assumes that reads with the same tag are grouped together (consecutive in the file). Also, this tool
    |calls each end of a pair independently, and does not jointly call bases that overlap within a pair.  Insertion or
    |deletion errors in the reads are not considered in the consensus model.
  """,
  group = ClpGroups.Umi)
class CallMolecularConsensusReads
( @arg(flag="i", doc="The input SAM or BAM file.") val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file to write consensus reads.") val output: PathToBam,
  @arg(flag="r", doc="Output SAM or BAM file to write reads not used.") val rejects: PathToBam,
  @arg(flag="t", doc="The SAM attribute with the unique molecule tag.") val tag: String = DefaultTag,
  @arg(flag="p", doc="The Prefix all consensus read names") val readNamePrefix: Option[String] = None,
  @arg(flag="R", doc="The new read group ID for all the consensus reads.") val readGroupId: String = "A",
  @arg(flag="1", doc="The Phred-scaled error rate for an error prior to the UMIs being integrated.") val errorRatePreUmi: PhredScore = DefaultErrorRatePreUmi,
  @arg(flag="2", doc="The Phred-scaled error rate for an error post the UMIs have been integrated.") val errorRatePostUmi: PhredScore = DefaultErrorRatePostUmi,
  @arg(flag="q", doc="Cap the maximum base quality in the input (after shifting).") val maxBaseQuality: Double = 40,
  @arg(flag="s", doc="Subtract this base quality from the input base qualities (prior to capping).") val baseQualityShift: PhredScore = DefaultBaseQualityShift,
  @arg(flag="N", doc="Mask (make 'N') consensus bases with quality less than this threshold.") val minConsensusBaseQuality: PhredScore = DefaultMinConsensusBaseQuality,
  @arg(flag="M", doc="The minimum number of reads to produce a consensus base.") val minReads: Int = DefaultMinReads,
  @arg(flag="Q", doc="The minimum mean base quality across a consensus base to output.") val minMeanConsensusBaseQuality: PhredScore = DefaultMinMeanConsensusBaseQuality,
  @arg(flag="P", doc="Require a consensus call for both ends of a pair if true.") val requireConsensusForBothPairs: Boolean = DefaultRequireConsensusForBothPairs
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Seq(output, rejects).foreach(Io.assertCanWriteFile(_))
  if (tag.length != 2)      throw new ValidationException("attribute must be of length 2")
  if (errorRatePreUmi < 0)  throw new ValidationException("Phred-scaled error rate pre UMI must be >= 0")
  if (errorRatePostUmi < 0) throw new ValidationException("Phred-scaled error rate post UMI must be >= 0")

  /** Main method that does the work of reading input files, creating the consensus reads, and writing the output file. */
  override def execute(): Unit = {
    val in  = SamReaderFactory.make().open(input.toFile)
    val out = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, null)
    val rej = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, rejects.toFile, null)
    // TODO: metrics...

    val options = new ConsensusCallerOptions(
      tag                          = tag,
      errorRatePreUmi              = errorRatePreUmi,
      errorRatePostUmi             = errorRatePostUmi,
      maxBaseQuality               = maxBaseQuality,
      baseQualityShift             = baseQualityShift,
      minConsensusBaseQuality      = minConsensusBaseQuality,
      minReads                     = minReads,
      minMeanConsensusBaseQuality  = minMeanConsensusBaseQuality,
      requireConsensusForBothPairs = requireConsensusForBothPairs
    )

    val progress = new ProgressLogger(logger, unit=10000)
    val consensusCaller = new ConsensusCaller(
      input          = in.iterator().asScala,
      header         = in.getFileHeader,
      readNamePrefix = readNamePrefix,
      readGroupId    = readGroupId,
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
