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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.{LazyLogging, LogLevel, Logger}
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}

@clp(description =
  """
    |Calls consensus sequences from reads with the same unique molecular tag.
    |
    |Reads with the same unique molecular tag are examined base-by-base to assess the likelihood of each base in the
    |source molecule.  The likelihood model is as follows:
    |
    |1. First, the base qualities are adjusted. The base qualities are assumed to represent the probability of a
    |   sequencing error (i.e. the sequencer observed the wrong base present on the cluster/flowcell/well). The base
    |   quality scores are converted to probabilities incorporating a probability representing the chance of an error
    |   from the time the unique molecular tags were integrated to just prior to sequencing.  The resulting probability
    |   is the error rate of all processes from right after integrating the molecular tag through to the end of
    |   sequencing.
    |2. Next, a consensus sequence is called for all reads with the same unique molecular tag base-by-base.  For a
    |   given base position in the reads, the likelihoods that an A, C, G, or T is the base for the underlying
    |   source molecule respectively are computed by multiplying the likelihood of each read observing the base
    |   position being considered.  The probability of error (from 1.) is used when the observed base does not match
    |   the hypothesized base for the underlying source molecule, while one minus that probability is used otherwise.
    |   The computed likelihoods are normalized by dividing them by the sum of all four likelihoods to produce a
    |   posterior probability, namely the probability that the source molecule was an A, C, G, or T from just after
    |   integrating molecular tag through to sequencing, given the observations.  The base with the maximum posterior
    |   probability as the consensus call, and the posterior probability is used as its raw base quality.
    |3. Finally, the consensus raw base quality is modified by incorporating the probability of an error prior to
    |   integrating the unique molecular tags.  Therefore, the probability used for the final consensus base
    |   quality is the posterior probability of the source molecule having the consensus base given the observed
    |   reads with the same molecular tag, all the way from sample extraction and through sample and library
    |   preparation, through preparing the library for sequencing (e.g. amplification, target selection), and finally,
    |   through sequencing.
    |
    |This tool assumes that reads with the same tag are grouped together (consecutive in the file). Also, this tool
    |calls each end of a pair independently, and does not jointly call bases that overlap within a pair.  Insertion or
    |deletion errors in the reads are not considered in the consensus model.
    |
    |Particular attention should be paid to setting the `--min-reads` parameter as this can have a dramatic effect on
    |both results and runtime.  For libraries with low duplication rates (e.g. 100-300X exomes libraries) in which it
    |is desirable to retain singleton reads while making consensus reads from sets of duplicates, `--min-reads=1` is
    |appropriate.  For libraries with high duplication rates where it is desirable to only produce consensus reads
    |supported by 2+ reads to allow error correction, `--min-reads=2` or higher is appropriate.  After generation,
    |consensus reads can be further filtered using the _FilterConsensusReads_ tool.  As such it is always safe to run
    |with `--min-reads=1` and filter later, but filtering at this step can improve performance significantly.
    |
    |Consensus reads have a number of additional optional tags set in the resulting BAM file.  The tags break down into
    |those that are single-valued per read:
    |
    |```
    |consensus depth      [cD] (int)  : the maximum depth of raw reads at any point in the consensus read
    |consensus min depth  [cM] (int)  : the minimum depth of raw reads at any point in the consensus read
    |consensus error rate [cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls
    |```
    |
    |And those that have a value per base:
    |
    |```
    |consensus depth  [cd] (short[]): the count of bases contributing to the consensus read at each position
    |consensus errors [ce] (short[]): the number of bases from raw reads disagreeing with the final consensus base
    |```
    |
    |The per base depths and errors are both capped at 32,767. In all cases no-calls (`N`s) and bases below the
    |`--min-input-base-quality` are not counted in tag value calculations.
  """,
  group = ClpGroups.Umi)
class CallMolecularConsensusReads
(@arg(flag='i', doc="The input SAM or BAM file.") val input: PathToBam,
 @arg(flag='o', doc="Output SAM or BAM file to write consensus reads.") val output: PathToBam,
 @arg(flag='r', doc="Optional output SAM or BAM file to write reads not used.") val rejects: Option[PathToBam] = None,
 @arg(flag='t', doc="The SAM attribute with the unique molecule tag.") val tag: String = DefaultTag,
 @arg(flag='p', doc="The Prefix all consensus read names") val readNamePrefix: Option[String] = None,
 @arg(flag='R', doc="The new read group ID for all the consensus reads.") val readGroupId: String = "A",
 @arg(flag='1', doc="The Phred-scaled error rate for an error prior to the UMIs being integrated.") val errorRatePreUmi: PhredScore = DefaultErrorRatePreUmi,
 @arg(flag='2', doc="The Phred-scaled error rate for an error post the UMIs have been integrated.") val errorRatePostUmi: PhredScore = DefaultErrorRatePostUmi,
 @arg(flag='m', doc="Ignore bases in raw reads that have Q below this value.") val minInputBaseQuality: PhredScore = DefaultMinInputBaseQuality,
 @arg(flag='N', doc=
   """
     |Deprecated: will be removed in future versions; use FilterConsensusReads to filter consensus bases on
     |quality instead. Mask (make 'N') consensus bases with quality less than this threshold.
   """)
 val minConsensusBaseQuality: PhredScore = 2.toByte,
 @arg(flag='M', doc="The minimum number of reads to produce a consensus base.") val minReads: Int,
 @arg(doc="""
            |The maximum number of reads to use when building a consensus. If more than this many reads are
            |present in a tag family, the family is randomly downsampled to exactly max-reads reads.
          """)
 val maxReads: Option[Int] = None,
 @arg(flag='B', doc="If true produce tags on consensus reads that contain per-base information.") val outputPerBaseTags: Boolean = DefaultProducePerBaseTags,
 @arg(flag='S', doc="The sort order of the output, if `:none:` then the same as the input.") val sortOrder: Option[SamOrder] = Some(SamOrder.Queryname),
 @arg(flag='D', doc="Turn on debug logging.") val debug: Boolean = false,
 @arg(doc="The number of threads to use while consensus calling.") val threads: Int = 1
) extends FgBioTool with LazyLogging {

  if (debug) Logger.level = LogLevel.Debug

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  rejects.foreach(Io.assertCanWriteFile(_))

  if (tag.length != 2)      throw new ValidationException("attribute must be of length 2")
  if (errorRatePreUmi < 0)  throw new ValidationException("Phred-scaled error rate pre UMI must be >= 0")
  if (errorRatePostUmi < 0) throw new ValidationException("Phred-scaled error rate post UMI must be >= 0")

  /** Main method that does the work of reading input files, creating the consensus reads, and writing the output file. */
  override def execute(): Unit = {
    val in  = SamSource(input)
    UmiConsensusCaller.checkSortOrder(in.header, input, logger.warning, fail)
    val rej = rejects.map(r => SamWriter(r, in.header))

    // The output file is unmapped, so for now let's clear out the sequence dictionary & PGs
    val outHeader = UmiConsensusCaller.outputHeader(in.header, readGroupId, sortOrder)
    val out       = SamWriter(output, outHeader, sort=sortOrder)

    val options = new VanillaUmiConsensusCallerOptions(
      tag                          = tag,
      errorRatePreUmi              = errorRatePreUmi,
      errorRatePostUmi             = errorRatePostUmi,
      minInputBaseQuality          = minInputBaseQuality,
      minConsensusBaseQuality      = minConsensusBaseQuality,
      minReads                     = minReads,
      maxReads                     = maxReads.getOrElse(VanillaUmiConsensusCallerOptions.DefaultMaxReads),
      producePerBaseTags           = outputPerBaseTags
    )

    val caller = new VanillaUmiConsensusCaller(
      readNamePrefix = readNamePrefix.getOrElse(UmiConsensusCaller.makePrefixFromSamHeader(in.header)),
      readGroupId    = readGroupId,
      options        = options,
      rejects        = rej
    )

    val iterator = new ConsensusCallingIterator(in.iterator, caller, Some(ProgressLogger(logger, unit=5e5.toInt)), threads=threads)
    out ++= iterator

    in.safelyClose()
    out.close()
    rej.foreach(_.close())
    caller.logStatistics(logger)
  }
}
