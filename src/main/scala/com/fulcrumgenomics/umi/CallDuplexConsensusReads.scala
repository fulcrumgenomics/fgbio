/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.sopt.clp
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.umi.VanillaUmiConsensusCallerOptions._
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.ProgressLogger

@clp(description =
  """
    |Calls duplex consensus sequences from reads generated from the same _double-stranded_ source molecule. Prior
    |to running this tool, read must have been grouped with `GroupReadsByUmi` using the `paired` strategy. Doing
    |so will apply (by default) MI tags to all reads of the form `*/A` and `*/B` where the /A and /B suffixes
    |with the same identifier denote reads that are derived from opposite strands of the same source duplex molecule.
    |
    |Reads from the same unique molecule are first partitioned by source strand and assembled into single
    |strand consensus molecules as described by CallMolecularConsensusReads.  Subsequently, for molecules that
    |have at least one observation of each strand, duplex consensus reads are assembled by combining the evidence
    |from the two single strand consensus reads.
    |
    |Because of the nature of duplex sequencing, this tool does not support fragment reads - if found in the
    |input they are _ignored_.  Similarly, read pairs for which consensus reads cannot be generated for one or
    |other read (R1 or R2) are omitted from the output.
    |
    |Consensus reads have a number of additional optional tags set in the resulting BAM file.  The tag names follow
    |a pattern where the first letter (a, b or c) denotes that the tag applies to the first single strand consensus (a),
    |second single-strand consensus (b) or the final duplex consensus (c).  The second letter is intended to capture
    |the meaning of the tag (e.g. d=depth, m=min depth, e=errors/error-rate) and is upper case for values that are
    |one per read and lower case for values that are one per base.
    |
    |The tags break down into those that are single-valued per read:
    |
    |```
    |consensus depth      [aD,bD,cD] (int)  : the maximum depth of raw reads at any point in the consensus reads
    |consensus min depth  [aM,bM,cM] (int)  : the minimum depth of raw reads at any point in the consensus reads
    |consensus error rate [aE,bE,cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls
    |```
    |
    |And those that have a value per base (duplex values are not generated, but can be generated by summing):
    |
    |```
    |consensus depth  [ad,bd] (short[]): the count of bases contributing to each single-strand consensus read at each position
    |consensus errors [ae,be] (short[]): the count of bases from raw reads disagreeing with the final single-strand consensus base
    |consensus errors [ac,bc] (string): the single-strand consensus bases
    |consensus errors [aq,bq] (string): the single-strand consensus qualities
    |```
    |
    |The per base depths and errors are both capped at 32,767. In all cases no-calls (Ns) and bases below the
    |min-input-base-quality are not counted in tag value calculations.
    |
    |The --min-reads option can take 1-3 values similar to `FilterConsensusReads`. For example:
    |
    |```
    |CallDuplexConsensusReads ... --min-reads 10 5 3
    |```
    |
    |If fewer than three values are supplied, the last value is repeated (i.e. `5 4` -> `5 4 4` and `1` -> `1 1 1`.  The
    |first value applies to the final consensus read, the second value to one single-strand consensus, and the last
    |value to the other single-strand consensus. It is required that if values two and three differ,
    |the _more stringent value comes earlier_.
  """,
  group = ClpGroups.Umi)
class CallDuplexConsensusReads
(@arg(flag='i', doc="The input SAM or BAM file.") val input: PathToBam,
 @arg(flag='o', doc="Output SAM or BAM file to write consensus reads.") val output: PathToBam,
 @arg(flag='p', doc="The prefix all consensus read names") val readNamePrefix: Option[String] = None,
 @arg(flag='R', doc="The new read group ID for all the consensus reads.") val readGroupId: String = "A",
 @arg(flag='1', doc="The Phred-scaled error rate for an error prior to the UMIs being integrated.") val errorRatePreUmi: PhredScore = DefaultErrorRatePreUmi,
 @arg(flag='2', doc="The Phred-scaled error rate for an error post the UMIs have been integrated.") val errorRatePostUmi: PhredScore = DefaultErrorRatePostUmi,
 @arg(flag='m', doc="Ignore bases in raw reads that have Q below this value.") val minInputBaseQuality: PhredScore = DefaultMinInputBaseQuality,
 @arg(flag='t', doc="If true, quality trim input reads in addition to masking low Q bases.") val trim: Boolean = false,
 @arg(flag='S', doc="The sort order of the output, if `:none:` then the same as the input.") val sortOrder: Option[SamOrder] = Some(SamOrder.Queryname),
 @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of input reads to a consensus read.") val minReads: Seq[Int] = Seq(1),
 @arg(doc="""
            |The maximum number of reads to use when building a single-strandconsensus. If more than this many reads are
            |present in a tag family, the family is randomly downsampled to exactly max-reads reads.
          """)
 val maxReads: Option[Int] = None,
 @arg(doc="The number of threads to use while consensus calling.") val threads: Int = 1,
) extends FgBioTool with LazyLogging {

  private val maxRecordsInRamPerThread = 128000

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  validate(errorRatePreUmi  > 0, "Phred-scaled error rate pre UMI must be > 0")
  validate(errorRatePostUmi > 0, "Phred-scaled error rate post UMI must be > 0")

  override def execute(): Unit = {
    val in  = SamSource(input)
    UmiConsensusCaller.checkSortOrder(in.header, input, logger.warning, fail)

    // The output file is unmapped, so for now let's clear out the sequence dictionary & PGs
    val outHeader = UmiConsensusCaller.outputHeader(in.header, readGroupId, sortOrder)
    val out = SamWriter(output, outHeader, sort=sortOrder)

    val caller = new DuplexConsensusCaller(
      readNamePrefix      = readNamePrefix.getOrElse(UmiConsensusCaller.makePrefixFromSamHeader(in.header)),
      readGroupId         = readGroupId,
      minInputBaseQuality = minInputBaseQuality,
      trim                = trim,
      errorRatePreUmi     = errorRatePreUmi,
      errorRatePostUmi    = errorRatePostUmi,
      minReads            = minReads,
      maxReads            = maxReads.getOrElse(VanillaUmiConsensusCallerOptions.DefaultMaxReads),
    )
    val progress = ProgressLogger(logger, unit=1000000)
    val iterator = new ConsensusCallingIterator(in.toIterator, caller, Some(progress), threads, maxRecordsInRamPerThread)
    out ++= iterator
    progress.logLast()

    in.safelyClose()
    out.close()
    caller.logStatistics(logger)
  }
}
