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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathToBam, PathToFasta, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.collection.ParIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.commons.util.Threads.IterableThreadLocal
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger}
import enumeratum.EnumEntry

import scala.collection.immutable


@clp(group = ClpGroups.SamOrBam, description=
  """
    |Consensus calls overlapping bases in read pairs.
    |
    |## Inputs and Outputs
    |
    |In order to correctly correct reads by template, the input BAM must be either `queryname` sorted or `query` grouped.
    |The sort can be done in streaming fashion with:
    |
    |```
    |samtools sort -n -u in.bam | fgbio CallOverlappingConsensusBases -i /dev/stdin ...
    |```
    |
    |The output sort order may be specified with `--sort-order`.  If not given, then the output will be in the same
    |order as input.
    |
    |The reference FASTA must be given so that any existing `NM`, `UQ` and `MD` tags can be repaired.
    |
    |## Correction
    |
    |Only mapped read pairs with overlapping bases will be eligible for correction.
    |
    |Each read base from the read and its mate that map to same position in the reference will be used to create
    |a consensus base as follows:
    |
    |1. If the base agree, then the chosen agreement strategy (`--agreement-strategy`) will be used.
    |2. If the base disagree, then the chosen disagreement strategy (`--disagreement-strategy`) will be used.
    |
    |The agreement strategies are as follows:
    |
    |* Consensus:   Call the consensus base and return a new base quality that is the sum of the two base qualities.
    |* MaxQual:     Call the consensus base and return a new base quality that is the maximum of the two base qualities.
    |* PassThrough: Leave the bases and base qualities unchanged.
    |
    |In the context of disagreement strategies, masking a base will make the base an "N" with base quality phred-value "2".
    |The disagreement strategies are as follows:
    |
    |* MaskBoth:      Mask both bases.
    |* MaskLowerQual: Mask the base with the lowest base quality, with the other base unchanged.  If the base qualities
    |                 are the same, mask both bases.
    |* Consensus:     Consensus call the base.  If the base qualities are the same, mask both bases.  Otherwise, call the
    |                 base with the highest base quality and return a new base quality that is the difference between the
    |                 highest and lowest base quality.
 |  """)
class CallOverlappingConsensusBases
(@arg(flag='i', doc="Input SAM or BAM file of aligned reads.") val input: PathToBam,
 @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Output metrics file.") val metrics: FilePath,
 @arg(flag='r', doc="Reference sequence fasta file.") val ref: PathToFasta,
 @arg(doc="The number of threads to use while consensus calling.") val threads: Int = 1,
 @arg(flag='S', doc="The sort order of the output. If not given, output will be in the same order as input if the input.")
 val sortOrder: Option[SamOrder] = None,
 @arg(doc="The strategy to consensus call when both bases agree.  See the usage for more details")
 val agreementStrategy: AgreementStrategy = AgreementStrategy.Consensus,
 @arg(doc="The strategy to consensus call when both bases disagree.  See the usage for more details")
 val disagreementStrategy: DisagreementStrategy = DisagreementStrategy.Consensus
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  private case class ThreadData
  (caller: OverlappingBasesConsensusCaller             = new OverlappingBasesConsensusCaller(agreementStrategy=agreementStrategy, disagreementStrategy=disagreementStrategy),
   templateMetric: CallOverlappingConsensusBasesMetric = CallOverlappingConsensusBasesMetric(kind=CountKind.Templates),
   basesMetric: CallOverlappingConsensusBasesMetric    = CallOverlappingConsensusBasesMetric(kind=CountKind.Bases)
  )

  override def execute(): Unit = {
    val source           = SamSource(input)
    val outSort          = sortOrder.flatMap { order => if (SamOrder(source.header).contains(order)) None else Some(order) }
    val writer           = Bams.nmUqMdTagRegeneratingWriter(writer=SamWriter(output, source.header.clone(), sort=outSort), ref=ref)
    val progress         = new ProgressLogger(logger)
    val templateIterator = Bams.templateIterator(source)
    val threadData       = new IterableThreadLocal(() => ThreadData())

    // Require queryname sorted or query grouped
    Bams.requireQueryGrouped(header=source.header, toolName="CallOverlappingConsensusBases")

    ParIterator(templateIterator, threads=threads)
      .map { template =>
        val threadDatum = threadData.get()
        threadDatum.synchronized {
          // update metrics
          threadDatum.templateMetric.total += 1
          threadDatum.basesMetric.total += template.primaryReads.map(_.length).sum
          // corrects
          val stats          = threadDatum.caller.call(template)
          val correctedBases = stats.r1CorrectedBases + stats.r2CorrectedBases
          if (stats.overlappingBases > 0) {
            threadDatum.templateMetric.overlapping += 1
            threadDatum.basesMetric.overlapping += stats.overlappingBases
            if (correctedBases > 0) {
              threadDatum.templateMetric.corrected += 1
              threadDatum.basesMetric.corrected += correctedBases
            }
          }
        }
        template.allReads.foreach(progress.record)
        template
      }.toAsync(ParIterator.DefaultChunkSize * 8).foreach { template =>
      writer ++= template.allReads
    }
    progress.logLast()
    source.safelyClose()
    writer.close()

    val templatesMetric = CallOverlappingConsensusBasesMetric(kind=CountKind.Templates)
    val basesMetric     = CallOverlappingConsensusBasesMetric(kind=CountKind.Bases)
    threadData.foreach { datum =>
      templatesMetric += datum.templateMetric
      basesMetric     += datum.basesMetric
    }

    Metric.write(metrics, templatesMetric, basesMetric)
  }
}

/** Collects the the number of reads or bases that were examined, had overlap, and were corrected as part of
  * the [[CallOverlappingConsensusBases]] tool.
  *
  * @param kind template if the counts are per template, bases if counts are in units of bases.
  * @param total the total number of templates (bases) examined
  * @param overlapping the total number of templates (bases) that were overlapping
  * @param corrected the total number of templates (bases) that were corrected.
  */
case class CallOverlappingConsensusBasesMetric
(
  kind: CountKind,
  var total: Long = 0,
  var overlapping: Long = 0,
  var corrected: Long = 0,
) extends Metric {
  def +=(other: CallOverlappingConsensusBasesMetric): CallOverlappingConsensusBasesMetric = {
    require(this.kind == other.kind)
    this.total += other.total
    this.overlapping += other.overlapping
    this.corrected += other.corrected
    this
  }
}

sealed trait CountKind extends EnumEntry

/** Enumeration for the type of counts in [[CallOverlappingConsensusBasesMetric]]. */
object CountKind extends FgBioEnum[CountKind] {
  case object Templates extends CountKind
  case object Bases extends CountKind

  override def values: immutable.IndexedSeq[CountKind] = findValues
}