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

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathToBam, SafelyClosable}
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
    |In order to correctly correct reads by template, the input BAM must be either `queryname` sorted or `query` grouped.  The
    |sort can be done in streaming fashion with:
    |
    |```
    |fgbio --compression 0 SortBam -i in.bam -o out.bam -s queryname | fgbio CallOverlappingConsensusBases -i /dev/stdin ...
    |```
    |
    |The output sort order may be specified with `--sort-order`.  If not given, then the output will be in the same
    |order as input.
    |
    |## Correction
    |
    |Only mapped read pairs with overlapping bases will be eligible for correction.
    |
    |Each read base from the read and its mate that map to same position in the reference will be used to create
    |a consensus base as follows:
    |
    |1. if the read and mate bases are the same, the consensus base is that base with the base quality equal to the sum
    |   of the two base qualities.  The base quality can be the maximum base quality of the two base qualities if
    |   `--max-qual-on-agreement` is used.
    |2. if the read and mate bases differ, then the base with the highest associated base quality will be the consensus
    |   call.  If the read and mate have the same base quality, then the output base quality will be 2.  Otherwise,
    |   the base quality will be the difference between the larger and smaller base quality.  The
    |   `--only-mask-disagreements` option overrides behavior and sets all differing bases to `N` with a base quality of
    |   2.
  """)
class CallOverlappingConsensusBases
(@arg(flag='i', doc="Input SAM or BAM file of aligned reads.") val input: PathToBam,
 @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Output metrics file.") val metrics: FilePath,
 @arg(doc="The number of threads to use while consensus calling.") val threads: Int = 1,
 @arg(flag='S', doc="The sort order of the output. If not given, output will be in the same order as input if the input.")
  val sortOrder: Option[SamOrder] = None,
 @arg(doc="""If the read and mate bases disagree at a given reference position, true to mask (make 'N') the read and mate
             |bases, otherwise pick the base with the highest base quality and return a base quality that's the difference
             |between the higher and lower base qualities.""")

  val maskDisagreements: Boolean = false,
 @arg(doc= """If the read and mate bases agree at a given reference position, true to for the resulting base quality
              |to be the maximum base quality, otherwise the sum of the base qualities.""")
  val maxQualOnAgreement: Boolean = false
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  private case class ThreadData
  (caller: OverlappingBasesConsensusCaller             = new OverlappingBasesConsensusCaller(maskDisagreements=maskDisagreements, maxQualOnAgreement=maxQualOnAgreement),
   templateMetric: CallOverlappingConsensusBasesMetric = CallOverlappingConsensusBasesMetric(kind=CountKind.Templates),
   basesMetric: CallOverlappingConsensusBasesMetric    = CallOverlappingConsensusBasesMetric(kind=CountKind.Bases)
  )

  override def execute(): Unit = {
    val source           = SamSource(input)
    val outSort          = sortOrder.flatMap { order => if (SamOrder(source.header).contains(order)) None else Some(order) }
    val writer           = SamWriter(output, source.header, sort=outSort)
    val progress         = new ProgressLogger(logger)
    val templateIterator = Bams.templateIterator(source)
    val threadData       = new IterableThreadLocal(() => ThreadData())

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