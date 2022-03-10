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

package com.fulcrumgenomics.personal.nhomer

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
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
    |Corrects overlapping bases in FR read pairs.
    |
    |To examine reads from the same template, the input is resorted if neither queryname sorted nor query grouped.
  """)
class CorrectOverlappingBases
( @arg(flag='i', doc="Input SAM or BAM file of aligned reads.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='m', doc="Output metrics file.") val metrics: FilePath,
  @arg(doc="The number of threads to use while consensus calling.") val threads: Int = 1
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  private case class ThreadData(caller: OverlappingBasesConsensusCaller, templateMetric: CorrectOverlappingBasesMetric, basesMetric: CorrectOverlappingBasesMetric)
  private object ThreadData {
    def apply(): ThreadData = ThreadData(
      caller         = new OverlappingBasesConsensusCaller(),
      templateMetric = CorrectOverlappingBasesMetric(tpe=CountType.Templates),
      basesMetric    = CorrectOverlappingBasesMetric(tpe=CountType.Bases)
    )
  }

  override def execute(): Unit = {
    val source           = SamSource(input)
    val writer           = SamWriter(output, source.header)
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
          val stats = threadDatum.caller.correct(template)
          val overlappingBases = stats.r1OverlappingBases + stats.r2OverlappingBases
          val correctedBases   = stats.r1CorrectedBases + stats.r2CorrectedBases
          if (overlappingBases > 0) {
            threadDatum.templateMetric.overlapping += 1
            threadDatum.basesMetric.overlapping += overlappingBases
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

    val templatesMetric = CorrectOverlappingBasesMetric(tpe=CountType.Templates)
    val basesMetric     = CorrectOverlappingBasesMetric(tpe=CountType.Bases)
    threadData.foreach { datum =>
      templatesMetric += datum.templateMetric
      basesMetric     += datum.basesMetric
    }

    Metric.write(metrics, templatesMetric, basesMetric)
  }

}

case class CorrectOverlappingBasesMetric
(
  tpe: CountType,
  var total: Long = 0,
  var overlapping: Long = 0,
  var corrected: Long = 0,
) extends Metric {
  def +=(other: CorrectOverlappingBasesMetric): CorrectOverlappingBasesMetric = {
    require(this.tpe == other.tpe)
    this.total += other.total
    this.overlapping += other.overlapping
    this.corrected += other.corrected
    this
  }
}

sealed trait CountType extends EnumEntry

object CountType extends FgBioEnum[CountType] {
  case object Templates extends CountType
  case object Bases extends CountType

  override def values: immutable.IndexedSeq[CountType] = findValues
}