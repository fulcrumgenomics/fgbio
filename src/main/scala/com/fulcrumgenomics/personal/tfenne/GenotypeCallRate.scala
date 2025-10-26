/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToVcf}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.personal.tfenne.GenotypeCallRate.GenotypingMetric
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}

import scala.collection.mutable

object GenotypeCallRate {
  case class GenotypingMetric(val sample: String,
                              var total_sites: Long = 0,
                              var no_calls: Long = 0,
                              var low_gq: Long = 0,
                              var call_rate: Double = 0.0) extends Metric {

    def calculateDerived(): Unit = {
      this.call_rate = 1 - ((no_calls + low_gq) / total_sites.toDouble)
    }
  }
}

@clp(group=ClpGroups.Personal, description=
  """
    |Calculates per-sample genotype call rate over one or more VCFs.
  """)
class GenotypeCallRate(
  @arg(flag='i', doc="Input VCFs.") val input: Seq[PathToVcf],
  @arg(flag='o', doc="Output stats.") val output: FilePath,
  @arg(flag='q', doc="Minimum genotype quality to count as a call") val minGq: Int = 20)
  extends FgBioTool {

  /** All tools should implement this method. */
  override def execute(): Unit = {
    val allStats = mutable.HashMap[String, GenotypingMetric]()
    val progress = ProgressLogger(logger, noun="variant record", unit=10e6.toInt)

    input.foreach { vcf =>
      logger.info(s"Processing ${vcf}")
      val in = VcfSource(vcf)
      in.header.samples.filterNot(allStats.contains).foreach(s => allStats(s) = GenotypingMetric(s))

      in
        .filter(v => v.filters.isEmpty || v.filters == Variant.PassingFilters)
        .foreach { rec =>
          rec.genotypes.foreach { case (sample, gt) =>
            val stats = allStats(sample)
            stats.total_sites += 1

            if (gt.isNoCall) {
              stats.no_calls += 1
            }
            else if (gt.getOrElse[Int]("GQ", 0) < minGq) {
              stats.low_gq += 1
            }

            progress.record(rec)
          }
        }
    }

    val sortedStats = allStats.values.toIndexedSeq.sortBy(_.sample)
    sortedStats.foreach(_.calculateDerived())
    Metric.write(output,  sortedStats)
  }
}
