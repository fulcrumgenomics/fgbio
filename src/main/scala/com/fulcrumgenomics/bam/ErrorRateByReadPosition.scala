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
 *
 */

package com.fulcrumgenomics.bam

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Metric.Count
import com.fulcrumgenomics.util.{Metric, ProgressLogger, Rscript}
import com.fulcrumgenomics.vcf.{ByIntervalListVariantContextIterator, VariantMask}
import htsjdk.samtools.filter.{DuplicateReadFilter, FailsVendorReadQualityFilter, SamRecordFilter, SecondaryOrSupplementaryFilter}
import htsjdk.samtools.reference.ReferenceSequenceFileWalker
import htsjdk.samtools.util.{IntervalList, SamLocusIterator, SequenceUtil}
import htsjdk.samtools.{SAMRecord, SAMSequenceDictionary}
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.mutable
import scala.util.Failure

@clp(group=ClpGroups.SamOrBam, description =
  """
    |Calculates the error rate by read position on coordinate sorted mapped BAMs. The output file contains
    |a row per read (first of pair, second of pair and unpaired), per position in read, with the total number
    |of bases observed, the number of errors observed, the overall error rate, and the rate of each kind of
    |substitution error.
    |
    |Substitution types are collapsed based on the reference or expected base, with only six substitution
    |types being reported: `A>C`, `A>G`, `A>T`, `C>A`, `C>G` and `C>T`.  For example, `T>G` is grouped in
    |with `A>C`.
    |
    |Analysis can be restricted to a set of intervals via the `--intervals` option. Genomic positions can be
    |excluded from analysis by supplying a set of variants (either known variants in the sample or a catalog
    |of known variants such as dbSNP).  For data believed to have low error rates it is recommended to use
    |both the `--intervals` and `--variants` options to restrict analysis to only regions expected to be
    |homozygous reference in the data.
    |
    |The following are reads / bases are excluded from the analysis:
    |
    |- Unmapped reads
    |- Reads marked as failing vendor quality
    |- Reads marked as duplicates (unless `--include-duplicates` is specified)
    |- Secondary and supplemental records
    |- Soft-clipped bases in records
    |- Reads with MAPQ < `--min-mapping-quality` (default: 20)
    |- Bases with base quality < `--min-base-quality` (default: 0)
    |- Bases where either the read base or the reference base is non-ACGT
    |
    |An output text file is generated with the extension `.error_rate_by_read_position.txt`
    |
    |If R's `Rscript` utility is on the path and `ggplot2` is installed in the R distribution then a PDF
    |of error rate plots will also be generated with extension `.error_rate_by_read_position.pdf`.
  """
  )
class ErrorRateByReadPosition
( @arg(flag='i', doc="Input BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output metrics prefix. If not given, will use the input BAM basename.") val output: Option[PathPrefix] = None,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: PathToFasta,
  @arg(flag='v', doc="Optional file of variant sites to ignore.") val variants: Option[PathToVcf] = None,
  @arg(flag='l', doc="Optional list of intervals to restrict analysis to.") val intervals: Option[PathToIntervals] = None,
  @arg(flag='d', doc="Include duplicate reads, otherwise ignore.") val includeDuplicates: Boolean = false,
  @arg(flag='m', doc="The minimum mapping quality for a read to be included.") val minMappingQuality: Int = 20,
  @arg(flag='q', doc="The minimum base quality for a base to be included.") val minBaseQuality: Int = 0
) extends FgBioTool with LazyLogging {

  private val ScriptPath = "com/fulcrumgenomics/bam/ErrorRateByReadPosition.R"

  Io.assertReadable(input)
  Io.assertReadable(ref)
  variants.foreach(Io.assertReadable)
  intervals.foreach(Io.assertReadable)
  output.foreach(out => Io.assertCanWriteFile(out))

  def execute(): Unit = {
    val in          = SamSource(input, ref=Some(ref))
    val metrics     = computeMetrics(in)
    val prefix      = output.getOrElse(PathUtil.removeExtension(input))
    val metricsPath = PathUtil.pathTo(prefix + ErrorRateByReadPositionMetric.FileExtension)
    val plotPath    = PathUtil.pathTo(prefix + ErrorRateByReadPositionMetric.PlotExtension)
    Metric.write(metricsPath, metrics=metrics)

    // And try plotting
    if (metrics.isEmpty) {
      logger.warning("No metrics generated (is your BAM empty?). Plots will not be generated.")
    }
    else {
      val description = plotDescription(in, input)
      Rscript.execIfAvailable(ScriptPath, metricsPath.toString, plotPath.toString, description) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => Unit
      }
    }
  }

  /** Computes the metrics from the BAM file. */
  private[bam] def computeMetrics(in: SamSource= SamSource(input, ref=Some(ref))): Seq[ErrorRateByReadPositionMetric] = {
    val progress = ProgressLogger(logger, verb="Processed", noun="loci", unit=50000)
    val ilist    = this.intervals.map(p => IntervalList.fromFile(p.toFile).uniqued(false))

    val refWalker     = new ReferenceSequenceFileWalker(this.ref.toFile)
    val locusIterator = buildSamLocusIterator(in, ilist).iterator()
    val variantMask   = buildVariantMask(variants, ilist, refWalker.getSequenceDictionary)

    val counters = List.tabulate(3)(_ => mutable.Map[Int, ObsCounter]())

    locusIterator.filterNot(l => variantMask.isVariant(l.getSequenceIndex, l.getPosition)).foreach { locus =>
      val refBase = SequenceUtil.upperCase(refWalker.get(locus.getSequenceIndex).getBases()(locus.getPosition - 1))
      if (SequenceUtil.isValidBase(refBase)) {
        locus.getRecordAndOffsets.iterator().filter(r => SequenceUtil.isValidBase(r.getReadBase)).foreach { rec =>

          val readBase = SequenceUtil.upperCase(rec.getReadBase)
          val readNum  = readNumber(rec.getRecord)
          val cycle    = if (rec.getRecord.getReadNegativeStrandFlag) rec.getRecord.getReadLength - rec.getOffset else rec.getOffset + 1
          val counter  = counters(readNum).getOrElseUpdate(cycle, new ObsCounter(readNum, cycle))

          if (refBase == 'A' || refBase == 'C') counter.count(refBase.toChar, readBase.toChar)
          else counter.count(SequenceUtil.complement(refBase).toChar, SequenceUtil.complement(readBase).toChar)
        }
      }

      progress.record(locus.getSequenceName, locus.getPosition)
    }
    in.safelyClose()

    // Ensure that all the maps have values up to the max
    for ((counter, readNum) <- counters.filter(_.nonEmpty).zipWithIndex; cycle <- 1 to counter.keys.max) {
      counter.getOrElseUpdate(cycle, new ObsCounter(readNum, cycle))
    }

    counters.flatMap(c => c.values.map(_.toMetric)).sortBy(m => (m.read_number, m.position))
  }

  /** Generates a SamLocusIterator that will traverse over the relevant parts of the BAM. */
  private def buildSamLocusIterator(in: SamSource, intervals: Option[IntervalList]): SamLocusIterator = {
    val iterator = intervals match {
      case None => new SamLocusIterator(in.toSamReader)
      case Some(i) => new SamLocusIterator(in.toSamReader, i)
    }

    val filters = new util.ArrayList[SamRecordFilter]()
    filters.add(new SecondaryOrSupplementaryFilter)
    filters.add(new FailsVendorReadQualityFilter)
    if (!this.includeDuplicates) filters.add(new DuplicateReadFilter)
    iterator.setSamFilters(filters)

    iterator.setMappingQualityScoreCutoff(this.minMappingQuality)
    iterator.setQualityScoreCutoff(this.minBaseQuality)
    iterator.setEmitUncoveredLoci(false)
    iterator.setIncludeIndels(false)
    iterator
  }

  /** Generates a variant mask object. */
  private def buildVariantMask(variants: Option[PathToVcf], intervals: Option[IntervalList], dict: SAMSequenceDictionary): VariantMask = {
    this.variants match {
      case None =>
        new VariantMask(Iterator.empty, dict)
      case Some(path) =>
        val reader = new VCFFileReader(path.toFile, this.intervals.isDefined)
        intervals match {
          case None => new VariantMask(reader.iterator(), dict)
          case Some(i) => new VariantMask(ByIntervalListVariantContextIterator(reader, i), dict)
        }
    }
  }

  /** Compute the read number for the record. */
  def readNumber(rec: SAMRecord): Int = if (!rec.getReadPairedFlag) 0 else if (rec.getFirstOfPairFlag) 1 else 2
}

private class ObsCounter(readNumber: Int, position: Int) {
  var a_ref_obs: Long  = 0
  var c_ref_obs: Long  = 0
  var a_to_c_obs: Long = 0
  var a_to_g_obs: Long = 0
  var a_to_t_obs: Long = 0
  var c_to_a_obs: Long = 0
  var c_to_g_obs: Long = 0
  var c_to_t_obs: Long = 0

  /** Increments the appropriate counter. */
  def count(ref: Char, read: Char): Unit = (ref, read) match {
    case ('A', 'A') => a_ref_obs  += 1
    case ('A', 'C') => a_to_c_obs += 1
    case ('A', 'G') => a_to_g_obs += 1
    case ('A', 'T') => a_to_t_obs += 1
    case ('C', 'A') => c_to_a_obs += 1
    case ('C', 'C') => c_ref_obs  += 1
    case ('C', 'G') => c_to_g_obs += 1
    case ('C', 'T') => c_to_t_obs += 1
    case _          => unreachable("Should never be invoked with a case other than the above.")
  }


  def toMetric: ErrorRateByReadPositionMetric = {
    val totalA = a_ref_obs + a_to_c_obs + a_to_g_obs + a_to_t_obs
    val totalC = c_ref_obs + c_to_a_obs + c_to_g_obs + c_to_t_obs
    val total  = totalA + totalC
    val errors = total - a_ref_obs - c_ref_obs

    new ErrorRateByReadPositionMetric(
      read_number = readNumber,
      position    = position,
      bases_total = total,
      errors      = errors,
      error_rate  = if (total == 0) 0 else errors / total.toDouble,
      a_to_c_error_rate = if (totalA == 0) 0 else a_to_c_obs / totalA.toDouble,
      a_to_g_error_rate = if (totalA == 0) 0 else a_to_g_obs / totalA.toDouble,
      a_to_t_error_rate = if (totalA == 0) 0 else a_to_t_obs / totalA.toDouble,
      c_to_a_error_rate = if (totalC == 0) 0 else c_to_a_obs / totalC.toDouble,
      c_to_g_error_rate = if (totalC == 0) 0 else c_to_g_obs / totalC.toDouble,
      c_to_t_error_rate = if (totalC == 0) 0 else c_to_t_obs / totalC.toDouble
    )
  }
}

object ErrorRateByReadPositionMetric {
  val FileExtension = ".error_rate_by_read_position.txt"
  val PlotExtension = ".error_rate_by_read_position.pdf"
}

/**
  * Metrics produced by `ErrorRateByReadPosition` describing the number of base observations and
  * substitution errors at each position within each sequencing read.  Error rates are given for
  * the overall substitution error rate and also for each kind of substitution separately. Instead
  * of reporting 12 substitution rates, 6 are reported where complementary substitutions are grouped
  * together, e.g. `T>G` substitutions are reported as `A>C`.
  *
  * @param read_number The read number (0 for fragments, 1 for first of pair, 2 for second of pair).
  * @param position The position or cycle within the read (1-based).
  * @param bases_total The total number of bases observed at this position.
  * @param errors The total number of errors or non-reference basecalls observed at this position.
  * @param error_rate The overall error rate at position.
  * @param a_to_c_error_rate The rate of `A>C` (and `T>G`) errors at the position.
  * @param a_to_g_error_rate The rate of `A>G` (and `T>C`) errors at the position.
  * @param a_to_t_error_rate The rate of `A>T` (and `T>A`) errors at the position.
  * @param c_to_a_error_rate The rate of `C>A` (and `G>T`) errors at the position.
  * @param c_to_g_error_rate The rate of `C>G` (and `G>C`) errors at the position.
  * @param c_to_t_error_rate The rate of `C>T` (and `G>A`) errors at the position.
  */
case class ErrorRateByReadPositionMetric
( read_number: Int,
  position: Int,
  bases_total: Count,
  errors: Count,
  error_rate: Double,
  a_to_c_error_rate: Double,
  a_to_g_error_rate: Double,
  a_to_t_error_rate: Double,
  c_to_a_error_rate: Double,
  c_to_g_error_rate: Double,
  c_to_t_error_rate: Double
) extends Metric
