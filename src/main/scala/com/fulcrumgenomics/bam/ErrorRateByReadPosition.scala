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
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Metric.Count
import com.fulcrumgenomics.util.{Metric, ProgressLogger, Rscript, Sequences}
import com.fulcrumgenomics.vcf.{ByIntervalListVariantContextIterator, VariantMask}
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.filter.{DuplicateReadFilter, FailsVendorReadQualityFilter, SamRecordFilter, SecondaryOrSupplementaryFilter}
import htsjdk.samtools.reference.ReferenceSequenceFileWalker
import htsjdk.samtools.util.{IntervalList, SamLocusIterator, SequenceUtil}
import htsjdk.variant.vcf.VCFFileReader

import scala.annotation.switch
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
  @arg(flag='q', doc="The minimum base quality for a base to be included.") val minBaseQuality: Int = 0,
  @arg(doc="Collapse substitution types based on the reference or expected base, with only six substitution" +
    " types being reported: `A>C`, `A>G`, `A>T`, `C>A`, `C>G` and `C>T`.For example, `T>G` is grouped in with `A>C`." +
    " Otherwise, all possible substitution types will be reported."
  ) val collapse: Boolean = true
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
    val metricsPath = PathUtil.pathTo(s"${prefix}${ErrorRateByReadPositionMetric.FileExtension}")
    val plotPath    = PathUtil.pathTo(s"${prefix}${ErrorRateByReadPositionMetric.PlotExtension}")
    Metric.write(metricsPath, metrics=metrics)

    // And try plotting
    if (metrics.isEmpty) {
      logger.warning("No metrics generated (is your BAM empty?). Plots will not be generated.")
    }
    else {
      val description = plotDescription(in, input)
      Rscript.execIfAvailable(ScriptPath, metricsPath.toString, plotPath.toString, description) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => ()
      }
    }
  }

  /** Computes the metrics from the BAM file. */
  private[bam] def computeMetrics(in: SamSource= SamSource(input, ref=Some(ref))): Seq[ErrorRateByReadPositionMetric] = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    val progress = ProgressLogger(logger, verb="Processed", noun="loci", unit=50000)
    val ilist    = this.intervals.map(p => IntervalList.fromFile(p.toFile).uniqued(false))

    val refWalker     = new ReferenceSequenceFileWalker(this.ref.toFile)
    val locusIterator = buildSamLocusIterator(in, ilist).iterator()
    val variantMask   = buildVariantMask(variants, ilist, refWalker.getSequenceDictionary.fromSam)

    val counters = List.tabulate(3)(_ => mutable.Map[Int, ObsCounter]())

    locusIterator.filterNot(l => variantMask.isVariant(l.getSequenceIndex, l.getPosition)).foreach { locus =>
      val refBase = SequenceUtil.upperCase(refWalker.get(locus.getSequenceIndex).getBases()(locus.getPosition - 1))
      if (SequenceUtil.isValidBase(refBase)) {
        locus.getRecordAndOffsets.iterator().filter(r => SequenceUtil.isValidBase(r.getReadBase)).foreach { rec =>

          val readBase = SequenceUtil.upperCase(rec.getReadBase)
          val readNum  = readNumber(rec.getRecord)
          val cycle    = if (rec.getRecord.getReadNegativeStrandFlag) rec.getRecord.getReadLength - rec.getOffset else rec.getOffset + 1
          val counter  = counters(readNum).getOrElseUpdate(cycle, new ObsCounter(readNum, cycle, collapse))

          if (!collapse || refBase == 'A' || refBase == 'C') counter.count(refBase.toChar, readBase.toChar)
          else counter.count(SequenceUtil.complement(refBase).toChar, SequenceUtil.complement(readBase).toChar)
        }
      }

      progress.record(locus.getSequenceName, locus.getPosition)
    }
    in.safelyClose()

    // Ensure that all the maps have values up to the max
    for ((counter, readNum) <- counters.filter(_.nonEmpty).zipWithIndex; cycle <- 1 to counter.keys.max) {
      counter.getOrElseUpdate(cycle, new ObsCounter(readNum, cycle, collapse))
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
  private def buildVariantMask(variants: Option[PathToVcf], intervals: Option[IntervalList], dict: SequenceDictionary): VariantMask = {
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

/** Counter for substitutions for a given read number and position in that read.
  *
  * @param readNumber the read number
  * @param position the position in the read
  * @param collapse true to collapse substitution types, false otherwise
  */
private class ObsCounter(readNumber: Int, position: Int, collapse: Boolean) {
  private val Dna: Seq[Char]             = IndexedSeq('A', 'C', 'G', 'T')
  private val counts: Array[Array[Long]] = new Array[Long](4).map(_ => new Array[Long](4))

  /** Maps A/C/G/T to an index 0-3 and throws an exception for any other base. */
  @inline private def index(base: Char): Int = (base: @switch) match {
    case 'A' => 0
    case 'C' => 1
    case 'G' => 2
    case 'T' => 3
    case _   => throw new IllegalArgumentException(s"Invalid base: $base")
  }

  @inline def count(ref: Char, read: Char): Unit =
    try { this.counts(index(ref))(index(read)) += 1 } catch { case _: IllegalArgumentException => () }

  @inline def countOf(ref: Char, read: Char): Long =
    try { this.counts(index(ref))(index(read)) } catch { case _: IllegalArgumentException => 0 }

  @inline private def total(ref: Char): Long    = this.counts(index(ref)).sum
  @inline private def totalRef(ref: Char): Long = countOf(ref, ref)

  @inline private def errorRate(ref: Char, read: Char): Double = {
    val totalRef = this.total(ref)
    if (totalRef == 0) 0 else countOf(ref=ref, read=read) / totalRef.toDouble
  }

  @inline private def getErrorRate(ref: Char, read: Char): Option[Double] = {
    if (!collapse || ref == 'A' || ref == 'C') Some(errorRate(ref=ref, read=read))
    else None
  }

  def toMetric: ErrorRateByReadPositionMetric = {
    // Double check that we have zero counts for G>X and T>X when we collapse substitution types
    if (collapse) {
      Seq('G', 'T').foreach { ref =>
        Dna.foreach { read =>
          require(this.countOf(ref=ref, read=read) == 0,
            s"Bug: collapse is true but found a non-zero count for $ref>$read")
        }
      }
    }

    val total  = Dna.map(ref => this.total(ref=ref)).sum
    val errors = total - Dna.map(ref => this.totalRef(ref=ref)).sum
    new ErrorRateByReadPositionMetric(
      read_number = readNumber,
      position    = position,
      bases_total = total,
      errors      = errors,
      error_rate  = if (total == 0) 0 else errors / total.toDouble,
      a_to_c_error_rate = errorRate('A', 'C'),
      a_to_g_error_rate = errorRate('A', 'G'),
      a_to_t_error_rate = errorRate('A', 'T'),
      c_to_a_error_rate = errorRate('C', 'A'),
      c_to_g_error_rate = errorRate('C', 'G'),
      c_to_t_error_rate = errorRate('C', 'T'),
      g_to_a_error_rate = getErrorRate('G', 'A'),
      g_to_c_error_rate = getErrorRate('G', 'C'),
      g_to_t_error_rate = getErrorRate('G', 'T'),
      t_to_a_error_rate = getErrorRate('T', 'A'),
      t_to_c_error_rate = getErrorRate('T', 'C'),
      t_to_g_error_rate = getErrorRate('T', 'G'),
      collapsed         = collapse
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
  * the overall substitution error rate and also for each kind of substitution separately.
  *
  * If `collapsed` is `true`, then complementary substitutions are grouped together into the first 6 error rates.
  * e.g. `T>G` substitutions are reported as `A>C`.  Otherwise, all 12 substitution rates are reported.
  *
  * @param read_number The read number (0 for fragments, 1 for first of pair, 2 for second of pair).
  * @param position The position or cycle within the read (1-based).
  * @param bases_total The total number of bases observed at this position.
  * @param errors The total number of errors or non-reference basecalls observed at this position.
  * @param error_rate The overall error rate at position.
  * @param a_to_c_error_rate The rate of `A>C` (and `T>G` when collapsed) errors at the position.
  * @param a_to_g_error_rate The rate of `A>G` (and `T>C` when collapsed) errors at the position.
  * @param a_to_t_error_rate The rate of `A>T` (and `T>A` when collapsed) errors at the position.
  * @param c_to_a_error_rate The rate of `C>A` (and `G>T` when collapsed) errors at the position.
  * @param c_to_g_error_rate The rate of `C>G` (and `G>C` when collapsed) errors at the position.
  * @param c_to_t_error_rate The rate of `C>T` (and `G>A` when collapsed) errors at the position.
  * @param g_to_a_error_rate The rate of `G>A` errors at the position.
  * @param g_to_c_error_rate The rate of `G>C` errors at the position.
  * @param g_to_t_error_rate The rate of `G>T` errors at the position.
  * @param t_to_a_error_rate The rate of `T>A` errors at the position.
  * @param t_to_c_error_rate The rate of `T>C` errors at the position.
  * @param t_to_g_error_rate The rate of `T>T` errors at the position.
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
  c_to_t_error_rate: Double,
  g_to_a_error_rate: Option[Double] = None,
  g_to_c_error_rate: Option[Double] = None,
  g_to_t_error_rate: Option[Double] = None,
  t_to_a_error_rate: Option[Double] = None,
  t_to_c_error_rate: Option[Double] = None,
  t_to_g_error_rate: Option[Double] = None,
  collapsed: Boolean = true
) extends Metric
