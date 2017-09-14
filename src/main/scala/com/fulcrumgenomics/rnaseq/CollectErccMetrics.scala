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
 *
 */

package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathPrefix, PathToBam, SafelyClosable, plotDescription}
import com.fulcrumgenomics.bam.api._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{DelimitedDataParser, LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger, Rscript}
import enumeratum.EnumEntry
import org.apache.commons.math3.exception.MathIllegalArgumentException
import org.apache.commons.math3.stat.correlation.{PearsonsCorrelation, SpearmansCorrelation}
import org.apache.commons.math3.stat.regression.SimpleRegression

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.util.Failure

sealed trait ErccMixture extends EnumEntry
object ErccMixture extends FgBioEnum[ErccMixture] {
  def values: IndexedSeq[ErccMixture] = findValues
  case object Mix1 extends ErccMixture
  case object Mix2 extends ErccMixture
}

object CollectErccMetrics {
  private[rnaseq] val StandardMixtureMetadataPath = "/com/fulcrumgenomics/rnaseq/ErccStandardMixtures.txt"
}

@clp(group = ClpGroups.RnaSeq, description=
  """
    |Collects metrics for ERCC spike-ins for RNA-Seq experiments.
    |
    |Currently calculates per-transcript ERCC metrics and summarizes dose response, but does not calculate fold-change
    |response.
    |
    |The input BAM should contain reads mapped to a reference containing the ERCC transcripts.  The reference may have
    |additional contigs, for example, when concatenating the sample's reference genome and the ERCC transcripts.  The
    |BAM should have sequence lines in the header matching the ERCC ids (ex. ERCC-00130 or ERCC-00004).
    |
    |The standard ERCC transcripts, including their unique IDs and concentrations, are taken from
    |[here](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).  The second column lists the ERCC
    |transcript identifier which should be present as a sequence line in the BAM header, while columns four and five
    |give the concentration of the transcript for mixtures #1 and #2.
    |
    |The choice of mixture to use can be specified with the `--mixture-name` option (either `Mix1` or `Mix2`), or with
    |a file containing a custom list of transcripts and concentrations using the `--custom-mixture` option as follows.
    |The custom mixture file should be tab-delimited file containing the following columns:
    |  1. ERCC ID - each ERCC ID should match a contig/reference-sequence name in the input SAM/BAM header.
    |  2. Concentration - the concentration (in `attomoles/ul`).
    |The custom mixture file should contain a header line with names for each column, though the actual values will be ignored.
    |
    |Three outputs will be produced:
    |  1. <output>.ercc_summary_metrics.txt - summary statistics for total # of reads mapping to the ERCC transcripts and dose
    |                                response metrics.
    |  2. <output>.ercc_detailed_metrics.txt - gives a per-ERCC-transcript expected concentration and observed fragment count.
    |  3. <output>.ercc_metrics.pdf - plots the expected concentration versus the observed fragment count.
    |
    |Secondary andsupplementary reads will be ignored.  A read pair mapping to an ERCC transcript is counted only if both
    |ends of the pair map to the same ERCC transcript.  A minimum mapping quality can be required for reads to be
    |counted as mapped to an ERCC transcript.
  """)
class CollectErccMetrics
(@arg(flag='i', doc="Input SAM or BAM file of aligned reads in coordinate order.") val input: PathToBam,
 @arg(flag='o', doc="Output prefix.") val output: PathPrefix,
 @arg(flag='m', doc="The name of the standard ERCC mixture.", mutex=Array("customMixture")) val mixtureName: Option[ErccMixture] = None,
 @arg(          doc="Tab-delimited file containing ERCC IDs and expected concentrations.", mutex=Array("mixtureName")) val customMixture: Option[FilePath] = None,
 @arg(flag='c', doc="Minimum # of counts required to include an ERCC transcript.") val minTranscriptCount: Int = 3,
 @arg(flag='M', doc="The minimum mapping quality") val minimumMappingQuality: Int = 10
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  customMixture.foreach(Io.assertReadable)
  Io.assertCanWriteFile(output)

  if (mixtureName.isEmpty && customMixture.isEmpty) {
    throw new ValidationException("Either --mixture-name or --mixture-metadata must be provided")
  }

  private val ScriptPath = "com/fulcrumgenomics/rnaseq/CollectErccMetrics.R"

  override def execute(): Unit = {
    val in       = SamSource(input)
    val progress = ProgressLogger(logger, unit=2e6.toInt)
    val erccData = customMixture match {
      case Some(path) => DelimitedDataParser(path, '\t').map { row => row[String](0) -> row[Double](1) }.toMap
      case None       =>
        val stream    = getClass.getResourceAsStream(CollectErccMetrics.StandardMixtureMetadataPath)
        val lines     = Source.fromInputStream(stream).withClose(() => stream.close()).getLines.filterNot(_.startsWith("#"))
        val mixColumn = if (mixtureName.get == ErccMixture.Mix1) 3 else 4
        new DelimitedDataParser(lines, '\t').map { row => row[String](1) -> row[Double](mixColumn) }.toMap
    }

    // check for no ERCC data
    if (erccData.isEmpty) {
      val description = customMixture.map("in '" + _).getOrElse(mixtureName.map("for '" + _).get) + "'"
      throw new IllegalArgumentException(s"No ERCC data given " + description)
    }

    // require that the SAM/BAM header contains all the ERCC transcripts
    val missingErccTranscripts = erccData.keySet -- in.header.getSequenceDictionary.getSequences.map(_.getSequenceName)
    if (missingErccTranscripts.nonEmpty) {
      throw new IllegalArgumentException(s"Missing ERCC transcripts in the SAM/BAM header: ${missingErccTranscripts.mkString(", ")}")
    }

    var totalReads: Long     = 0 // the total # of reads in the BAM file
    var erccReads: Long      = 0 // the total # of reads that are counted as mapped to an ERCC transcript
    var erccTemplates: Long  = 0 // the total # of read pairs (or single end reads) that are counted as mapped to an ERCC transcript
    val erccTemplatesCounter = new SimpleCounter[String]() // the count of read pairs (or single end reads) mapped to a given ERCC transcript
    in.iterator
      .filterNot { r => progress.record(r); r.secondary || r.supplementary } // no secondary or supplementary
      .filter { r =>
        totalReads += 1
        // consider only records mapped to an ERCC transcript
        r.mapped && r.mapq >= minimumMappingQuality && erccData.contains(r.refName)
      }
      .filter { r =>
        erccReads += 1
        // consider only single-end records, or paired records where both ends map to the sam ERCC transcript
        !r.paired || (r.mateMapped && r.refIndex == r.mateRefIndex && r.firstOfPair)

      }
      .foreach { r =>
        erccTemplates += 1
        erccTemplatesCounter.count(r.refName)
      }
    in.safelyClose()

    // Transform to log2 counts to compute some simple correlation and linear regression metrics.
    val counts            = ArrayBuffer[Double]()
    val normalizedCounts  = ArrayBuffer[Double]()
    val concentrations    = ArrayBuffer[Double]()
    val maximumErccLength = erccData.keys.map { name => in.header.getSequence(name).getSequenceLength }.max
    erccData.iterator.foreach { case (name, concentration) =>
      val count = erccTemplatesCounter.countOf(name)
      if (minTranscriptCount <= count) {
        val erccLength = in.header.getSequence(name).getSequenceLength
        counts += log2(count)
        normalizedCounts += log2(count * maximumErccLength / erccLength.toDouble) // normalize by the ERCC transcript length
        concentrations += log2(concentration)
      }
    }

    // Rocket science
    val simpleRegression = new SimpleRegression()
    normalizedCounts.zip(concentrations).foreach { case (cnt, conc) => simpleRegression.addData(cnt, conc) }

    val erccMetrics = ErccSummaryMetrics(
      total_reads                = totalReads,
      ercc_reads                 = erccReads,
      fraction_ercc_reads        = if (totalReads == 0) 0 else erccReads / totalReads.toDouble,
      ercc_templates             = erccTemplates,
      total_transcripts          = erccTemplatesCounter.iterator.length,
      passing_filter_transcripts = erccTemplatesCounter.iterator.count(_._2 >= minTranscriptCount),
      pearsons_correlation       = noneIfNaNOrInsufficientData { new PearsonsCorrelation().correlation(normalizedCounts.toArray, concentrations.toArray) },
      spearmans_correlation      = noneIfNaNOrInsufficientData { new SpearmansCorrelation().correlation(normalizedCounts.toArray, concentrations.toArray) },
      intercept                  = noneIfNaN { simpleRegression.getIntercept },
      slope                      = noneIfNaN { simpleRegression.getSlope },
      r_squared                   = noneIfNaN { simpleRegression.getRSquare }
    )

    // Write the output
    {
      def f(ext: String): FilePath = PathUtil.pathTo(output + ext)

      val summaryPath  = f(".ercc_summary_metrics.txt")
      val detailedPath = f(".ercc_detailed_metrics.txt")
      val plotPath     = f(".ercc_plot.pdf")

      // summary metrics file
      Metric.write[ErccSummaryMetrics](summaryPath, erccMetrics)

      // detailed (per-ERCC-transcript) file, written in the order defined in the Sam header
      val metrics = in.header.getSequenceDictionary
        .getSequences
        .map(_.getSequenceName)
        .filter { name => erccData.contains(name) }
        .map { name =>
          val count           = erccTemplatesCounter.countOf(name)
          val normalizedCount = count / maximumErccLength.toDouble
          val concentration   = erccData(name)
          ErccDetailedMetric(name=name, concentration=concentration, count=count, normalized_count=normalizedCount)
        }.toSeq
      Metric.write[ErccDetailedMetric](detailedPath, metrics)

      // And try plotting
      if (metrics.isEmpty || totalReads == 0 || simpleRegression.getN <= 1) {
        logger.warning("No metrics generated (is your BAM empty?). Plots will not be generated.")
      }
      else {
        Rscript.execIfAvailable(ScriptPath, detailedPath.toString, plotPath.toString, plotDescription(in, input), minTranscriptCount.toString) match {
          case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
          case _ => Unit
        }
      }
    }
  }

  /** The natural log of 2! */
  private val lnOf2 = scala.math.log(2) // natural log of 2

  /** Computes the log base 2. */
  private def log2(x: Double): Double = scala.math.log(x) / lnOf2

  /** Use this method to wrap the calculation of correlations such that if insufficient data are given, None is returned. */
  private def noneIfNaNOrInsufficientData(f: => Double): Option[Double] = {
    try {
      noneIfNaN(f)
    } catch {
      case _: MathIllegalArgumentException => None
    }
  }

  /** Converts all `NaN` values to None. */
  private def noneIfNaN(d: Double): Option[Double] = if (d.isNaN) None else Some(d)
}

/**
  * Metrics produced by `CollectErccMetrics` describing various per-transcript metrics related to the spike-in of ERCC
  * (External RNA Controls Consortium) into an RNA-Seq experiment.  One metric per ERCC transcript will be present.
  *
  * @param name the name (or ID) of the ERCC transcript.
  * @param concentration the expected concentration as input to `CollectErccMetrics`.
  * @param count the observed count of the number of read pairs (or single end reads) .
  * @param normalized_count the observed count of the number of read pairs (or single end reads) normalized by the ERCC transcript length.
  */
case class ErccDetailedMetric
(
  name: String,
  concentration: Double,
  count: Long,
  normalized_count: Double
) extends Metric

/**
  * Metrics produced by `CollectErccMetrics` describing various summary metrics related to the spike-in of ERCC
  * (External RNA Controls Consortium) into an RNA-Seq experiment.
  *
  * The correlation coefficients and linear regression are calculated based on the log2 observed read pair count normalized
  * by ERCC transcript length versus the log2 expected concentration.
  *
  * @param total_reads the total number of reads considered.
  * @param ercc_reads the total number of reads mapping to an ERCC transcript.
  * @param fraction_ercc_reads the fraction of total reads that map to an ERCC transcript.
  * @param ercc_templates the total number of read pairs (or single end reads) mapping to an ERCC transcript.
  * @param total_transcripts the total number of ERCC transcripts with at least one read observed.
  * @param passing_filter_transcripts the total number of ERCC transcripts with at least the user-set minimum # of reads observed.
  * @param pearsons_correlation Pearson's correlation coefficient for correlation of concentration and normalized counts.
  * @param spearmans_correlation Spearman's correlation coefficient for correlation of concentration and normalized counts.
  * @param intercept the intercept of the linear regression.
  * @param slope the slope of the linear regression.
  * @param r_squared the r-squared of the linear regression.
  */
case class ErccSummaryMetrics
(
  total_reads: Long,
  ercc_reads: Long,
  fraction_ercc_reads: Double,
  ercc_templates: Long,
  total_transcripts: Int,
  passing_filter_transcripts: Int,
  pearsons_correlation: Option[Double],
  spearmans_correlation: Option[Double],
  intercept: Option[Double],
  slope: Option[Double],
  r_squared: Option[Double]
) extends Metric