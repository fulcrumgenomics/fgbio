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
 */

package com.fulcrumgenomics.bam

import java.lang.Math.{max, min}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util.{Io, Metric, Sequences}
import com.fulcrumgenomics.vcf.ByIntervalListVariantContextIterator
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import htsjdk.samtools.util._
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Examines sequence data generated from a pooled sample and estimates the fraction of sequence data
    |coming from each constituent sample. Uses a VCF of known genotypes for the samples within the
    |mixture along with a BAM of sequencing data derived from the pool.  Performs a multiple regression
    |for the alternative allele fractions at each SNP locus, using as inputs the individual sample's genotypes.
    |Only SNPs that are bi-allelic within the pooled samples are used.
    |
    |Various filtering parameters can be used to control which loci are used:
    |
    |- _--intervals_ will restrict analysis to variants within the described intervals
    |- _--min-genotype-quality_ will filter out any site with any genotype with GQ < n
    |- _--min-mean-sample-coverage_ requires that the coverage of a site in the BAM be >= `min-mean-sample-coverage * n_samples`
    |- _--min-mapping-quality_ filters out reads in the BAM with MQ < n
    |- _--min-base-quality_ filters out bases in the BAM with Q < n
  """)
class EstimatePoolingFractions
(@arg(flag='v', doc="VCF of individual sample genotypes.")  val vcf: PathToVcf,
 @arg(flag='b', doc="Path to BAM file of sequencing data.") val bam: PathToBam,
 @arg(flag='o', doc="Output file to write with pooling fractions.") val output: FilePath,
 @arg(flag='l', minElements=0, doc="Zero or more set of regions to restrict analysis to.") val intervals: Seq[PathToIntervals] = Seq.empty,
 @arg(flag='s', minElements=0, doc="Optional subset of samples from VCF to use.") val samples: Seq[String] = Seq.empty,
 @arg(flag='n', doc="Non-autosomal chromosomes to avoid.") val nonAutosomes: Seq[String] = Sequences.CommonNonAutosomalContigNames,
 @arg(flag='g', doc="Minimum genotype quality. Use -1 to disable.") val minGenotypeQuality: Int = 30,
 @arg(flag='c', doc="Minimum (sequencing coverage @ SNP site / n_samples).") val minMeanSampleCoverage: Int = 6,
 @arg(flag='m', doc="Minimum mapping quality.") val minMappingQuality: Int = 20,
 @arg(flag='q', doc="Minimum base quality.") val minBaseQuality:Int = 5
) extends FgBioTool with LazyLogging {
  Io.assertReadable(vcf :: bam :: intervals.toList)

  private val Ci99Width = 2.58 // Width of a 99% confidence interval in units of std err

  /* Class to hold information about a single locus. */
  case class Locus(chrom: String, pos: Int, ref: Char, alt: Char, expectedSampleFractions: Array[Double], var observedFraction: Option[Double] = None)

  override def execute(): Unit = {
    val vcfReader   = new VCFFileReader(vcf.toFile)
    val sampleNames = pickSamplesToUse(vcfReader)
    val intervals   = loadIntervals

    // Get the expected fractions from the VCF
    val vcfIterator = constructVcfIterator(vcfReader, intervals, sampleNames)
    val loci = vcfIterator.filterNot(v => this.nonAutosomes.contains(v.getContig)).map { v => Locus(
      chrom = v.getContig,
      pos = v.getStart,
      ref = v.getReference.getBaseString.charAt(0),
      alt = v.getAlternateAllele(0).getBaseString.charAt(0),
      expectedSampleFractions = sampleNames.map { s => val gt = v.getGenotype(s); if (gt.isHomRef) 0 else if (gt.isHet) 0.5 else 1.0 }
    )}.toArray

    logger.info(s"Loaded ${loci.length} bi-allelic SNPs from VCF.")

    val coveredLoci = fillObserveredFractionAndFilter(loci, this.minMeanSampleCoverage * sampleNames.length)

    logger.info(s"Regressing on ${coveredLoci.length} of ${loci.length} that met coverage requirements.")
    val regression = new OLSMultipleLinearRegression
    regression.setNoIntercept(true) // Intercept should be at 0!
    regression.newSampleData(
      coveredLoci.map(_.observedFraction.getOrElse(unreachable("observed fraction must be defined"))),
      coveredLoci.map(_.expectedSampleFractions)
    )

    val regressionParams = regression.estimateRegressionParameters()
    val total            = regressionParams.sum
    val fractions        = regressionParams.map(_ / total)
    val stderrs          = regression.estimateRegressionParametersStandardErrors().map(_ / total)
    logger.info(s"R^2 = ${regression.calculateRSquared()}")
    logger.info(s"Sum of regression parameters = ${total}")

    val metrics = sampleNames.toSeq.zipWithIndex.map { case (sample, index) =>
      val sites      = coveredLoci.count(l => l.expectedSampleFractions(index) > 0)
      val singletons = coveredLoci.count(l => l.expectedSampleFractions(index) > 0 && l.expectedSampleFractions.sum == l.expectedSampleFractions(index))
      new PoolingFractionMetric(
        sample             = sample,
        variant_sites      = sites,
        singletons         = singletons,
        estimated_fraction = fractions(index),
        standard_error     = stderrs(index),
        ci99_low           = max(0, fractions(index) - stderrs(index)*Ci99Width),
        ci99_high          = min(1, fractions(index) + stderrs(index)*Ci99Width))
    }

    Metric.write(output, metrics)

    if (regression.estimateRegressionParameters().exists(_ < 0)) {
      logger.error("#################################################################################")
      logger.error("# One or more samples is estimated to have fraction < 0. This is likely due to  #")
      logger.error("# incorrect samples being used, insufficient coverage and/or too few SNPs.      #")
      logger.error("#################################################################################")
      fail(1)
    }
  }

  /** Verify a provided sample list, or if none provided retrieve the set of samples from the VCF. */
  private def pickSamplesToUse(vcfReader: VCFFileReader): Array[String] = {
    if (samples.nonEmpty) {
      val samplesInVcf = vcfReader.getFileHeader.getSampleNamesInOrder.toIterator.toSet
      val missingSamples = samples.filterNot(samplesInVcf.contains)
      if (missingSamples.nonEmpty) fail(s"Samples not present in VCF: ${missingSamples.mkString(", ")}")
      else samples.toArray.sorted
    }
    else {
      vcfReader.getFileHeader.getSampleNamesInOrder.toIterator.toSeq.toArray.sorted // toSeq.toArray is necessary cos util.ArrayList.toArray() exists
    }
  }

  /** Loads up and merges all the interval lists provided. Returns None if no intervals were specified. */
  private def loadIntervals: Option[IntervalList] = this.intervals match {
    case head +: tail =>
      val list = IntervalList.fromFile(head.toFile)
      val dict = list.getHeader.getSequenceDictionary
      tail.foreach { f =>
        val tmp = IntervalList.fromFile(f.toFile)
        if (!SequenceUtil.areSequenceDictionariesEqual(dict, tmp.getHeader.getSequenceDictionary))
          fail(s"Sequence dictionaries differ between ${this.intervals.head} and ${f}")
        else
          list.addall(tmp.getIntervals)
      }
      Some(list.uniqued(false))
    case _ => None
  }

  /** Generates an iterator over non-filtered bi-allelic SNPs where all the required samples are genotyped. */
  def constructVcfIterator(in: VCFFileReader, intervals: Option[IntervalList], samples: Array[String]): Iterator[VariantContext] = {
    val vcfIterator: Iterator[VariantContext] = intervals match {
      case None     => in.toIterator
      case Some(is) => ByIntervalListVariantContextIterator(in, is)
    }

    val samplesAsUtilSet = CollectionUtil.makeSet(samples:_*)

    vcfIterator
      .filterNot(_.isFiltered)
      .map(_.subContextFromSamples(samplesAsUtilSet, true))
      .filter(v => v.isSNP && v.isBiallelic && !v.isMonomorphicInSamples)
      .filter(_.getNoCallCount == 0)
      .filter(v => v.getGenotypesOrderedByName.toIterator.forall(gt => gt.getGQ >= this.minGenotypeQuality))
  }

  /** Constructs a SamLocusIterator that will visit every locus in the input. */
  def constructBamIterator(loci: Traversable[Locus]): Iterator[LocusInfo] = {
    val in = SamSource(this.bam)
    val intervals = new IntervalList(in.header)
    loci.foreach(l => intervals.add(new Interval(l.chrom, l.pos, l.pos)))
    val iterator = new SamLocusIterator(in.toSamReader, intervals.uniqued())
    iterator.setEmitUncoveredLoci(true)
    iterator.setIncludeNonPfReads(false)
    iterator.setMappingQualityScoreCutoff(this.minMappingQuality)
    iterator.setQualityScoreCutoff(this.minBaseQuality)
    javaIteratorAsScalaIterator(iterator)
  }

  /**
    * Fills in the observedFraction field for each locus that meets coverage and then returns
    * the subset of loci that met coverage.
    */
  def fillObserveredFractionAndFilter(loci: Array[Locus], minCoverage: Int): Array[Locus] = {
    val locusIterator = constructBamIterator(loci)
    locusIterator.zip(loci.iterator).foreach { case (locusInfo, locus) =>
      if (locusInfo.getSequenceName != locus.chrom || locusInfo.getPosition != locus.pos) fail("VCF and BAM iterators out of sync.")

      // A gross coverage check here to avoid a lot of work; better check below
      if (locusInfo.getRecordAndOffsets.size() > minCoverage) {
        val counts = BaseCounts(locusInfo)
        val (ref, alt) = (counts(locus.ref), counts(locus.alt))

        // Somewhat redundant with check above, but this protects against a large fraction
        // of Ns or other alleles, and also against a large proportion of overlapping reads
        if (ref + alt >= minCoverage) {
          locus.observedFraction = Some(alt / (ref + alt).toDouble)
        }
      }
    }

    loci.filter(_.observedFraction.isDefined)
  }
}

/**
  * Metrics produced by `EstimatePoolingFractions` to quantify the estimated proportion of a sample
  * mixture that is attributable to a specific sample with a known set of genotypes.
  *
  * @param sample The name of the sample within the pool being reported on.
  * @param variant_sites How many sites were examined at which the reported sample is known to be variant.
  * @param singletons How many of the variant sites were sites at which only this sample was variant.
  * @param estimated_fraction The estimated fraction of the pool that comes from this sample.
  * @param standard_error The standard error of the estimated fraction.
  * @param ci99_low  The lower bound of the 99% confidence interval for the estimated fraction.
  * @param ci99_high The upper bound of the 99% confidence interval for the estimated fraction.
  */
case class PoolingFractionMetric(sample: String,
                                 variant_sites: Count,
                                 singletons: Count,
                                 estimated_fraction: Proportion,
                                 standard_error: Double,
                                 ci99_low: Proportion,
                                 ci99_high: Proportion
                                ) extends Metric
