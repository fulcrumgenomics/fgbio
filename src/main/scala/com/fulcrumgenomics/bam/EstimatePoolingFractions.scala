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
import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util.{Io, Metric, Sequences}
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import htsjdk.samtools.util._
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Examines sequence data generated from a pooled sample and estimates the fraction of sequence data
    |coming from each constituent sample. Uses a VCF of known genotypes for the samples within the
    |mixture along with a BAM of sequencing data derived from the pool.  Performs a multiple regression
    |for the alternative allele fractions at each SNP locus, using as inputs the individual sample's genotypes.
    |Only SNPs that are bi-allelic within the pooled samples are used.
    |
    |Each sample's contribution of REF vs. ALT alleles at each site is derived in one of two ways. If
    |the sample's genotype in the VCF has the `AF` attribute then the value from that field will be used.  If the
    |genotype has no AF attribute then the contribution is estimated based on the genotype (e.g. 0/0 will be 100%
    |ref, 0/1 will be 50% ref and 50% alt, etc.).
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
 @arg(flag='q', doc="Minimum base quality.") val minBaseQuality:Int = 5,
 @arg(doc="Examine input reads by sample given in each read's read group.") bySample: Boolean = false
) extends FgBioTool with LazyLogging {
  Io.assertReadable(vcf :: bam :: intervals.toList)

  private val Ci99Width = 2.58 // Width of a 99% confidence interval in units of std err

  private val AllReadGroupsName: String = "all"

  /* Class to hold information about a single locus. */
  case class Locus(chrom: String,
                   pos: Int,
                   ref: Char,
                   alt: Char,
                   expectedSampleFractions: Array[Double],
                   var observedFraction: Map[String, Double] = Map.empty)

  override def execute(): Unit = {
    val vcfReader   = VcfSource(vcf)
    val sampleNames = pickSamplesToUse(vcfReader)
    val intervals   = loadIntervals

    // Get the expected fractions from the VCF
    val vcfIterator = constructVcfIterator(vcfReader, intervals, sampleNames)
    val loci = vcfIterator.filterNot(v => this.nonAutosomes.contains(v.chrom)).map { v => Locus(
      chrom = v.chrom,
      pos   = v.pos,
      ref   = v.alleles.ref.bases.charAt(0),
      alt   = v.alleles.alts.head.value.charAt(0),
      expectedSampleFractions = sampleNames.map { s =>
        val gt = v.gt(s)
        gt.get[IndexedSeq[Float]]("AF") match {
          case None      => if (gt.isHomRef) 0 else if (gt.isHet) 0.5 else 1.0
          case Some(afs) => afs(0)
        }
      }
    )}.toArray

    logger.info(s"Loaded ${loci.length} bi-allelic SNPs from VCF.")

    fillObserveredFractionAndFilter(loci, this.minMeanSampleCoverage * sampleNames.length)

    val observedSamples = loci.flatMap { locus => locus.observedFraction.keySet }.distinct.sorted
    logger.info(f"Regressing on ${observedSamples.length}%,d input samples.")

    val regression = new OLSMultipleLinearRegression
    regression.setNoIntercept(true) // Intercept should be at 0!

    val metrics = observedSamples.flatMap { observedSample =>
      logger.info(f"Examining $observedSample")
      val (observedFractions, lociExpectedSampleFractions) = loci.flatMap { locus =>
        locus.observedFraction.get(observedSample).map { observedFraction =>
          (observedFraction, locus.expectedSampleFractions)
        }
      }.unzip
      logger.info(f"Regressing on ${observedFractions.length}%,d of ${loci.length}%,d loci that met coverage requirements.")
      regression.newSampleData(
        observedFractions,
        lociExpectedSampleFractions
      )

      val regressionParams = regression.estimateRegressionParameters()
      val total            = regressionParams.sum
      val fractions        = regressionParams.map(_ / total)
      val stderrs          = regression.estimateRegressionParametersStandardErrors().map(_ / total)
      logger.info(s"R^2 = ${regression.calculateRSquared()}")
      logger.info(s"Sum of regression parameters = ${total}")

      if (regression.estimateRegressionParameters().exists(_ < 0)) {
        logger.error("#################################################################################")
        logger.error("# One or more samples is estimated to have fraction < 0. This is likely due to  #")
        logger.error("# incorrect samples being used, insufficient coverage and/or too few SNPs.      #")
        logger.error("#################################################################################")
        fail(1)
      }

      sampleNames.toSeq.zipWithIndex.map { case (pool_sample, index) =>
        val sites      = lociExpectedSampleFractions.count(expectedSampleFractions => expectedSampleFractions(index) > 0)
        val singletons = lociExpectedSampleFractions.count { expectedSampleFractions =>
          expectedSampleFractions(index) > 0 && expectedSampleFractions.sum == expectedSampleFractions(index)
        }
        PoolingFractionMetric(
          observed_sample       = observedSample,
          pool_sample        = pool_sample,
          variant_sites      = sites,
          singletons         = singletons,
          estimated_fraction = fractions(index),
          standard_error     = stderrs(index),
          ci99_low           = max(0, fractions(index) - stderrs(index)*Ci99Width),
          ci99_high          = min(1, fractions(index) + stderrs(index)*Ci99Width))
      }
    }

    logger.info("Writing metrics")
    Metric.write(output, metrics)
  }

  /** Verify a provided sample list, or if none provided retrieve the set of samples from the VCF. */
  private def pickSamplesToUse(vcfIn: VcfSource): Array[String] = {
    if (this.samples.isEmpty) vcfIn.header.samples.toArray else {
      val samplesInVcf   = vcfIn.header.samples
      val missingSamples = samples.toSet.diff(samplesInVcf.toSet)
      if (missingSamples.nonEmpty) fail(s"Samples not present in VCF: ${missingSamples.mkString(", ")}")
      else samples.toArray.sorted
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
  def constructVcfIterator(in: VcfSource, intervals: Option[IntervalList], samples: Seq[String]): Iterator[Variant] = {
    val iterator: Iterator[Variant] = intervals match {
      case None     => in.iterator
      case Some(is) => is.flatMap(i => in.query(i.getContig, i.getStart, i.getEnd))
    }

    iterator
      .filter(v => v.filters.isEmpty || v.filters == Variant.PassingFilters)
      .filter(v => v.alleles.size == 2 && v.alleles.forall(a => a.value.length == 1))  // Just biallelic SNPs
      .filter(v => samples.map(v.gt).forall(gt => gt.isFullyCalled && (this.minGenotypeQuality <= 0 || gt.get[Int]("GQ").exists(_ >= this.minGenotypeQuality))))
      .map   (v => v.copy(genotypes=v.genotypes.filter { case (s, _) => samples.contains(s)} ))
      .filter(v => v.gts.flatMap(_.calls).toSet.size > 1)  // Not monomorphic
  }

  /** Constructs a SamLocusIterator that will visit every locus in the input. */
  def constructBamIterator(loci: Iterable[Locus]): Iterator[LocusInfo] = {
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

  /** Computes the observed fraction of the alternate allele at the given locus*/
  private def getObservedFraction(recordAndOffsets: Seq[SamLocusIterator.RecordAndOffset],
                                  locus: Locus,
                                  minCoverage: Int): Option[Double] = {
    if (recordAndOffsets.length < minCoverage) None else {
      val counts = BaseCounts(recordAndOffsets)
      val (ref, alt) = (counts(locus.ref), counts(locus.alt))

      // Somewhat redundant with check above, but this protects against a large fraction
      // of Ns or other alleles, and also against a large proportion of overlapping reads
      if (ref + alt < minCoverage) None else {
        Some(alt / (ref + alt).toDouble)
      }
    }
  }

  /**
    * Fills in the observedFraction field for each locus that meets coverage and then returns
    * the subset of loci that met coverage.
    */
  def fillObserveredFractionAndFilter(loci: Array[Locus], minCoverage: Int): Unit = {
    val locusIterator = constructBamIterator(loci)
    locusIterator.zip(loci.iterator).foreach { case (locusInfo, locus) =>
      if (locusInfo.getSequenceName != locus.chrom || locusInfo.getPosition != locus.pos) fail("VCF and BAM iterators out of sync.")

      if (bySample) {
        locus.observedFraction = locusInfo.getRecordAndOffsets.toSeq
          .groupBy(_.getRecord.getReadGroup.getSample)
          .flatMap { case (sample, recordAndOffsets) =>
            val observedFraction = getObservedFraction(
              recordAndOffsets = recordAndOffsets,
              locus            = locus,
              minCoverage      = minCoverage
            )
            observedFraction.map(frac => sample -> frac)
         }
      }
      else {
        val observedFraction = getObservedFraction(
          recordAndOffsets = locusInfo.getRecordAndOffsets.toSeq,
          locus            = locus,
          minCoverage      = minCoverage
        )
        observedFraction.foreach { frac =>
          locus.observedFraction = Map(AllReadGroupsName -> frac)
        }
      }
    }
  }
}

/**
  * Metrics produced by `EstimatePoolingFractions` to quantify the estimated proportion of a sample
  * mixture that is attributable to a specific sample with a known set of genotypes.
  *
  * @param observed_sample The name of the input sample as reported in the read group, or "all" if all read groups are
  *                     being treated as a single input sample.
  * @param pool_sample The name of the sample within the pool being reported on.
  * @param variant_sites How many sites were examined at which the reported sample is known to be variant.
  * @param singletons How many of the variant sites were sites at which only this sample was variant.
  * @param estimated_fraction The estimated fraction of the pool that comes from this sample.
  * @param standard_error The standard error of the estimated fraction.
  * @param ci99_low  The lower bound of the 99% confidence interval for the estimated fraction.
  * @param ci99_high The upper bound of the 99% confidence interval for the estimated fraction.
  */
case class PoolingFractionMetric(observed_sample: String,
                                 pool_sample: String,
                                 variant_sites: Count,
                                 singletons: Count,
                                 estimated_fraction: Proportion,
                                 standard_error: Double,
                                 ci99_low: Proportion,
                                 ci99_high: Proportion
                                ) extends Metric
