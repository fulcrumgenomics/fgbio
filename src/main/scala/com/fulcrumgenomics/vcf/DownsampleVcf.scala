package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api.{Allele, Genotype, Variant, VcfCount, VcfFieldType, VcfFormatHeader, VcfHeader, VcfSource, VcfWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.vcf.DownsampleVcf.{downsampleAndRegenotype, winnowVariants}

import scala.math.log10
import scala.util.Random
import scala.tools.nsc.doc.html.HtmlTags

object DownsampleVcf extends LazyLogging {
  /** Removes variants that are within a specified distance from a previous variant.
   * The end position of the current variant is compared with the start position of the following variant.
   * @param variants   an iterator of the variants to process
   * @param windowSize the interval (exclusive) in which to check for additional variants.
   *                   windowSize considers the distance between the end position of a variant
   *                   with the start position of the following variant
   * @param dict       a sequencing dictionary to get contig ordering
   * @return a new iterator of variants with just the variant entries we want to keep
   */
  def winnowVariants(variants: Iterator[Variant], windowSize: Int, dict: SequenceDictionary): Iterator[Variant] = {
    require(windowSize >= 0, f"the windowSize ($windowSize) is negative")
    new Iterator[Variant] {
      private val iter = variants.bufferBetter

      def hasNext: Boolean = iter.hasNext

      private def isInOrder(current: Variant, next: Variant, currentIndex: Int, nextIndex: Int): Boolean = {
        (currentIndex < nextIndex) || (currentIndex == nextIndex && current.end <= next.pos)
      }

      def next(): Variant = {
        val current = iter.next()
        val currentIndex = dict(current.chrom).index
        iter.dropWhile { next: Variant =>
          val nextIndex = dict(next.chrom).index
          require(
            isInOrder(current, next, currentIndex, nextIndex),
            f"variants out of order; ${current.chrom}:${current.pos} > ${next.chrom}:${next.pos}")

          currentIndex == nextIndex && next.pos - current.end < windowSize
        }
        current
      }
    }
  }

  /** Downsamples variants by randomly sampling the total allele depths at the given proportion.
   * @param oldAds     an indexed seq of the original allele depths
   * @param proportion the proportion to use for downsampling,
   *                   calculated using total base count from the index and a target base count
   * @return a new IndexedSeq of allele depths of the same length as `oldAds`
   */
  def downsampleADs(oldAds: IterableOnce[Int], proportion: Double, random: Random): IndexedSeq[Int] = {
    require(proportion <= 1, f"proportion must be less than 1: proportion = ${proportion}")
    oldAds.iterator.toIndexedSeq.map(s => Range(0, s).iterator.map(_ => random.nextDouble()).count(_ < proportion))
  }

  /**
   * Re-genotypes a variant for each sample after downsampling the allele counts based on the given 
   * per-sample proportions.
   * @param variant the variant to downsample and re-genotype
   * @param proportions proportion to downsample the allele counts for each sample prior to re-genotyping
   * @param random random number generator for downsampling
   * @param epsilon the sequencing error rate for genotyping
   * @return a new variant with updated genotypes, downsampled ADs, and recomputed PLs
   */
  def downsampleAndRegenotype(variant: Variant,
                              proportions: Map[String, Double],
                              random: Random, epsilon: Double=0.01): Variant = {
    try {
      variant.copy(genotypes=variant.genotypes.map { case (sample, gt) =>
        val proportion = proportions(sample)
        sample -> downsampleAndRegenotype(gt=gt, proportion=proportion, random=random, epsilon=epsilon)
      })
    } catch {
      case e: MatchError => throw new Exception(
        "processing " + variant.id +  " at " + variant.chrom + ":" + variant.pos + "-" + variant.end, e
      )
    }
  }

  /**
   * Re-genotypes a sample after downsampling the allele counts based on the given proportion.
   * @param gt the genotype to downsample
   * @param proportion proportion to downsample the allele count prior to re-genotyping
   * @param random random number generator for downsampling
   * @param epsilon the sequencing error rate for genotyping
   * @return a new Genotype with updated allele depths, PLs, and genotype
   */
  def downsampleAndRegenotype(gt: Genotype, proportion: Double, random: Random, epsilon: Double): Genotype = {
    val oldAds = gt.getOrElse[IndexedSeq[Int]]("AD", throw new Exception(s"AD tag not found for sample ${gt.sample}"))
    val newAds = downsampleADs(oldAds, proportion, random)
    val likelihoods = Likelihoods(newAds)
    val pls = likelihoods.pls
    val calls = likelihoods.mostLikelyCall(gt.alleles.toSeq)
    gt.copy(attrs=Map("PL" -> pls, "AD" -> newAds, "DP" -> newAds.sum), calls=calls)
  }

  object Likelihoods {
    /**Converts a sequence of log-likelihoods to phred-scale by 1) multiplying each by -10, 2)
     * subtracting from each the min value so the smallest value is 0, and 3) rounding to the
     * nearest integer.
     */
    def logToPhredLikelihoods(logLikelihoods: IndexedSeq[Double]): IndexedSeq[Int] = {
        val rawPL = logLikelihoods.map(gl => gl * -10)
        val minPL = rawPL.min
        rawPL.map(pl => (pl - minPL).round.toInt)
    }

    /** Computes the likelihoods for each possible biallelic genotype.
     * @param alleleDepthA the reference allele depth
     * @param alleleDepthB the alternate allele depth
     * @param epsilon      the error rate for genotyping
     * @return a new `Likelihood` that has the likelihoods of AA, AB, and BB
     */
    def biallelic(alleleDepthA: Int, alleleDepthB: Int, epsilon: Double = 0.01): IndexedSeq[Double] = {
      val aGivenAA = log10(1 - epsilon)
      val aGivenBB = log10(epsilon)
      val aGivenAB = log10((1 - epsilon) / 2)

      val rawGlAA = ((alleleDepthA * aGivenAA) + (alleleDepthB * aGivenBB))
      val rawGlBB = ((alleleDepthA * aGivenBB) + (alleleDepthB * aGivenAA))
      val rawGlAB = ((alleleDepthA + alleleDepthB) * aGivenAB)

      IndexedSeq(rawGlAA, rawGlAB, rawGlBB)
    }

    /** Computes the likelihoods for each possible genotype given a sequence of read depths for any
     * number of alleles.
     * @param alleleDepths the sequence of allele depths in the order specified in the VCF
     * @param epsilon      the error rate for genotyping
     * @return a new `Likelihood` that has the log likelihoods of all possible genotypes in the
     * order specified in VFC spec for the GL/PL tags.
     */
    def generalized(alleleDepths: IndexedSeq[Int], epsilon: Double = 0.01): IndexedSeq[Double] = {
      val numAlleles = alleleDepths.length
      // probabilities associated with each possible genotype for a pair of alleles
      val logProbs: Array[Double] = Array(
        math.log10(epsilon),
        math.log10((1 - epsilon) / 2),
        math.log10(1 - epsilon)
      )
      // compute genotype log-likelihoods      
      (0 until numAlleles).flatMap(b =>
        (0 to b).map(a =>
          (0 until numAlleles).map(allele =>
            logProbs(Array(a, b).count(_ == allele)) * alleleDepths(allele)
          ).sum
        )
      )
    }

    def apply(alleleDepths: IndexedSeq[Int], epsilon: Double = 0.01): Likelihoods = {
      val numAlleles = alleleDepths.length
      require(numAlleles >= 2, "at least two alleles are required to calculate genotype likelihoods")
      Likelihoods(numAlleles, generalized(alleleDepths, epsilon))
    }
  }

  /** Stores the log10(likelihoods) for all possible genotypes.
   * @param numAlleles the number of alleles the variant has
   * @param genotypeLikelihoods sequence of GLs in the order specified in the VCF spec
   */
  case class Likelihoods(numAlleles: Int, genotypeLikelihoods: IndexedSeq[Double]) {
    /**
      * Returns the likelihoods as a list of phred-scaled integers (i.e, the value of the PL tag).
      * @return a list of phred-scaled likelihooodS for AA, AB, BB.
      */
    def pls: IndexedSeq[Int] = {
      Likelihoods.logToPhredLikelihoods(genotypeLikelihoods)
    }

    def mostLikelyGenotype: Option[(Int, Int)] = {
      val minIndexes = pls.zipWithIndex.filter(pair => pair._1 == 0)
      minIndexes.length match {
        case 0 => throw new RuntimeException("expected the most likely PL to have a value of 0.0")
        case 1 => {
          val genotypes = 
            for (b <- 0 until numAlleles; a <- 0 to b)
            yield (a, b)
          Some(genotypes(minIndexes.head._2))
        }
        case _ => None  // if multiple genotypes are most likely, don't make a call
      }
    }

    def mostLikelyCall(alleles: Seq[Allele]): IndexedSeq[Allele] = {
      mostLikelyGenotype match {
        case None => IndexedSeq(NoCallAllele, NoCallAllele)
        case Some((a, b)) => IndexedSeq(alleles(a), alleles(b))
      }
    }
  }
}

@clp(group=ClpGroups.VcfOrBcf, description =
  """
    |Re-genotypes a VCF after downsampling the allele counts.
    |
    |The input VCF must have at least one sample.
    |
    |If the input VCF contains a single sample, the downsampling target may be specified as a
    |proportion of the original read depth using `--proportion=(0..1)`, or as the combination of
    |the original and target _number of sequenced bases_ (`--originalBases` and
    |`--downsampleToBases`). For multi-sample VCFs, the downsampling target must be specified using
    |`--downsampleToBases`, and a metadata file with the total number of sequenced bases per sample
    |is required as well. The metadata file must follow the
    |[[https://www.internationalgenome.org/category/meta-data/] 1000 Genomes index format], but the
    |only required columns are `SAMPLE_NAME` and `BASE_COUNT`. A propportion for each sample is
    |calculated by dividing the _target number of sequenced bases_ by the _original number of
    |sequenced bases_.
    |
    |The tool first (optionally) winnows the VCF file to remove variants within a distance to each
    |other specified by `--window-size` (the default value of `0` disables winnowing). Next, each
    |sample at each variant is examined independently. The allele depths per-genotype are randoml
    |downsampled given the proportion. The downsampled allele depths are then used to re-compute
    |allele likelhoods and produce a new genotype.
    |
    |The tool outputs a downsampled VCF file with the winnowed variants removed, and with the
    |genotype calls and `DP`, `AD`, and `PL` tags updated for each sample at each retained variant.
""")
class DownsampleVcf
(@arg(flag='i', doc="The vcf to downsample.") val input: PathToVcf,
  @arg(flag='p', doc="Proportion of bases to retain (for single-sample VCF).") val proportion: Option[Double] = None,
  @arg(flag='b', doc="Original number of bases (for single-sample VCF).") val originalBases: Option[Double] = None,
  @arg(flag='m', doc="Index file with bases per sample.") val metadata: Option[FilePath] = None,
  @arg(flag='n', doc="Target number of bases to downsample to.") val downsampleToBases: Option[Double] = None,
  @arg(flag='o', doc="Output file name.") val output: PathToVcf,
  @arg(flag='w', doc="Winnowing window size.") val windowSize: Int = 0,
  @arg(flag='e', doc="Sequencing Error rate for genotyping.") val epsilon: Double = 0.01,
  @arg(flag='c', doc="True to write out no-calls.") val writeNoCall: Boolean = false,
  @arg(flag='s', doc="Random seed value.") val seed: Int = 42,
  ) extends FgBioTool {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  require(windowSize >= 0, "window size must be greater than or equal to zero")
  require(0 <= epsilon && epsilon <= 1, "epsilon/error rate must be between 0 and 1")
  (proportion, originalBases, metadata, downsampleToBases) match {
    case (Some(x), None, None, None) => 
      require(x > 0, "proportion must be greater than 0")
      require(x < 1, "proportion must be less than 1")
    case (None, Some(original), None, Some(target)) =>
      require(original > 0, "originalBases must be greater than zero")
      require(target > 0, "target base count must be greater than zero")
    case (None, None, Some(metadata), Some(target)) =>
      Io.assertReadable(metadata)
      require(target > 0, "target base count must be greater than zero")
    case (None, _, _, None) =>
      throw new IllegalArgumentException(
        "exactly one of proportion or downsampleToBases must be specified"
      )
    case _ =>
      throw new IllegalArgumentException(
        "exactly one of proportion, originalBases, or metadata must be specified"
      )
  }
  
  override def execute(): Unit = {
    val vcf = VcfSource(input)
    val proportions = (
      (proportion, originalBases, metadata, downsampleToBases) match {
        case (Some(x), None, None, None) =>
          require(vcf.header.samples.length == 1, "--original-bases requires a single-sample VCF")
          LazyList(vcf.header.samples.head -> x)
        case (None, Some(original), None, Some(target)) =>
          require(vcf.header.samples.length == 1, "--original-bases requires a single-sample VCF")
          LazyList(vcf.header.samples.head -> math.min(target / original, 1.0))
        case (None, None, Some(metadata), Some(target)) =>
          Sample.read(metadata)
            .filter(s => vcf.header.samples.contains(s.SAMPLE_NAME))
            .map(sample => sample.SAMPLE_NAME -> math.min(target / sample.BASE_COUNT.toDouble, 1.0))
        case _ =>
          throw new RuntimeException("unexpected parameter combination")
      }
    ).toMap
    proportions.foreach { case (s, p) => logger.info(f"Downsampling $s with proportion ${p}%.4f") }
    
    val inputProgress = ProgressLogger(logger, noun="variants read")
    val inputVariants = ProgressLogger.ProgressLoggingIterator(vcf.iterator).progress(inputProgress)
    val winnowed = if (windowSize > 0) {
      val winnowed = winnowVariants(inputVariants, windowSize=windowSize, dict=vcf.header.dict)
      val winnowedProgress = ProgressLogger(logger, noun="variants retained")
      ProgressLogger.ProgressLoggingIterator(winnowed).progress(winnowedProgress)
    } else {
      inputVariants
    }
    val outputVcf = VcfWriter(path=output, header=buildOutputHeader(vcf.header))

    val progress = ProgressLogger(logger, noun="variants written")
    val random = new Random(seed)
    winnowed.foreach { v =>
      val ds = downsampleAndRegenotype(v, proportions=proportions, random=random, epsilon=epsilon)
      if (writeNoCall || !ds.gts.forall(g => g.isNoCall)) {
        outputVcf += ds
        progress.record(ds)
      }
    }
    
    progress.logLast()
    vcf.safelyClose()
    outputVcf.close()
  }

  private def buildOutputHeader(in: VcfHeader): VcfHeader = {
    val fmts = Seq.newBuilder[VcfFormatHeader]
    fmts ++= in.formats

    if (!in.format.contains("AD")) {
      fmts += VcfFormatHeader(id="AD", count=VcfCount.OnePerAllele, kind=VcfFieldType.Integer, description="Per allele depths.")
    }

    if (!in.format.contains("DP")) {
      fmts += VcfFormatHeader(id="DP", count=VcfCount.Fixed(1), kind=VcfFieldType.Integer, description="Total depth across alleles.")
    }

    if (!in.format.contains("PL")) {
      fmts += VcfFormatHeader(id="PL", count=VcfCount.OnePerGenotype, kind=VcfFieldType.Integer, description="Per genotype phred scaled likelihoods.")
    }

    in.copy(formats=fmts.result())
  }
}

private object Sample {
  /** Load a set of samples from the 1KG metadata file. */
  def read(path: FilePath): Seq[Sample] = {
    val lines = Io.readLines(path).dropWhile(_.startsWith("##")).map(line => line.dropWhile(_ == '#'))
    Metric.read[Sample](lines=lines)
  }
}

case class Sample(ENA_FILE_PATH: String = ".",
                  MD5SUM: String = ".",
                  RUN_ID: String = ".",
                  STUDY_ID: String = ".",
                  STUDY_NAME: String = ".",
                  CENTER_NAME: String = ".",
                  SUBMISSION_ID: String = ".",
                  SUBMISSION_DATE: String = ".",
                  SAMPLE_ID: String = ".",
                  SAMPLE_NAME: String,
                  POPULATION: String = ".",
                  EXPERIMENT_ID: String = ".",
                  INSTRUMENT_PLATFORM: String = ".",
                  INSTRUMENT_MODEL: String = ".",
                  LIBRARY_NAME: String = ".",
                  RUN_NAME: String = ".",
                  INSERT_SIZE: String = ".",
                  LIBRARY_LAYOUT: String = ".",
                  PAIRED_FASTQ: String = ".",
                  READ_COUNT: String = ".",
                  BASE_COUNT: Long,
                  ANALYSIS_GROUP: String = ".") extends Metric
