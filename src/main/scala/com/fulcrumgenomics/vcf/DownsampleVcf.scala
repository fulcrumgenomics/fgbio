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

object DownsampleVcf extends LazyLogging {
  /** Removes variants that are within a specified distance from a previous variant
   * The end position of the current variant is compared with the start position of the following variant
   *
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

      def isInOrder(current: Variant, next: Variant, currentIndex: Int, nextIndex: Int): Boolean = {
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

  /** Downsamples variants using Allele Depths
   *
   * @param oldAds     an indexed seq of the original allele depths
   * @param proportion the proportion to use for downsampling,
   *                   calculated using total base count from the index and a target base count
   * @return a new IndexedSeq of allele depths of the same length as `oldAds`
   */
  def downsampleADs(oldAds: IndexedSeq[Int], proportion: Double, random: Random): IndexedSeq[Int] = {
    require(proportion <= 1, f"proportion must be less than 1: proportion = ${proportion}")
    oldAds.map(s => Range(0, s).iterator.map(_ => random.nextDouble()).count(_ < proportion))
  }

  /**
   * Does the downsampling on a Variant
   * @param variant the variant with the genotype to downsample
   * @param proportions a map of downsampling target proportions for each sample
   * @param random random number generator for downsampling
   * @param epsilon the error rate for genotyping
   * @return a new variant with updated genotypes
   */
  // Returns a new variant that has downsampled ADs, recomputed PLs and updated genotypes
  def downsampleAndRegenotype(variant: Variant, proportions: Map[String, Double], random: Random, epsilon: Double = 0.01): Variant = {
    try {
      variant.copy(genotypes = variant.genotypes.map { case (sample, gt) =>
        val proportion = proportions(sample)
        sample -> downsampleAndRegenotype(gt = gt, proportion = proportion, random = random, epsilon = epsilon)
      })
    } catch {
      case e: MatchError => throw new Exception(
        "processing " + variant.id +  " at " + variant.chrom + ":" + variant.pos + "-" + variant.end, e
      )
    }
  }

  /**
   * Does the downsampling on a Genotype
   * @param gt the genotype to downsample
   * @param proportion the proportion to use for downsampling allele depths
   * @param random random number generator for downsampling
   * @param epsilon the error rate for genotyping
   * @return a new Genotype with updated allele depths, PLs and genotype
   */
  def downsampleAndRegenotype(gt: Genotype, proportion: Double, random: Random, epsilon: Double): Genotype = {
    val oldAds = gt[IndexedSeq[Int]]("AD")
    val newAds = downsampleADs(oldAds, proportion, random)
    val Seq(aa, ab, bb) = computePls(newAds)
    val Seq(alleleA, alleleB) = gt.alleles.toSeq

    val calls = {
      if (aa == 0 && ab == 0 && bb == 0) IndexedSeq(NoCallAllele, NoCallAllele)
      else if (aa < ab && aa < bb) IndexedSeq(alleleA, alleleA)
      else if (bb < ab && bb < aa) IndexedSeq(alleleB, alleleB)
      else IndexedSeq(alleleA, alleleB)
    }
    gt.copy(attrs = Map("PL" -> IndexedSeq(aa, ab, bb), "AD" -> newAds, "DP" -> newAds.sum), calls = calls)
  }

  /**
   * Compute the genotype likelihoods given the allele depths.
   * @param ads The allele depths to generate likelihoods from
   * @return a list of three likelihoods
   */
  def computePls(ads: IndexedSeq[Int]): IndexedSeq[Int] = {
    val likelihoods = Likelihoods(ads(0), ads(1))
    IndexedSeq(likelihoods.aa.round.toInt, likelihoods.ab.round.toInt, likelihoods.bb.round.toInt)
  }


  object Likelihoods {
    /** Computes the likelihoods for each possible genotype.
     *
     * @param alleleDepthA the reference allele depth
     * @param alleleDepthB the alternate allele depth
     * @param epsilon      the error rate for genotyping
     * @return a new `Likelihood` that has the likelihoods of AA, AB, and BB
     */
    def apply(alleleDepthA: Int, alleleDepthB: Int, epsilon: Double = 0.01): Likelihoods = {
      val aGivenAA = log10(1 - epsilon)
      val aGivenBB = log10(epsilon)
      val aGivenAB = log10((1 - epsilon) / 2)

      val rawGlAA = ((alleleDepthA * aGivenAA) + (alleleDepthB * aGivenBB)) * -10
      val rawGlBB = ((alleleDepthA * aGivenBB) + (alleleDepthB * aGivenAA)) * -10
      val rawGlAB = ((alleleDepthA + alleleDepthB) * aGivenAB) * -10

      val minGL = math.min(math.min(rawGlAA, rawGlAB), rawGlBB)

      Likelihoods(
        aa = rawGlAA - minGL,
        ab = rawGlAB - minGL,
        bb = rawGlBB - minGL
      )
    }
  }

  /** Stores the log10(likelihoods) for all possible bi-allelic genotypes.
   *
   * @param aa likelihood of AA
   * @param ab likelihood of AB
   * @param bb likelihood of BB
   */
  case class Likelihoods(aa: Double, ab: Double, bb: Double) {
    def pls = IndexedSeq(aa.round.toInt, ab.round.toInt, bb.round.toInt)
  }
}

  @clp(group = ClpGroups.VcfOrBcf, description =
    """
      |DownsampleVcf takes a vcf file and metadata with sequencing info and
      |1. winnows the vcf to remove variants within a specified distance to each other,
      |2. downsamples the variants using the provided allele depths and target base count by
      |   re-computing/downsampling the allele depths for the new target base count
      |   and re-computing the genotypes based on the new allele depths
      |and writes a new downsampled vcf file.
      |For single-sample VCFs, the metadata file can be omitted, and instead you can specify originalBases.
  """)
  class DownsampleVcf
  (@arg(flag = 'i', doc = "The vcf to downsample.") val input: PathToVcf,
   @arg(flag = 'm', doc = "Index file with bases per sample.") val metadata: Option[FilePath] = None,
   @arg(flag = 'b', doc = "Original number of bases (for single-sample VCF)") val originalBases: Option[Double] = None,
   @arg(flag = 'n', doc = "Target number of bases to downsample to.") val downsampleToBases: Double,
   @arg(flag = 'o', doc = "Output file name.") val output: PathToVcf,
   @arg(flag = 'w', doc = "Winnowing window size.") val windowSize: Int = 150,
   @arg(flag = 'e', doc = "Error rate for genotyping.") val epsilon: Double = 0.01,
   @arg(flag = 'c', doc = "True to write out no-calls.") val writeNoCall: Boolean = false)
    extends FgBioTool {
    Io.assertReadable(input)
    Io.assertReadable(metadata)
    Io.assertCanWriteFile(output)
    require(downsampleToBases > 0, "target base count must be greater than zero")
    require(windowSize >= 0, "window size must be greater than or equal to zero")
    require(0 <= epsilon && epsilon <= 1, "epsilon/error rate must be between 0 and 1")
    originalBases match {
      case Some(x) =>
        require(x > 0, "originalBases must be greater than zero")
        require(metadata.isEmpty, "Must pass either originalBases (for single-sample VCF) or metadata, not both")
      case None =>
        require(metadata.isDefined, "Must pass either originalBases (for single-sample VCF) or metadata, not both")
    }

    override def execute(): Unit = {
      val vcf = VcfSource(input)
      val progress = ProgressLogger(logger, noun = "variants")
      val proportions = (
        originalBases match {
          case Some(x) =>
            require(vcf.header.samples.length == 1, "--original-bases requires a single-sample VCF")
            LazyList(vcf.header.samples.head -> math.min(downsampleToBases / x, 1.0))
          case _ =>
            Sample.read(metadata.getOrElse(throw new RuntimeException))
              .filter(s => vcf.header.samples.contains(s.SAMPLE_NAME))
              .map(sample => sample.SAMPLE_NAME -> math.min(downsampleToBases / sample.BASE_COUNT.toDouble, 1.0))
        }
      ).toMap
      proportions.foreach { case (s, p) => logger.info(f"Downsampling $s with proportion ${p}%.4f") }

      val winnowed = if (windowSize > 0) winnowVariants(vcf.iterator, windowSize = windowSize, dict = vcf.header.dict) else vcf.iterator
      val outputVcf = VcfWriter(path = output, header = buildOutputHeader(vcf.header))

      val random = new Random(42)
      winnowed.foreach { v =>
        val ds = downsampleAndRegenotype(v, proportions = proportions, random = random, epsilon = epsilon)
        if (writeNoCall) {
          outputVcf += ds
          progress.record(ds)
        }
        else if (!ds.gts.forall(g => g.isNoCall)) {
          outputVcf += ds
          progress.record(ds)
        }
      }
      
      progress.logLast()
      vcf.safelyClose()
      outputVcf.close()
    }

    def buildOutputHeader(in: VcfHeader): VcfHeader = {
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

      in.copy(formats = fmts.result())
    }
  }

object Sample {
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
