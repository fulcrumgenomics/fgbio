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

package com.fulcrumgenomics.vcf

import java.nio.file.Paths
import java.util
import java.util.Comparator

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{GenomicSpan, Metric, NumericCounter, ProgressLogger}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.vcf.PhaseCigarOp.PhaseCigarOp
import dagr.commons.io.{Io, PathUtil}
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.util.{IntervalList, OverlapDetector}
import htsjdk.variant.variantcontext.{Genotype, GenotypeBuilder, VariantContext, VariantContextBuilder}
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.annotation.tailrec
import scala.collection.mutable.ListBuffer
import scala.collection.JavaConverters._

@clp(
  description =
    """
      |Assess the accuracy of phasing for a set of variants.
      |
      |All phased genotypes should be annotated with the "PS" (phase set) FORMAT tag, which by convention is the
      |position of the first variant in the phase set (see the VCF specification).  Furthermore, the alleles of a phased
      |genotype should use the '|' separator instead of the '/' separator, where the latter indicates the genotype is
      |unphased.
      |
      |The input VCFs are assumed to be single sample: the genotype from the first sample is used.
      |
      |Only bi-allelic heterozygous SNPs are considered.
      |
      |The input known phased variants can be subsetted using the known interval list, for example to keep only variants
      |from high-confidence regions.
      |
      |If the intervals argument is supplied, only the set of chromosomes specified will be analyzed.  Note that the full
      |chromosome will be analyzed and start/stop positions will be ignored.
    """,
  group=ClpGroups.VcfOrBcf
)
class AssessPhasing
( @arg(flag="c", doc="The VCF with called phased variants.") val calledVcf: PathToVcf,
  @arg(flag="t", doc="The VCF with known phased variants.") val truthVcf: PathToVcf,
  @arg(flag="o", doc="The output prefix for all output files.") val output: PathPrefix,
  @arg(flag="k", doc="The interval list over which known phased variants should be kept.") val knownIntervals: Option[PathToIntervals] = None,
  @arg(flag="m", doc="Allow missing fields in the VCF header.") val allowMissingFieldsInVcfHeader: Boolean = true,
  @arg(flag="s", doc="Skip sites where the truth and call are both called but do not share the same alleles.") val skipMismatchingAlleles: Boolean = true,
  @arg(flag="l", doc="Analyze only the given chromosomes in the interval list.  The entire chromosome will be analyzed (start and end ignored).") val intervals: Option[PathToIntervals] = None,
  @arg(flag="b", doc="Remove enclosed phased blocks and truncate overlapping blocks.") val modifyBlocks: Boolean = true,
  @arg(flag="d", doc="Output a VCF with the called variants annotated by if their phase matches the truth") val debugVcf: Boolean = false
) extends FgBioTool with LazyLogging {
  import AssessPhasing.{CalledSampleName, TruthSampleName}

  Io.assertReadable(Seq(calledVcf, truthVcf))
  knownIntervals.foreach(Io.assertReadable)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Setup
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // check the sequence dictionaries.
    val dict = {
      val calledReader = new VCFFileReader(calledVcf.toFile, true)
      val truthReader = new VCFFileReader(truthVcf.toFile, true)
      calledReader.getFileHeader.getSequenceDictionary.assertSameDictionary(truthReader.getFileHeader.getSequenceDictionary)
      calledReader.close()
      truthReader.close()
      calledReader.getFileHeader.getSequenceDictionary
    }

    val knownIntervalList        = knownIntervals.map { intv => IntervalList.fromFile(intv.toFile).uniqued() }
    val calledBlockLengthCounter = new NumericCounter[Long]()
    val truthBlockLengthCounter  = new NumericCounter[Long]()
    val metric                   = new AssessPhasingMetric
    val writer = if (!debugVcf) None else {
      val path    = Paths.get(output.toString + AssessPhasing.AnnotatedVcfExtension)
      val reader  = new VCFFileReader(calledVcf.toFile, true)
      val header  = reader.getFileHeader
      val builder = new VariantContextWriterBuilder()
        .setOutputFile(path.toFile)
        .setReferenceDictionary(header.getSequenceDictionary)
        .setOption(Options.INDEX_ON_THE_FLY)
        .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
      val writer: VariantContextWriter = builder.build
      val headerLines: util.Set[VCFHeaderLine] = new util.HashSet[VCFHeaderLine](header.getMetaDataInSortedOrder)
      headerLines.add(AssessPhasing.PhaseConcordanceFormatHeaderLine) // add the new format header lines
      VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS) // add standard header lines
      writer.writeHeader(new VCFHeader(headerLines, List("call", "truth").asJava))
      reader.safelyClose()
      Some(writer)
    }
    val chromosomes = intervals.map { intv =>
      val intervalList = IntervalList.fromFile(intv.toFile)
      // Developer Note: warn the user if the supplied intervals do not span entire chromosomes.
      intervalList.getIntervals.find { interval =>
        interval.getStart != 1 || interval.getEnd != dict.getSequence(interval.getContig).getSequenceLength
      }.foreach { interval =>
        logger.warning(s"Interval list (--intervals) given with intervals that do not span entire chromosomes (ex. '$interval').  Start/end will be ignored and entire chromosome analyzed.")
      }

      intervalList.getIntervals.map { i => i.getContig }.toSet
    }

    // NB: could parallelize!
    dict.getSequences
      .iterator
      .filter { sequence => chromosomes match {
          case Some(set) => set.contains(sequence.getSequenceName)
          case None      => true
        }
      }
      .foreach { sequence =>
      executeContig(
        dict                     = dict,
        contig                   = sequence.getSequenceName,
        contigLength             = sequence.getSequenceLength,
        knownIntervalList        = knownIntervalList,
        metric                   = metric,
        calledBlockLengthCounter = calledBlockLengthCounter,
        truthBlockLengthCounter  = truthBlockLengthCounter,
        writer                   = writer
      )
    }

    val calledBlockLengthMetrics = calledBlockLengthCounter.map { case (length, count) => new PhaseBlockLengthMetric(dataset=CalledSampleName, length=length, count=count) }
    val truthBlockLengthMetrics  = truthBlockLengthCounter.map { case (length, count) => new PhaseBlockLengthMetric(dataset=TruthSampleName, length=length, count=count) }
    val blockLengthMetrics       = (calledBlockLengthMetrics ++ truthBlockLengthMetrics).toSeq

    val calledAssemblyStats = AssemblyStatistics(calledBlockLengthCounter)
    val truthAssemblyStats  = AssemblyStatistics(truthBlockLengthCounter)

    metric.mean_called_block_length   = calledBlockLengthCounter.mean()
    metric.median_called_block_length = calledBlockLengthCounter.median()
    metric.stddev_called_block_length = calledBlockLengthCounter.stddev(m=metric.mean_called_block_length)
    metric.n50_called_block_length    = calledAssemblyStats.n50
    metric.n90_called_block_length    = calledAssemblyStats.n90
    metric.l50_called                 = calledAssemblyStats.l50
    metric.mean_truth_block_length    = truthBlockLengthCounter.mean()
    metric.median_truth_block_length  = truthBlockLengthCounter.median()
    metric.stddev_truth_block_length  = truthBlockLengthCounter.stddev(m=metric.mean_truth_block_length)
    metric.n50_truth_block_length     = truthAssemblyStats.n50
    metric.n90_truth_block_length     = truthAssemblyStats.n90
    metric.l50_truth                  = truthAssemblyStats.l50

    metric.finalizeMetric()

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Output the metrics and finish up
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    logger.info("Outputting")

    Metric.write(path=PathUtil.pathTo(output + AssessPhasingMetric.MetricExtension), metric=metric)
    Metric.write(PathUtil.pathTo(output + PhaseBlockLengthMetric.MetricExtension), blockLengthMetrics)

    writer.foreach(_.close())
  }

  private def executeContig(dict: SAMSequenceDictionary,
                            contig: String,
                            contigLength: Int,
                            knownIntervalList: Option[IntervalList],
                            metric: AssessPhasingMetric,
                            calledBlockLengthCounter: NumericCounter[Long],
                            truthBlockLengthCounter: NumericCounter[Long],
                            writer: Option[VariantContextWriter] = None
                           ): Unit = {
    logger.info(s"Assessing $contig")

    val intervalListForContig = knownIntervalList.map { oldList =>
      val newList = new IntervalList(oldList.getHeader)
      newList.addall(oldList.getIntervals.filter { _.getContig == contig }.toJavaList)
      newList
    }

    // get the phased blocks
    logger.info("Getting the called phase blocks")
    val calledPhaseBlockDetector = {
      val calledReader = new VCFFileReader(calledVcf.toFile, true)
      val detector = PhaseBlock.buildOverlapDetector(
        iterator     = toVariantContextIterator(calledReader, contig, contigLength),
        dict         = dict,
        modifyBlocks = modifyBlocks
      )
      calledReader.close()
      detector
    }
    logger.info("Getting the known phase blocks")
    val truthPhaseBlockDetector = {
      val truthReader = new VCFFileReader(truthVcf.toFile, true)
      val detector = PhaseBlock.buildOverlapDetector(
        iterator     = toVariantContextIterator(truthReader, contig, contigLength, intervalList=intervalListForContig),
        dict         = dict,
        modifyBlocks = modifyBlocks
      )
      truthReader.close()
      detector
    }

    // get an iterator of the pairs
    val calledReader   = new VCFFileReader(calledVcf.toFile, true)
    val truthReader    = new VCFFileReader(truthVcf.toFile, true)
    val pairedIterator = JointVariantContextIterator(
      iters = Seq(toVariantContextIterator(truthReader, contig, contigLength, intervalList=intervalListForContig), toVariantContextIterator(calledReader, contig, contigLength)),
      dict = dict
    ).map { case Seq(left, right) => (left, right) }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create the phasing cigar
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    logger.info("Creating the phasing CIGAR")
    val cigar = PhaseCigar(
      pairedIterator           = pairedIterator,
      truthPhaseBlockDetector  = truthPhaseBlockDetector,
      calledPhaseBlockDetector = calledPhaseBlockDetector,
      metric                   = metric,
      skipMismatchingAlleles   = skipMismatchingAlleles,
      writer                   = writer
    )

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get the number of short and long switch errors, and other metrics
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    logger.info("Computing short switch errors")
    metric.num_short_switch_errors = cigar.toShortSwitchErrorIndices().length

    // To find the # of long switch errors, we need to ignore runs of consecutive indices.
    logger.info("Computing long switch errors")
    val (numLongSwitchErrors, numLongSwitchSites) = cigar.toLongSwitchErrorsAndSites()
    metric.num_long_switch_errors += numLongSwitchErrors
    metric.num_switch_sites  += numLongSwitchSites

    // Use the calculation described here: http://dx.doi.org/10.1038%2Fng.3119
    logger.info("Computing Illumina switch errors")
    val illuminaSwitchErrors                = cigar.toIlluminaSwitchErrors()
    metric.num_illumina_point_switch_errors += illuminaSwitchErrors.numPointErrors
    metric.num_illumina_long_switch_errors  += illuminaSwitchErrors.numLongSwitchErrors
    metric.num_illumina_switch_sites   += illuminaSwitchErrors.numSites

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get the phase blocks
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    logger.info("Computing block metrics")

    calledPhaseBlockDetector.getAll.toIterator.foreach { block => calledBlockLengthCounter.count(block.length) }
    truthPhaseBlockDetector.getAll.toIterator.foreach { block => truthBlockLengthCounter.count(block.length) }

    calledReader.safelyClose()
    truthReader.safelyClose()

    logger.info(s"Completed $contig")
  }

  private def toVariantContextIterator(reader: VCFFileReader,
                                       contig: String,
                                       contigLength: Int,
                                       intervalList: Option[IntervalList] = None): Iterator[VariantContext] = {
    val sampleName = reader.getFileHeader.getSampleNamesInOrder.iterator().next()
    (intervalList match {
      case Some(intv) =>
        ByIntervalListVariantContextIterator(reader.iterator().toIterator, intv, dict=reader.getFileHeader.getSequenceDictionary)
      case None =>
        reader.query(contig, 1, contigLength).toIterator
    })
    .map(_.subContextFromSample(sampleName))
    .filter(v => v.isSNP && v.isBiallelic && v.getGenotype(sampleName).isHet)
  }
}

object AssessPhasing {
  /** The output sample name for the called variants in the debug VCF. */
  val CalledSampleName = "call"
  /** The output sample name for the truth/known variants in the debug VCF. */
  val TruthSampleName = "truth"

  val AnnotatedVcfExtension = ".assess_phasing.vcf.gz"

  val PhaseConcordanceFormatTag = "PHASE_CONC"
  val PhaseConcordanceFormatDescription = "The phase concordance (Match or Mismatch) determined by fgbio's AssessPhasing"
  val PhaseConcordanceFormatHeaderLine = new VCFFormatHeaderLine(PhaseConcordanceFormatTag, 1, VCFHeaderLineType.Integer, PhaseConcordanceFormatDescription)

  def getPhasingSetId(ctx: VariantContext): Int = {
    ctx.getGenotype(0).getAttributeAsInt("PS", -1) // Integer.valueOf(ctx.getGenotype(0).getExtendedAttribute("PS", -1).toString)
  }
}

private object AssemblyStatistics {
  type BlockLength = Long
  def apply(blockLengthCounter: NumericCounter[BlockLength]): AssemblyStatistics = {
    if (blockLengthCounter.isEmpty) return AssemblyStatistics(0, 0, 0)

    val numBases = blockLengthCounter.map { case (length, count) => length * count }.sum
    val fiftyPercent  = numBases / 2.0
    val ninetyPercent = numBases * 0.9

    val blockLengthsAndCounts = blockLengthCounter.toSeq.sortBy(_._1)
    var blockLengthSum = 0d
    var n50: BlockLength = 0
    var n90: BlockLength = 0
    var l50: Long = 0
    forloop (blockLengthsAndCounts.length - 1) (0 <= _) (_ - 1) { i =>
      val (length, count)   = blockLengthsAndCounts(i)

      // L50
      if (blockLengthSum < fiftyPercent) { // only if L50 has not been set yet
        val remainingSum      = fiftyPercent - blockLengthSum
        val numBlocksRequired = Math.ceil(remainingSum / length).toInt
        if (numBlocksRequired < count) {
          l50 += numBlocksRequired
        }
        else {
          l50 += count
        }
      }

      // ** IMPORTANT ** add the sum of bases in blocks of this size *after* updating the L50.
      blockLengthSum += length * count

      // N50
      if (n50 == 0 && blockLengthSum >= fiftyPercent) {
        n50 = length
      }

      // N90
      if (n90 == 0 && blockLengthSum >= ninetyPercent) {
        n90 = length
      }
    }

    AssemblyStatistics(n50=n50, n90=n90, l50=l50)
  }
}

/**
  * @param n50 the longest block length such that the bases covered by all blocks this length and longer are at least
  *            50% of the # of bases covered by all blocks.
  * @param n90 the longest block length such that the bases covered by all blocks this length and longer are at least
  *            90% of the # of bases covered by all blocks.
  * @param l50 the smallest number of blocks such that the sum of the lengths of the blocks is >= 50% of the sum of
  *            the lengths of all blocks.
  */
private case class AssemblyStatistics(n50: Long, n90: Long, l50: Long)

object PhaseBlockLengthMetric {
  val MetricExtension = PathUtil.replaceExtension(Paths.get(AssessPhasingMetric.MetricExtension), ".block_lengths.txt").toString
}

/** Provides the number of phased blocks of a given length.
  *
  * @param dataset the name of the dataset (ex. truth or call)
  * @param length the length of the phased block
  * @param count the number of phased blocks of the given length.
  */
case class PhaseBlockLengthMetric
( dataset: String,
  length: Long = 0,
  count: Long = 0
) extends Metric

object PhaseBlock extends LazyLogging {
  import scala.collection.mutable

  /** Creates an overlap detector for blocks of phased variants.  Variants from the same block are found using the
    * "PS" tag.  The modify blocks option resolves overlapping blocks
    */
  private[vcf] def buildOverlapDetector(iterator: Iterator[VariantContext], dict: SAMSequenceDictionary, modifyBlocks: Boolean = true): OverlapDetector[PhaseBlock] = {
    val detector = new OverlapDetector[PhaseBlock](0, 0)
    val progress = new ProgressLogger(logger)

    // create the blocks
    val phaseBlocks = mutable.HashMap[Int, PhaseBlock]()
    while (iterator.hasNext) {
      val ctx = iterator.next()
      val phaseSetId: Int = AssessPhasing.getPhasingSetId(ctx)
      if (phaseSetId > 0) {
        val phaseBlock = phaseBlocks.get(phaseSetId) match {
          case Some(block) =>
            require(block.getContig == ctx.getContig)
            require(block.getStart <= ctx.getStart)
            block.copy(end=Math.max(block.getEnd, ctx.getEnd)) // extend the block
          case None        =>
            new PhaseBlock(contig=ctx.getContig, start=ctx.getStart, end=ctx.getEnd)
        }
        phaseBlocks.put(phaseSetId, phaseBlock)
      }
      progress.record(ctx.getContig, ctx.getStart)
    }

    val blocksIn = new util.TreeSet[PhaseBlock](new Comparator[PhaseBlock] {
      /** Compares the two blocks based on start position, then returns the shorter block. */
      def compare(a: PhaseBlock, b: PhaseBlock): Int = {
        if (a.start < b.start) -1
        else if (a.start == b.start) b.length - a.length
        else 1
      }
    })

    phaseBlocks.values.foreach(blocksIn.add)
    logger.info(s"Found ${phaseBlocks.size} phase block")


    // Make sure we do not have any overlapping phase blocks
    // - if block #1 is enclosed in block #2, keep only block #2
    // - otherwise, if block #1 and #2 overlap, truncate the smaller one
    // The loop below will compare two blocks at a time: the two blocks with the smallest start positions.  If the do
    // not overlap, the first block (smaller position) is kept.  If one is fully-contained in the other, then the
    // enclosed block is discarded.  If they overlap, the smaller block is truncated such that they do not overlap.  In
    // the case the start position is changed due to truncation, then that block must be re-inserted into the list of
    // input blocks to guarantee we compare the two blocks with the smallest start positions.
    val blocksOut = ListBuffer[PhaseBlock]()
    var left = blocksIn.pollFirst()
    while (!blocksIn.isEmpty) {
      val right = blocksIn.pollFirst()
      require(left.getStart <= right.getStart, s"left: $left right: $right")
      // At this point, left has the smallest start position, and right as the next smallest start position.
      if (left.overlaps(right)) { // do they have any overlap?
        require(modifyBlocks, s"Block $left overlaps $right")
        if (left.encloses(right)) { // the right is enclosed in the left, so keep the left block only
          logger.info(s"Removing $right enclosed in $left")
          // do nothing (keep left == left) because left may overlap subsequent blocks too!
        }
        else if (right.encloses(left)) { // the left is enclosed in the right, so keep the right block only
          logger.info(s"Removing $left enclosed in $right")
          left = right
        }
        else if (left.length < right.length) { // the left is smaller, so truncate left
          logger.info(s"Truncating $left which overlaps $right");
          blocksOut.append(left.copy(end=right.getStart-1))
          left = right
        }
        else { // the right is smaller, so truncate right
          logger.info(s"Truncating $right which overlaps $left");
          blocksOut.append(left)
          // Since the start position of the right block is changed, and we are going to use it in the next iteration,
          // then re-insert it into the set and poll.
          blocksIn.add(right.copy(start=left.getEnd+1))
          left = blocksIn.pollFirst()
        }
      }
      else {
        blocksOut.append(left)
        left = right
      }
    }
    if (left != null) blocksOut.append(left)
    blocksOut.foreach { block => logger.info(s"Keeping $block") }

    // add any remaining
    detector.addAll(blocksOut.toList.asJava, blocksOut.toList.asJava)
    phaseBlocks.clear()

    detector
  }
}

case class PhaseBlock private (contig: String, start: Int, end: Int) extends GenomicSpan

object AssessPhasingMetric {
  val MetricExtension = ".assess_phasing_metrics.txt"
}

/** Some counts about phasing
  *
  * @param num_called the number of variants called.
  * @param num_phased the number of variants called with phase.
  * @param num_truth the number of variants with known truth genotypes.
  * @param num_truth_phased the number of variants with known truth genotypes with phase.
  * @param num_called_with_truth_phased the number of variants called that had a known phased genotype.
  * @param num_phased_with_truth_phased the number of variants called with phase that had a known phased genotype.
  * @param num_truth_phased_in_called_block the number of known phased variants that were in a called phased block.
  * @param num_both_phased_in_called_block the number of called phase variants that had a known phased genotype in a called phased block.
  * @param num_short_switch_errors the number of short switch errors (isolated switch errors).
  * @param num_long_switch_errors the number of long switch errors (# of runs of consecutive switch errors).
  * @param num_switch_sites the number of sites that could be (short or long) switch errors (i.e. the # of sites with both known and called phased variants).
  * @param num_illumina_point_switch_errors the number of point switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).
  * @param num_illumina_long_switch_errors the number of long switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).
  * @param num_illumina_switch_sites the number of sites that could be (point or long) switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).
  * @param frac_phased the fraction of called variants with phase.
  * @param frac_phased_with_truth_phased the fraction of known phased variants called with phase.
  * @param frac_truth_phased_in_called_block the fraction of phased known genotypes in a called phased block.
  * @param frac_phased_with_truth_phased_in_called_block the fraction of called phased variants that had a known phased genotype in a called phased block.
  * @param short_accuracy 1 - (num_short_switch_errors / num_switch_sites)
  * @param long_accuracy 1 - (num_long_switch_errors / num_switch_sites)
  * @param illumina_point_accuracy 1 - (num_illumina_point_switch_errors / num_illumina_switch_sites )
  * @param illumina_long_accuracy 1 - (num_illumina_long_switch_errors / num_illumina_switch_sites )
  * @param mean_called_block_length the mean phased block length in the callset.
  * @param median_called_block_length the median phased block length in the callset.
  * @param stddev_called_block_length the standard deviation of the phased block length in the callset.
  * @param n50_called_block_length the N50 of the phased block length in the callset.
  * @param n90_called_block_length the N90 of the phased block length in the callset.
  * @param l50_called the L50  of the phased block length in the callset.
  * @param mean_truth_block_length the mean phased block length in the truth.
  * @param median_truth_block_length the median phased block length in the truth.
  * @param stddev_truth_block_length the standard deviation of the phased block length in the truth.
  * @param n50_truth_block_length the N50 of the phased block length in the truth.
  * @param n90_truth_block_length the N90 of the phased block length in the callset.
  * @param l50_truth the L50 of the phased block length in the callset.
  */
case class  AssessPhasingMetric
(
  var num_called: Long = 0,
  var num_phased: Long = 0,
  var num_truth: Long = 0,
  var num_truth_phased: Long = 0,
  var num_called_with_truth_phased: Long = 0,
  var num_phased_with_truth_phased: Long = 0,
  var num_truth_phased_in_called_block: Long = 0,
  var num_both_phased_in_called_block: Long = 0,
  var num_short_switch_errors: Long = 0,
  var num_long_switch_errors: Long = 0,
  var num_switch_sites: Long = 0,
  var num_illumina_point_switch_errors: Long = 0,
  var num_illumina_long_switch_errors: Long = 0,
  var num_illumina_switch_sites: Long = 0,
  var frac_phased: Double = 0,
  var frac_phased_with_truth_phased: Double = 0,
  var frac_truth_phased_in_called_block: Double = 0,
  var frac_phased_with_truth_phased_in_called_block: Double = 0,
  var short_accuracy: Double = 0,
  var long_accuracy: Double = 0,
  var illumina_point_accuracy: Double = 0,
  var illumina_long_accuracy: Double = 0,
  var mean_called_block_length: Double = 0,
  var median_called_block_length: Double = 0,
  var stddev_called_block_length: Double = 0,
  var n50_called_block_length: Double = 0,
  var n90_called_block_length: Double = 0,
  var l50_called: Double = 0,
  var mean_truth_block_length: Double = 0,
  var median_truth_block_length: Double = 0,
  var stddev_truth_block_length: Double = 0,
  var n50_truth_block_length: Double = 0,
  var n90_truth_block_length: Double = 0,
  var l50_truth: Double = 0
) extends Metric {
  private def divide(a: Double, b: Double): Double = {
    if (b == 0) 0
    else a / b
  }
  def finalizeMetric(): this.type = {
    this.frac_phased                                   = divide(this.num_phased, this.num_called)
    this.frac_phased_with_truth_phased                 = divide(this.num_phased_with_truth_phased, this.num_called_with_truth_phased)
    this.frac_truth_phased_in_called_block             = divide(this.num_truth_phased_in_called_block, this.num_truth_phased)
    this.frac_phased_with_truth_phased_in_called_block = divide(this.num_both_phased_in_called_block, this.num_truth_phased_in_called_block)
    this.short_accuracy                                = 1.0 - divide(this.num_short_switch_errors, this.num_switch_sites)
    this.long_accuracy                                 = 1.0 - divide(this.num_long_switch_errors, this.num_switch_sites)
    this.illumina_point_accuracy                       = 1.0 - divide(this.num_illumina_point_switch_errors, this.num_illumina_switch_sites )
    this.illumina_long_accuracy                        = 1.0 - divide(this.num_illumina_long_switch_errors, this.num_illumina_switch_sites )
    this
  }
}

/** The cigar Elements */
private object PhaseCigarOp extends Enumeration {
  type PhaseCigarOp = Value
  val Match,      // 0
  Mismatch,       // 1
  TruthOnly,      // 2
  CallOnly,       // 3
  TruthEnd,       // 4
  CallEnd,        // 5
  BothEnd = Value // 6
}

private[vcf] object PhaseCigar {
  import PhaseCigarOp._
  import AssessPhasing.{CalledSampleName, TruthSampleName}

  private type VCtx = VariantContext

  def apply(cigar: Seq[PhaseCigarOp]): PhaseCigar = new PhaseCigar(cigar)

  /** Creates a Cigar from an iterator over variant contexts from a truth and call sample.  The variants should be on a
    * single chromosome only.
    *
    * If `skipMismatchingAlleles`, then it skips sites where both truth and call have a variant call, but the alleles
    * disagree between the two.
    *
    * Examples of short and long errors:
    *   CALL :  MMMMXXXMMM   MMMM|FFFFFFF|MMMMMFMMMM|F|MMFFFMMM
    *   TRUTH:  MMMM   MMMMMMMMMM-MMMMMMM|MMMMMMMMMM|F|MMMMMMMM
    *             correct       long        short     long
    *   CIGAR:  00003330002220000511111116000001000060600111000
    * shows one known block split
    */
  def apply(pairedIterator: Iterator[(Option[VariantContext], Option[VariantContext])],
            truthPhaseBlockDetector: OverlapDetector[PhaseBlock],
            calledPhaseBlockDetector: OverlapDetector[PhaseBlock],
            metric: AssessPhasingMetric,
            skipMismatchingAlleles: Boolean,
            writer: Option[VariantContextWriter] = None,
            assumeFixedAlleleOrder: Boolean = false): PhaseCigar = {
    val iter = pairedIterator.filter {
      // ensure the alleles are the same when both truth and call are called
      case ((Some(t: VCtx), Some(c: VCtx))) =>
        !skipMismatchingAlleles || t.getAlleles.toSet == c.getAlleles.toSet
      case _ => true
    }.flatMap { case (t, c) => // collect metrics but only keep sites where either variant context (i.e. truth or call) is phased.
      val ctxs                  = Seq(t, c).flatten // NB: can be 1 or 2 contexts here
      val inTruthPhaseBlock     = ctxs.headOption.exists(truthPhaseBlockDetector.overlapsAny)
      val inCalledPhaseBlock    = ctxs.headOption.exists(calledPhaseBlockDetector.overlapsAny)
      val isTruthVariant        = t.isDefined
      val isCalledVariant       = c.isDefined
      val isTruthVariantPhased  = t.exists { ctx => AssessPhasing.getPhasingSetId(ctx) > 0 }
      val isCalledVariantPhased = c.exists { ctx => AssessPhasing.getPhasingSetId(ctx) > 0 }

      // Fill in the metrics
      // Basic metrics
      if (isCalledVariant) {
        metric.num_called += 1
        if (isCalledVariantPhased) metric.num_phased += 1
      }
      if (isTruthVariant) {
        metric.num_truth += 1
        if (isTruthVariantPhased) metric.num_truth_phased += 1
      }
      // Less fun  compute metrics
      if (isCalledVariant && isTruthVariantPhased) {
        metric.num_called_with_truth_phased += 1
        if (isCalledVariantPhased) metric.num_phased_with_truth_phased += 1
      }
      if (inCalledPhaseBlock) {
        if (isTruthVariantPhased) {
          metric.num_truth_phased_in_called_block += 1
          if (isCalledVariantPhased) metric.num_both_phased_in_called_block += 1
        }
      }

      if (isCalledVariantPhased || isTruthVariantPhased) {
        Some((t, c))
      }
      else None
    }

    // The assignment of maternal/paternal alleles in some cases is arbitrary (ex. without a pedigree).  Therefore,
    // the first variant in any block (truth or called) that has a corresponding phased variant in the other sample
    // (called or truth) is assumed to match and therefore determines the order of the alleles in the subsequent
    // variants, until an end of block is reached for either sample.

    var applyInvertMatch: Boolean       = false
    var invertMatch: Boolean            = false

    val cigarOps = new ListBuffer[PhaseCigarOp]()
    while (iter.hasNext) {
      val (t, c) = iter.next()
      val truthPhasedVariant    = t.filter { ctx => AssessPhasing.getPhasingSetId(ctx) > 0 }
      val calledPhasedVariant   = c.filter { ctx => AssessPhasing.getPhasingSetId(ctx) > 0 }

      // Get the end-block operator if we encountered one.
      val blockEndOp = contextsToBlockEndOperator(truthPhasedVariant, calledPhasedVariant)

      // If we reach the end of a block, we need to infer if we should invert the match/mismatches based on the next
      // site where both called and truth have a phased variant.
      if (!assumeFixedAlleleOrder && blockEndOp.nonEmpty) {
         applyInvertMatch = false
      }

      // Get the cigar operator for the current variant site.
      val matchingOp = contextsToMatchingOperator(truthPhasedVariant, calledPhasedVariant).map { op =>
        if (assumeFixedAlleleOrder || (op != Match && op != Mismatch)) { // nothing to do
          op
        }
        else {
          // Check if the invert status has been set, if not set it, otherwise, apply it
          if (!applyInvertMatch) { // the first element!
            invertMatch = op == Mismatch
            applyInvertMatch = true
            Match
          }
          else if (invertMatch) { // invert it
            if (op == Match) Mismatch else Match
          }
          else op // keep it the way it is
        }
      }


      // Add it to the current buffer of cigar operators
      if (cigarOps.isEmpty) {
        cigarOps.append(BothEnd)
        matchingOp.foreach(cigarOps.append(_))
      }
      else {
        Seq(blockEndOp, matchingOp).flatten.foreach(cigarOps.append(_))
      }

      // Write to the output debug variant file
      (matchingOp, t, c, writer) match {
        case (Some(op), Some(truthCtx), Some(calledCtx), Some(w)) =>
          val calledGenotype = {
            val builder = new GenotypeBuilder(calledCtx.getGenotype(0))
            builder.name(CalledSampleName)
            builder.attribute(AssessPhasing.PhaseConcordanceFormatTag, op.toString)
            builder.make()
          }
          val truthGenotype = {
            val builder = new GenotypeBuilder(truthCtx.getGenotype(0))
            builder.name(TruthSampleName)
            builder.make()
          }
          val ctxBuilder = new VariantContextBuilder(calledCtx)
          ctxBuilder.genotypes(calledGenotype, truthGenotype)
          w.add(ctxBuilder.make())
        case _ => Unit
      }
    }

    if (cigarOps.isEmpty) {
      cigarOps.append(BothEnd)
    }
    cigarOps.append(BothEnd)

    PhaseCigar(cigarOps.toIndexedSeq)
  }

  /** Returns an end operator if we have reached a new block. */
  private[vcf] def contextsToBlockEndOperator(truth: Option[VariantContext], call: Option[VariantContext]): Option[PhaseCigarOp] = (truth, call) match {
    case (None, Some(c: VCtx))          => if (isStartOfPhaseBlock(c)) Some(CallEnd) else None
    case (Some(t: VCtx), None)          => if (isStartOfPhaseBlock(t)) Some(TruthEnd) else None
    case (Some(t: VCtx), Some(c: VCtx)) => (isStartOfPhaseBlock(t), isStartOfPhaseBlock(c)) match {
      case (true, true)   => Some(BothEnd)
      case (false, true)  => Some(CallEnd)
      case (true, false)  => Some(TruthEnd)
      case (false, false) => None
    }
    case _                              => unreachable()
  }

  /** Returns the cigar for the two variant contexts. */
  private[vcf] def contextsToMatchingOperator(truth: Option[VariantContext], call: Option[VariantContext]): Option[PhaseCigarOp] = (truth, call) match {
    case (None, Some(c: VCtx))          => Some(CallOnly)
    case (Some(t: VCtx), None)          => Some(TruthOnly)
    case (Some(t: VCtx), Some(c: VCtx)) => Some(cigarTypeForVariantContexts(t, c))
    case _                              => unreachable()
  }

  /** True if the phasing set id is the same as the start position of the given variant, false otherwise. */
  def isStartOfPhaseBlock(ctx: VariantContext): Boolean = ctx.getStart == AssessPhasing.getPhasingSetId(ctx)

  /** Computes the cigar for two variant contexts.  Returns [[Match]] if they share the same alleles in the same order,
    * [[Mismatch]] otherwise.
    */
  private[vcf] def cigarTypeForVariantContexts(truth: VariantContext, call: VariantContext): PhaseCigarOp = {
    val truthAlleles  = truth.getGenotype(0).getAlleles.toSeq
    val calledAlleles = call.getGenotype(0).getAlleles.toSeq
    require(truthAlleles.length == calledAlleles.length)
    require(truthAlleles.length == 2)
    if (truthAlleles.head != calledAlleles.head || truthAlleles.last != calledAlleles.last) PhaseCigarOp.Mismatch
    else PhaseCigarOp.Match
  }

  case class IlluminaSwitchErrors(var numPointErrors: Int, var numLongSwitchErrors: Int, var numSites: Int) {
    require(numPointErrors <= numSites)
    require(numLongSwitchErrors <= numSites)
    def add(other: IlluminaSwitchErrors): this.type = {
      this.numPointErrors      += other.numPointErrors
      this.numLongSwitchErrors += other.numLongSwitchErrors
      this.numSites            += other.numSites
      this
    }
  }
}


private[vcf] class PhaseCigar private(val cigar: Seq[PhaseCigarOp]) {
  import PhaseCigar._
  import PhaseCigarOp._

  /** Partitions the cigar into multiple cigars, with each cigar being a contiguous block of phased variants. */
  private[vcf] def toPhasedBlocks(isTruth: Boolean = false): Seq[PhaseCigar] = {
    val endCigarTypes = if (isTruth) Set(TruthEnd, BothEnd) else Set(CallEnd, BothEnd)

    // Splits the cigar any time we reach a PhaseCigarOp in the endCigarTypes
    val (cigarsToReturn, lastCigar) = this.cigar.foldLeft((ListBuffer[Seq[PhaseCigarOp]](), ListBuffer[PhaseCigarOp]())) {
      case ((previousCigars: ListBuffer[Seq[PhaseCigarOp]], currentCigar: ListBuffer[PhaseCigarOp]), phaseCigarOp: PhaseCigarOp) =>
        phaseCigarOp match {
          case tpe if endCigarTypes.contains(tpe) => // split!
            if (currentCigar.nonEmpty) previousCigars.append(currentCigar.toList)
            (previousCigars, new ListBuffer[PhaseCigarOp]())
          case _ => // keep going
            currentCigar.append(phaseCigarOp)
            (previousCigars, currentCigar)
        }
    }
    // Make sure to ge the last one
    if (lastCigar.nonEmpty) cigarsToReturn.append(lastCigar)

    cigarsToReturn.map(PhaseCigar(_)).toIndexedSeq
  }


  /** Gets all the indices in the cigar where there is a short switch error.
    *
    * Short switches are sites where the alleles are misphased (flipped), and both upstream and downstream, either the
    * haplotype block ends, or the next variant has the same phase.
    */
  def toShortSwitchErrorIndices(): Seq[Int] = {
    // check if we can find a match or the end of a block at the given idx.  If we cannot, move to the next idx.
    @tailrec
    def checkNext(cigar: Seq[PhaseCigarOp], idx: Int, by: Int): Boolean = {
      if (idx < 0 || cigar.length <= idx) return true // no more
      cigar(idx) match {
        case Match | CallEnd | BothEnd | TruthEnd => true
        case Mismatch     => false
        case CallOnly | TruthOnly => checkNext(cigar=cigar, idx=idx+by, by=by)
      }
    }
    // find cigars that have a mismatch at the index, but then anchored on either side
    this.cigar.indices.filter { i => this.cigar(i) == Mismatch && checkNext(this.cigar, i-1, -1) && checkNext(this.cigar, i+1, 1) }.toList
  }

  /** Gets the number of long switch errors and the number of sites examined.
    *
    * @return a tuple of the number of long switch errors and the total number of sites examined. Long switches are a
    *         stretch of two or more consecutive sites where the alleles are mis-phased (flipped), ann both upstream and
    *         downstream, either the haplotype block ends, or the next variant has the same phase. The number of sites
    *         examined are places where there are both a phased truth and phased called variant.
    */
  def toLongSwitchErrorsAndSites(): (Int, Int) = {
    // Gets all the indices in the cigar where there is a long switch error. Currently returns the index for all cigar
    // mismatches that contribute to a long switch error, and tries to maximize such a stretch.
    val iter = cigar.iterator.bufferBetter
    var numLongSwitchErrors = 0
    while (iter.hasNext) {
      if (iter.head == Mismatch) {
        val len = iter.takeWhile(c => c == Mismatch || c == CallOnly || c == TruthOnly).count(_ == Mismatch)
        if (len > 1) numLongSwitchErrors += 1
      }
      iter.dropWhile(_ != Mismatch)
    }

    val numSites = this.cigar.count {
      case Mismatch | Match => true
      case _                => false
    }

    (numLongSwitchErrors, numSites)
  }

  /** Computes the number of point and long switch error rates as described in http://dx.doi.org/10.1038%2Fng.3119 */
  def toIlluminaSwitchErrors(pointPenalty: Int = -1, transitionPenalty: Int = -5): IlluminaSwitchErrors = {

    // 1. partition the cigar into call blocks.
    // 2. for each call block, run the HMM, and sum the # of long swith errors
    this.toPhasedBlocks(isTruth=false).map { subCigar =>
      subCigar.toIlluminaSwitchErrorsHmm(pointPenalty=pointPenalty, transitionPenalty=transitionPenalty)
    }.foldLeft(IlluminaSwitchErrors(0, 0, 0)) { case (acc, cur) => acc.add(cur) }
  }

  /** Computes the number of point and long switch error rates as described in http://dx.doi.org/10.1038%2Fng.3119
    *
    * The HMM with a simple scoring system and two hidden states, one for each parental haplotype.  The emission
    * probability is scored as follows: 0 for being on the correct haplotype, otherwise `pointPenalty`.  The transition
    * probability is scored as follows: 0 for staying on the same haplotype, otherwise `transitionPenalty`.
    *
    * We prefer a long switch error over a sequence of point errors (i.e. when both produce have the same score). */
  private[vcf] def toIlluminaSwitchErrorsHmm(pointPenalty: Int = -1, transitionPenalty: Int = -5): IlluminaSwitchErrors = {
    if (this.cigar.isEmpty) throw new IllegalArgumentException("cigar was empty")

    // Get the initial phase: true is on the top haplotype, false is the bottom haplotype.
    val phase = this.cigar.filter { phaseCigarOp => phaseCigarOp == Match || phaseCigarOp == Mismatch }.map { phaseCigarOp =>
      if (phaseCigarOp == Match) true
      else false
    }
    if (phase.isEmpty) return IlluminaSwitchErrors(0, 0, 0)

    // one score tuple for each site, and each tuple is the score for being on a given haplotype
    val scores = ListBuffer.range(0, phase.length, 1).map(_ => (0, 0))
    val from   = ListBuffer.range(0, phase.length, 1).map(_ => (true, true)) // true if we stayed

    def topEmission(idx: Int): Int = if (phase(idx)) 0 else pointPenalty
    def botEmission(idx: Int): Int = if (phase(idx)) pointPenalty else 0

    // we arbitrarily assume we are on the top haplotype (?)
    // run viterbi for i == 0
    scores(0) = (topEmission(0), botEmission(0))
    from(0) = (true, true) // not actually needed

    // run Viterbi for i > 0
    var i = 1
    while (i < scores.length) {
      // top haplotype
      val (topScore, topFrom) = {
        val transition = scores(i - 1)._2 + transitionPenalty + topEmission(i)
        val stay = scores(i - 1)._1 + topEmission(i) // no penalty to stay :)
        if (transition < stay) (stay, true) // '<' makes it prefer transition
        else (transition, false)
      }

      // bottom haplotype
      val (botScore, botFrom) = {
        val transition = scores(i - 1)._1 + transitionPenalty + botEmission(i)
        val stay = scores(i - 1)._2 + botEmission(i) // no penalty to stay :)
        if (transition < stay) (stay, true) // '<' makes it prefer transition
        else (transition, false)
      }

      // update
      scores(i) = (topScore, botScore)
      from(i)   = (topFrom, botFrom)
      //println(s"scores($i): " + scores(i) + s" from($i): " + from(i))
      i += 1
    }

    /** Computes the # of long switch errors in the HMM output phased list. */
    def numLongSwitchErrors(states: List[Int]): Int = {
      if (states.length == 1) 0
      else states.sliding(2).count { case Seq(first, last) => first != last }
    }

    /** Backtracks to find the haplotype states, in reverse order. */
    def backtrack(haplotypeInit: Int, stayInit: Boolean): List[Int] = {
      var haplotype = haplotypeInit
      var stay      = stayInit
      val states    = new ListBuffer[Int]()
      i = scores.length - 2
      while (0 <= i) {
        //println("i: " + i + " haplotype: " + haplotype + " stay: " + stay + s" from($i): " + from(i))
        states.append(haplotype)
        if (stay) {
          stay = if (haplotype == 0) from(i)._1 else from(i)._2
        }
        else {
          haplotype = 1 - haplotype
          stay = if (haplotype == 0) from(i)._1 else from(i)._2
        }
        i -= 1
      }
      states.append(haplotype)
      require(states.length == scores.length)
      // NB: states is the reverse
      states.reverse.toList
    }

    // backtrack
    val states = {
      if (scores.last._1 > scores.last._2) backtrack(0, from(scores.length-1)._1)
      else if (scores.last._1 == scores.last._2) {
        // choose the solution with the most # of long switch errors
        val statesTop = backtrack(0, from(scores.length-1)._1)
        val statesBot = backtrack(1, from(scores.length-1)._2)
        val topScore  = numLongSwitchErrors(statesTop)
        val botScore  = numLongSwitchErrors(statesBot)
        if (topScore >= botScore) statesTop else statesBot
      }
      else backtrack(1, from(scores.length-1)._2)
    }

    //println("haplotypes: "  + states.toSeq.mkString(", "))
    //println(s"states.length=${states.length} scores.length=${scores.length}")

    // count how many mismatches are accounted for by point switch errors
    val numPointErrors = {
      require(phase.length == states.length)
      phase.zip(states.map(_ == 0)).count { case (p, s) => p != s }
    }

    IlluminaSwitchErrors(numPointErrors, numLongSwitchErrors(states), states.length)
  }
}
