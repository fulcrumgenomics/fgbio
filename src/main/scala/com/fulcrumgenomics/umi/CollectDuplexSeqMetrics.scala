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
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.ConsensusTags.{MolecularId => MI, UmiBases => RX}
import com.fulcrumgenomics.util._
import htsjdk.samtools.util.{Interval, IntervalList, Murmur3, OverlapDetector}
import org.apache.commons.math3.distribution.BinomialDistribution

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.util.Failure

/**
  * Companion object for CollectDuplexSeqMetrics that contains various constants and types,
  * including all the various [[Metric]] sub-types produced by the program.
  */
object CollectDuplexSeqMetrics {
  // File extensions for all the files produced
  val FamilySizeMetricsExt: String       = ".family_sizes.txt"
  val DuplexFamilySizeMetricsExt: String = ".duplex_family_sizes.txt"
  val UmiMetricsExt: String              = ".umi_counts.txt"
  val YieldMetricsExt: String            = ".duplex_yield_metrics.txt"
  val PlotsExt: String                   = ".duplex_qc.pdf"

  private val PlottingScript = "com/fulcrumgenomics/umi/CollectDuplexSeqMetrics.R"

  /** Contains an AB and BA count of reads. */
  private case class Pair(ab: Int, ba: Int)

  /**
    * Metrics that captures information about family sizes for three kinds of families:
    *   1. CS or "Coordinate & Strand": families of reads that are grouped together by their 5' genomic
    *      positions and strands just as they are in traditional PCR duplicate marking
    *   2. SS or "Single Strand": single-strand families that are each subsets of a CS family create by
    *      also using the UMIs to partition the larger family, but not linking up families that are
    *      created from opposing strands of the same source molecule.
    *   3. DS or "Double Strand": families that are created by combining single-strand families that are from
    *      opposite strands of the same source molecule. This does NOT imply that all DS families are composed
    *      of reads from both strands; where only one strand of a source molecule is observed a DS family is
    *      still created.
    *
    * @param family_size the family size, i.e. the number of read pairs grouped together into a family
    * @param cs_count the count of families, of family size, when grouping just by coordinates and strand information
    * @param cs_fraction the fraction of all CS families where size == family_size
    * @param cs_fraction_gt_or_eq_size the fraction of all CS families where size >= family_size
    * @param ss_count the count of families, of family size, when also grouping by UMI to create single-strand families
    * @param ss_fraction the fraction of all SS families where size == family_size
    * @param ss_fraction_gt_or_eq_size the fraction of all SS families where size >= family_size
    * @param ds_count the count of families, of family size, when also grouping by UMI and merging single-strand
    *                 families from opposite strands of the same source molecule
    * @param ds_fraction the fraction of all DS families where size == family_size
    * @param ds_fraction_gt_or_eq_size the fraction of all DS families where size >= family_size
    */
  case class FamilySizeMetric(family_size: Int,
                              var cs_count: Long = 0,
                              var cs_fraction: Double = 0,
                              var cs_fraction_gt_or_eq_size: Double = 0,
                              var ss_count: Long = 0,
                              var ss_fraction: Double = 0,
                              var ss_fraction_gt_or_eq_size: Double = 0,
                              var ds_count: Long = 0,
                              var ds_fraction: Double = 0,
                              var ds_fraction_gt_or_eq_size: Double =0
                             ) extends Metric

  /**
    * Metrics that capture information specifically about double-stranded tag families. Since double-stranded
    * families have contributions from two different strands it is useful to see the distribution of families
    * across both strands.
    *
    * We refer to the two strands as "ab" and "ba" because we identify the two strands by observing the same pair of
    * UMIs (A and B) in opposite order (A->B vs B->A). Which strand is AB and which is BA is largely arbitrary, so
    * to make interpretation of the metrics simpler we use a definition here that for a given tag family
    * AB is the sub-family with more reads and BA is the tag family with fewer reads.
    *
    * @param ab_size The number of reads in the larger single-strand tag family for this double-strand tag family
    * @param ba_size The number of reads in the smaller single-strand tag family for this double-strand tag family
    * @param count The number of families with the A and B ss families of size ab_size and ba_size
    * @param fraction The fraction of all double-stranded tag families that have ab_size and ba_size
    * @param fraction_gt_or_eq_size The fraction of all double-stranded tag families that have
    *                               AB reads >= ab_size and BA reads >= ba_size
    */
  case class DuplexFamilySizeMetric(ab_size: Int,
                                    ba_size: Int,
                                    count: Long = 0,
                                    var fraction: Double = 0,
                                    var fraction_gt_or_eq_size: Double = 0
                                   ) extends Metric with Ordered[DuplexFamilySizeMetric] {

    override def compare(that: DuplexFamilySizeMetric): Int = {
      var retval = this.ab_size - that.ab_size
      if (retval == 0) retval = this.ba_size - that.ba_size
      retval
    }
  }

  /**
    * Metrics that are sampled at various levels of coverage, via random downsampling, during the construction
    * of duplex metrics.  The downsampling is done in such a way that the fractions are approximate, and not
    * exact, therefore the fraction field should only be interpreted as a guide and the read_pairs field used
    * to quantify how much data was used.
    *
    * @param fraction the approximate fraction of the full dataset that was used to generate the remaining values
    * @param read_pairs the number of read pairs upon which the remaining metrics are based
    * @param cs_families the number of CS (Coordinate & Strand) families present in the data
    * @param ss_families the number of SS (Single-Strand by UMI) families present in the data
    * @param ds_families the number of DS (Double-Strand by UMI) families present in the data
    * @param ds_duplexes the number of DS families that had the minimum number of observations on both strands to be
    *                    called duplexes (default = 1 read on each strand)
    * @param ds_fraction_duplexes the fraction of DS families that are duplexes (ds_duplexes / ds_families)
    * @param ds_fraction_duplexes_ideal the fraction of DS families that should be duplexes under an idealized model
    *                                   where each strand, A and B, have equal probability of being sampled, given
    *                                   the observed distribution of DS family sizes.
    */
  case class DuplexYieldMetric(fraction: Double,
                               read_pairs: Long,
                               cs_families: Long,
                               ss_families: Long,
                               ds_families: Long,
                               ds_duplexes: Long,
                               ds_fraction_duplexes: Double,
                               ds_fraction_duplexes_ideal: Double) extends Metric

  /**
    * Metrics describing the set of observed UMI sequences and the frequency of their observations.  The UMI
    * sequences reported may have been corrected using information within a double-stranded tag family.  For
    * example if a tag family is comprised of three read pairs with UMIs ACGT-TGGT ACGT-TGGT ACGT-TGGG that
    * a consensus UMI of will be generated ACGT-TGGT will be generated, and 3 raw observations counted for each
    * of ACGT and TGGT, and no observations counted for TGGG.
    *
    * @param umi the possibly-corrected UMI sequence
    * @param raw_observations the number of read pairs in the input BAM that observe the UMI (after correction)
    * @param raw_observations_with_errors the subset of raw-observations that underwent any correction
    * @param unique_observations the number of double-stranded tag families (i.e unique double-stranded molecules)
    *                            that observed the UMI
    * @param fraction_raw_observations the fraction of all raw observations that the UMI accounts for
    * @param fraction_unique_observations the fraction of all unique observations that the UMI accounts for
    */
  case class UmiMetric(umi: String,
                       var raw_observations: Long = 0,
                       var raw_observations_with_errors: Long = 0,
                       var unique_observations: Long = 0,
                       var fraction_raw_observations: Double = 0,
                       var fraction_unique_observations: Double = 0
                      ) extends Metric
}


@clp(group=ClpGroups.Umi, description=
  """
    |Collects a suite of metrics to QC duplex sequencing data.
    |
    |## Inputs
    |
    |The input to this tool must be a BAM file that is either:
    |
    |1. The exact BAM output by the `GroupReadsByUmi` tool (in the sort-order it was produced in)
    |2. A BAM file that has MI tags present on all reads (usually set by `GroupReadsByUmi` and has
    |   been sorted with `SortBam` into `TemplateCoordinate` order.
    |
    |Calculation of metrics may be restricted to a set of regions using the `--intervals` parameter. This
    |can significantly affect results as off-target reads in duplex sequencing experiments often have very
    |different properties than on-target reads due to the lack of enrichment.
    |
    |Several metrics are calculated related to the fraction of tag families that have duplex coverage. The
    |definition of "duplex" is controlled by the `--min-ab-reads` and `--min-ba-reads` parameters. The default
    |is to treat any tag family with at least one observation of each strand as a duplex, but this could be
    |made more stringent, e.g. by setting `--min-ab-reads=3 --min-ba-reads=3`.  If different thresholds are
    |used then `--min-ab-reads` must be the higher value.
    |
    |## Outputs
    |
    |The following output files produced:
    |
    |1. **<output>.umi_counts.txt**: metrics on the frequency of observations of UMIs within reads and tag families
    |2. **<output>.family_sizes.txt**: metrics on the frequency of different types of families of different sizes
    |3. **<output>.duplex_family_sizes.txt**: metrics on the frequency of duplex tag families by the number of
    |                                         observations from each strand
    |4. **<output>.duplex_yield_metrics.txt**: summary QC metrics produced using 5%, 10%, 15%...100% of the data
    |5. **<output>.duplx_qc.pdf**: a series of plots generated from the preceding metrics files for visualization
    |
    |Within the metrics files the prefixes `CS`, `SS` and `DS` are used to mean:
    |
    |* **CS**: tag families where membership is defined solely on matching genome coordinates and strand
    |* **SS**: single-stranded tag families where membership is defined by genome coordinates, strand and UMI;
    |          ie. 50/A and 50/B are considered different tag families.
    |* **DS**: double-stranded tag families where membership is collapsed across single-stranded tag families
    |          from the same double-stranded source molecule; i.e. 50/A and 50/B become one family
    |
    |## Requirements
    |
    |For plots to be generated R must be installed and the ggplot2 package installed with suggested
    |dependencies. Successfully executing the following in R will ensure a working installation:
    |
    |```R
    |install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)
    |```
  """)
class CollectDuplexSeqMetrics
( @arg(flag='i', doc="Input BAM file generated by GroupReadByUmi.") val input: PathToBam,
  @arg(flag='o', doc="Prefix of output files to write.") val output: PathPrefix,
  @arg(flag='l', doc="Optional set of intervals over which to restrict analysis.") val intervals: Option[PathToIntervals] = None,
  @arg(flag='d', doc="Description of dataset used to label plots. Defaults to sample/library.") val description: Option[String] = None,
  @arg(flag='a', doc="Minimum AB reads to call a tag family a 'duplex'.") val minAbReads: Int = 1,
  @arg(flag='b', doc="Minimum BA reads to call a tag family a 'duplex'.") val minBaReads: Int = 1,
  val generatePlots: Boolean = true // not a CLP arg - here to allow disabling of plots to speed up testing
) extends FgBioTool with LazyLogging {
  import CollectDuplexSeqMetrics._

  // Validate inputs
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertReadable)
  validate(minAbReads >= 1, "min-ab-reads must be >= 1")
  validate(minBaReads >= 0, "min-ba-reads must be >= 0")
  validate(minBaReads <= minAbReads, "min-ab-reads must be >= min-ba-reads")

  // Setup a whole bunch of counters for various things!
  private val dsLevels               = Range.inclusive(1, 20).toArray.map(_ * 0.05)
  private val startStopFamilyCounter = dsLevels.map(f => f -> new NumericCounter[Int]).toMap
  private val duplexFamilyCounter    = dsLevels.map(f => f -> new SimpleCounter[Pair]).toMap
  private val umiMetricsMap          = mutable.Map[String,UmiMetric]()

  // A consensus caller used to generate consensus UMI sequences
  private val consensusBuilder = new ConsensusCaller(errorRatePreLabeling=90.toByte, errorRatePostLabeling=90.toByte).builder()

  // A Murmur3 hash used to do the downsampling of the reads by generating an int hash of the read name
  // scaling it into the range 0-1 and then storing it a a transient attribute on the SAMRecord
  private val hasher         = new Murmur3(42)
  private val MaxIntAsDouble = Int.MaxValue.toDouble
  private val HashKey        = "_P".intern()

  override def execute(): Unit = {
    // Build the iterator we'll use based on whether or not we're restricting to a set of intervals
    val in = SamSource(input)
    val _filteredIterator = in.iterator.filter(r => r.paired && r.mapped && r.mateMapped && r.firstOfPair && !r.secondary && !r.supplementary)
    val iterator = intervals match {
      case None       => _filteredIterator
      case Some(path) =>
        val ilist    = IntervalList.fromFile(path.toFile).uniqued(false)
        val detector = new OverlapDetector[Interval](0,0)
        detector.addAll(ilist.getIntervals, ilist.getIntervals)
        _filteredIterator.filter { rec =>
          val (start, end) = if (rec.refIndex == rec.mateRefIndex) Bams.insertCoordinates(rec) else (rec.start, rec.end)
          detector.overlapsAny(new Interval(rec.refName, start, end))
        }
    }

    // Default the decription to something sensible if one wasn't provided
    val description = this.description.getOrElse {
      val samples   = in.readGroups.map(_.getSample).distinct
      val libraries = in.readGroups.map(_.getLibrary).distinct
      (samples, libraries) match {
        case (Seq(sample), Seq(library)) => s"$sample / $library"
        case _                           => input.getFileName.toString
      }
    }

    // Do a bunch of metrics collection
    collect(iterator)
    in.close()

    // Write the output files
    write(description)
  }

  /** Consumes from the iterator and collects information internal from which to generate metrics. */
  def collect(iterator: Iterator[SamRecord]): Unit = {
    val buffered = iterator.bufferBetter
    val progress = ProgressLogger(logger)

    while (buffered.hasNext) {
      val group = takeNextGroup(buffered)

      // Assign reads a random number between 0 and 1 inclusive based on their read name
      group.foreach { rec =>
        val intHash    = math.abs(this.hasher.hashUnencodedChars(rec.name))
        val doubleHash = intHash / MaxIntAsDouble
        rec.transientAttrs(HashKey) = doubleHash
      }

      // Update the counters
      this.dsLevels.foreach { fraction =>
        val downsampledGroup = group.filter(_.transientAttrs[Double](HashKey) <= fraction)

        if (downsampledGroup.nonEmpty) {
          this.startStopFamilyCounter(fraction).count(downsampledGroup.size)
          val dsGroups = downsampledGroup.groupBy(r => r[String](MI).takeWhile(_ != '/')).values.toSeq
          dsGroups.foreach { dsGroup =>
            // Family sizes first
            val ssGroups = dsGroup.groupBy(r => r[String](MI)).values.toSeq
            val counts = ssGroups.map(_.size).sortBy(x => -x)
            val ab = counts.head
            val ba = if (counts.length == 2) counts(1) else 0
            duplexFamilyCounter(fraction).count(Pair(ab, ba))

            // Then the UMIs
            if (fraction == 1.0) updateUmiMetrics(ssGroups)
          }
        }
      }

      // Record progress
      group.foreach(progress.record)
    }
  }

  /** Updates all the various metrics for the UMIs contained on the given records. */
  private def updateUmiMetrics(ssGroups: Seq[Seq[SamRecord]]): Unit = {
    val (ab, ba) = ssGroups match {
      case Seq(a, b) => (a, b)
      case Seq(a)    => (a, Seq.empty)
      case _         => unreachable(s"Found ${ssGroups.size} single strand families in a double strand family!")
    }

    val umi1s = new ArrayBuffer[String]
    val umi2s = new ArrayBuffer[String]

    ab.iterator.map(r => r[String](RX).split('-')).foreach { case Array(u1, u2) => umi1s += u1; umi2s += u2 }
    ba.iterator.map(r => r[String](RX).split('-')).foreach { case Array(u1, u2) => umi1s += u2; umi2s += u1 }

    Seq(umi1s, umi2s).foreach { umis =>
      val consensus = callConsensus(umis)
      if (consensus.contains('N')) {
        val x = 0
      }

      val metric    = this.umiMetricsMap.getOrElseUpdate(consensus, UmiMetric(umi=consensus))
      metric.raw_observations    += umis.size
      metric.unique_observations += 1
      metric.raw_observations_with_errors += umis.filterNot(_ == consensus).size
    }
  }

  /** Calls a simple consensus sequences from a set of sequences all the same length. */
  private def callConsensus(umis: Seq[String]) = {
    require(umis.nonEmpty, "Can't call consensus on an empty set of UMIs!")
    val buffer = new StringBuilder

    forloop(from=0, until=umis.head.length) { i =>
      this.consensusBuilder.reset()
      // The actual values for pError and pTruth don't matter since they're constant, but they are ln(0.01) and ln(0.99)
      umis.foreach(umi => this.consensusBuilder.add(umi.charAt(i).toByte, pError = -4.60517, pTruth = -0.01005034))
      val (base, qual) = this.consensusBuilder.call()
      buffer.append(base.toChar)
    }

    buffer.toString()
  }

  /** Generates the family size metrics from the current observations. */
  def familySizeMetrics: Seq[FamilySizeMetric] = {
    val map = mutable.Map[Int, FamilySizeMetric]()
    val startStopCounter = this.startStopFamilyCounter(1.0)
    val duplexCounter    = this.duplexFamilyCounter(1.0)

    // Add information about families grouped by genomic location alone
    startStopCounter.foreach { case (size, count) =>
      map.getOrElseUpdate(size, FamilySizeMetric(family_size=size)).cs_count += count
    }

    // Add information about the families grouped by ss and ds families
    duplexCounter.foreach { case (Pair(ab, ba), count) =>
      map.getOrElseUpdate(ab, FamilySizeMetric(family_size=ab)).ss_count += count
      if (ba > 0) map.getOrElseUpdate(ba, FamilySizeMetric(family_size=ba)).ss_count += count
      map.getOrElseUpdate(ab+ba, FamilySizeMetric(family_size=ab+ba)).ds_count += count
    }

    // Get a sorted Seq of metrics
    val metrics = map.values.toIndexedSeq.sortBy(_.family_size)

    // Fill in the fractions and cumulative fractions
    val csTotal = metrics.map(_.cs_count).sum.toDouble
    val ssTotal = metrics.map(_.ss_count).sum.toDouble
    val dsTotal = metrics.map(_.ds_count).sum.toDouble

    var csInverseCumulativeFraction = 0.0
    var ssInverseCumulativeFraction = 0.0
    var dsInverseCumulativeFraction = 0.0
    metrics.foreach { m =>
      m.cs_fraction = m.cs_count / csTotal
      m.ss_fraction = m.ss_count / ssTotal
      m.ds_fraction = m.ds_count / dsTotal

      m.cs_fraction_gt_or_eq_size = 1 - csInverseCumulativeFraction
      m.ss_fraction_gt_or_eq_size = 1 - ssInverseCumulativeFraction
      m.ds_fraction_gt_or_eq_size = 1 - dsInverseCumulativeFraction

      csInverseCumulativeFraction += m.cs_fraction
      ssInverseCumulativeFraction += m.ss_fraction
      dsInverseCumulativeFraction += m.ds_fraction
    }

    metrics
  }

  /** Generates the duplex family size metrics from the current observations. */
  def duplexFamilySizeMetrics: Seq[DuplexFamilySizeMetric] = {
    val metrics = this.duplexFamilyCounter(1.0)
      .map { case (Pair(ab, ba), count) => DuplexFamilySizeMetric(ab_size=ab, ba_size=ba, count=count) }
      .toIndexedSeq.sorted

    // Set the fractions
    val total = metrics.map(_.count).sum.toDouble
    metrics.foreach(m => m.fraction = m.count / total)

    // Set the cumulative fractions - there's probably a smarter way to do this!
    metrics.foreach { m =>
      val countGtOrEq = metrics.iterator.filter(n => n.ab_size >= m.ab_size && n.ba_size >= m.ba_size).map(_.count).sum
      m.fraction_gt_or_eq_size = countGtOrEq / total
    }

    metrics
  }

  /** Generates the duplex yield metrics from the current observations. */
  def yieldMetrics: Seq[DuplexYieldMetric] = {
    this.dsLevels.sorted.map { fraction =>
      val startStopCounter = this.startStopFamilyCounter(fraction)
      val duplexCounter    = this.duplexFamilyCounter(fraction)

      val dsFamilies       = duplexCounter.map { case (Pair(a,b), count) => count }.sum
      val countOfDuplexes  = duplexCounter.map {
        case (Pair(a,b), count) if a >= this.minAbReads && b >= this.minBaReads => count
        case _                                                                  => 0
      }.sum

      val countOfDuplexesIdeal = duplexCounter.map { case (pair, count) => count * pDuplexIdeal(pair.ab + pair.ba) }.sum


      new DuplexYieldMetric(
        fraction                   = fraction,
        read_pairs                 = startStopCounter.totalMass.toLong,
        cs_families                = startStopCounter.total,
        ss_families                = duplexCounter.map { case (Pair(a,b), count) => if (b>0) count*2 else count }.sum,
        ds_families                = dsFamilies,
        ds_duplexes                = countOfDuplexes,
        ds_fraction_duplexes       = countOfDuplexes / dsFamilies.toDouble,
        ds_fraction_duplexes_ideal = countOfDuplexesIdeal / dsFamilies.toDouble
      )
    }
  }

  /** Generates the UMI metrics from the current observations. */
  def umiMetrics: Seq[UmiMetric] = {
    val metrics     = this.umiMetricsMap.values.toIndexedSeq.sortBy(_.umi)
    val rawTotal    = metrics.map(_.raw_observations).sum.toDouble
    val uniqueTotal = metrics.map(_.unique_observations).sum.toDouble

    metrics.foreach { m =>
      m.fraction_raw_observations    = m.raw_observations / rawTotal
      m.fraction_unique_observations = m.unique_observations / uniqueTotal
    }

    metrics
  }

  /**
    * For a given family size/number of reads computes the probability that we would observed
    * at least minAb & minBa reads under a binomial sampling model with p=0.5.
    */
  def pDuplexIdeal(reads: Int): Double = {
    if (reads < this.minAbReads + this.minBaReads) {
      0
    }
    else {
      val minSsReads = math.min(this.minAbReads, this.minBaReads)

      // Here we can assert:
      //   a) reads >= minAbReads + minBaReads
      //   b) minSsReads <= minAbReads
      //   c) minSsReads <= minBaReads
      //   c) reads - minSsReads >= max(minAbReads, minBaReads)
      //   d) p(duplex == 1 | successes in [minSsReads..reads-minSsReads]
      val binom = new BinomialDistribution(reads, 0.5)
      Range.inclusive(minSsReads, reads - minSsReads)
        .map(successes => binom.probability(successes))
        .sum
    }
  }

  /**
    * Grabs the next group of records that all share the same start/stop/strand information. This can
    * and will contain reads with different MIs!
    */
  private def takeNextGroup(iterator: BetterBufferedIterator[SamRecord]): Seq[SamRecord] = {
    val rec = iterator.head
    val info = GroupReadsByUmi.ReadInfo(rec)
    iterator.takeWhile(rec => GroupReadsByUmi.ReadInfo(rec) == info).toIndexedSeq
  }

  /** Writes out all the metrics and plots. */
  private[umi] def write(description: String): Unit = {
    val Seq(fsPath, dfsPath, umiPath, yieldPath, pdfPath) = Seq(FamilySizeMetricsExt, DuplexFamilySizeMetricsExt, UmiMetricsExt, YieldMetricsExt, PlotsExt).map { ext =>
      output.getParent.resolve(output.getFileName + ext)
    }

    Metric.write(fsPath, familySizeMetrics)
    Metric.write(dfsPath, duplexFamilySizeMetrics)
    Metric.write(yieldPath, yieldMetrics)
    Metric.write(umiPath, umiMetrics)

    if (generatePlots) {
      Rscript.execIfAvailable(PlottingScript, fsPath.toString, dfsPath.toString, yieldPath.toString, umiPath.toString, pdfPath.toString, description) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => Unit
      }
    }
  }
}
