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
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
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
  val DuplexUmiMetricsExt: String        = ".duplex_umi_counts.txt"
  val YieldMetricsExt: String            = ".duplex_yield_metrics.txt"
  val PlotsExt: String                   = ".duplex_qc.pdf"

  private val PlottingScript = "com/fulcrumgenomics/umi/CollectDuplexSeqMetrics.R"

  /** Contains an AB and BA count of reads. */
  private case class Pair(ab: Int, ba: Int)

  /**
    * Metrics produced by `CollectDuplexSeqMetrics` to quantify the distribution of different kinds of read family
    * sizes.  Three kinds of families are described:
    *
    * 1. _CS_ or _Coordinate & Strand_: families of reads that are grouped together by their unclipped 5'
    *    genomic positions and strands just as they are in traditional PCR duplicate marking
    * 2. _SS_ or _Single Strand_: single-strand families that are each subsets of a CS family create by
    *    also using the UMIs to partition the larger family, but not linking up families that are
    *    created from opposing strands of the same source molecule.
    * 3. _DS_ or _Double Strand_: families that are created by combining single-strand families that are from
    *    opposite strands of the same source molecule. This does **not** imply that all DS families are composed
    *    of reads from both strands; where only one strand of a source molecule is observed a DS family is
    *    still counted.
    *
    * @param family_size The family size, i.e. the number of read pairs grouped together into a family.
    * @param cs_count The count of families with `size == family_size` when grouping just by coordinates and strand information.
    * @param cs_fraction The fraction of all _CS_ families where `size == family_size`.
    * @param cs_fraction_gt_or_eq_size The fraction of all _CS_ families where `size >= family_size`.
    * @param ss_count The count of families with `size == family_size` when also grouping by UMI to create single-strand families.
    * @param ss_fraction The fraction of all _SS_ families where `size == family_size`.
    * @param ss_fraction_gt_or_eq_size The fraction of all _SS_ families where `size >= family_size`.
    * @param ds_count The count of families with `size == family_size`when also grouping by UMI and merging single-strand
    *                 families from opposite strands of the same source molecule.
    * @param ds_fraction The fraction of all _DS_ families where `size == family_size`.
    * @param ds_fraction_gt_or_eq_size The fraction of all _DS_ families where `size >= family_size`.
    */
  case class FamilySizeMetric(family_size: Int,
                              var cs_count: Count = 0,
                              var cs_fraction: Proportion = 0,
                              var cs_fraction_gt_or_eq_size: Proportion = 0,
                              var ss_count: Count = 0,
                              var ss_fraction: Proportion = 0,
                              var ss_fraction_gt_or_eq_size: Proportion = 0,
                              var ds_count: Count = 0,
                              var ds_fraction: Proportion = 0,
                              var ds_fraction_gt_or_eq_size: Proportion =0
                             ) extends Metric

  /**
    * Metrics produced by `CollectDuplexSeqMetrics` to describe the distribution of double-stranded (duplex)
    * tag families in terms of the number of reads observed on each strand.
    *
    * We refer to the two strands as `ab` and `ba` because we identify the two strands by observing the same pair of
    * UMIs (A and B) in opposite order (A->B vs B->A). Which strand is `ab` and which is `ba` is largely arbitrary, so
    * to make interpretation of the metrics simpler we use a definition here that for a given tag family
    * `ab` is the sub-family with more reads and `ba` is the tag family with fewer reads.
    *
    * @param ab_size The number of reads in the `ab` sub-family (the larger sub-family) for this double-strand tag family.
    * @param ba_size The number of reads in the `ba` sub-family (the smaller sub-family) for this double-strand tag family.
    * @param count The number of families with the `ab` and `ba` single-strand families of size `ab_size` and `ba_size`.
    * @param fraction The fraction of all double-stranded tag families that have `ab_size` and `ba_size`.
    * @param fraction_gt_or_eq_size The fraction of all double-stranded tag families that have
    *                               `ab reads >= ab_size` and `ba reads >= ba_size`.
    */
  case class DuplexFamilySizeMetric(ab_size: Int,
                                    ba_size: Int,
                                    count: Count = 0,
                                    var fraction: Proportion = 0,
                                    var fraction_gt_or_eq_size: Proportion = 0
                                   ) extends Metric with Ordered[DuplexFamilySizeMetric] {

    /** Orders by ab_size, then ba_size. */
    override def compare(that: DuplexFamilySizeMetric): Int = {
      var retval = this.ab_size - that.ab_size
      if (retval == 0) retval = this.ba_size - that.ba_size
      retval
    }
  }

  /**
    * Metrics produced by `CollectDuplexSeqMetrics` that are sampled at various levels of coverage, via random
    * downsampling, during the construction of duplex metrics.  The downsampling is done in such a way that the
    * `fraction`s are approximate, and not exact, therefore the `fraction` field should only be interpreted as a guide
    * and the `read_pairs` field used to quantify how much data was used.
    *
    * See `FamilySizeMetric` for detailed definitions of `CS`, `SS` and `DS` as used below.
    *
    * @param fraction    The approximate fraction of the full dataset that was used to generate the remaining values.
    * @param read_pairs  The number of read pairs upon which the remaining metrics are based.
    * @param cs_families The number of _CS_ (Coordinate & Strand) families present in the data.
    * @param ss_families The number of _SS_ (Single-Strand by UMI) families present in the data.
    * @param ds_families The number of _DS_ (Double-Strand by UMI) families present in the data.
    * @param ds_duplexes The number of _DS_ families that had the minimum number of observations on both strands to be
    *                    called duplexes (default = 1 read on each strand).
    * @param ds_fraction_duplexes The fraction of _DS_ families that are duplexes (`ds_duplexes / ds_families`).
    * @param ds_fraction_duplexes_ideal The fraction of _DS_ families that should be duplexes under an idealized model
    *                                   where each strand, `A` and `B`, have equal probability of being sampled, given
    *                                   the observed distribution of _DS_ family sizes.
    */
  case class DuplexYieldMetric(fraction: Proportion,
                               read_pairs: Count,
                               cs_families: Count,
                               ss_families: Count,
                               ds_families: Count,
                               ds_duplexes: Count,
                               ds_fraction_duplexes: Proportion,
                               ds_fraction_duplexes_ideal: Proportion) extends Metric

  /**
    * Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed UMI sequences and the
    * frequency of their observations.  The UMI sequences reported may have been corrected using information
    * within a double-stranded tag family.  For example if a tag family is comprised of three read pairs with
    * UMIs `ACGT-TGGT`, `ACGT-TGGT`, and `ACGT-TGGG` then a consensus UMI of `ACGT-TGGT` will be generated,
    * and three raw observations counted for each of `ACGT` and `TGGT`, and no observations counted for `TGGG`.
    *
    * @param umi The UMI sequence, possibly-corrected.
    * @param raw_observations The number of read pairs in the input BAM that observe the UMI (after correction).
    * @param raw_observations_with_errors The subset of raw-observations that underwent any correction.
    * @param unique_observations The number of double-stranded tag families (i.e unique double-stranded molecules)
    *                            that observed the UMI.
    * @param fraction_raw_observations The fraction of all raw observations that the UMI accounts for.
    * @param fraction_unique_observations The fraction of all unique observations that the UMI accounts for.
    */
  case class UmiMetric(umi: String,
                       var raw_observations: Count = 0,
                       var raw_observations_with_errors: Count = 0,
                       var unique_observations: Count = 0,
                       var fraction_raw_observations: Proportion = 0,
                       var fraction_unique_observations: Proportion = 0
                      ) extends Metric

  /**
    * Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed duplex UMI sequences and the
    * frequency of their observations.  The UMI sequences reported may have been corrected using information
    * within a double-stranded tag family.  For example if a tag family is comprised of three read pairs with
    * UMIs `ACGT-TGGT`, `ACGT-TGGT`, and `ACGT-TGGG` then a consensus UMI of `ACGT-TGGT` will be generated.
    *
    * UMI pairs are normalized within a tag family so that observations are always reported as if they came
    * from a read pair with read 1 on the positive strand (F1R2). Another way to view this is that for FR or RF
    * read pairs, the duplex UMI reported is the UMI from the positive strand read followed by the UMI from the
    * negative strand read.  E.g. a read pair with UMI `AAAA-GGGG` and with R1 on the negative strand and R2 on
    * the positive strand, will be reported as `GGGG-AAAA`.
    *
    * @param umi The duplex UMI sequence, possibly-corrected.
    * @param raw_observations The number of read pairs in the input BAM that observe the duplex UMI (after correction).
    * @param raw_observations_with_errors The subset of raw observations that underwent any correction.
    * @param unique_observations The number of double-stranded tag families (i.e unique double-stranded molecules)
    *                            that observed the duplex UMI.
    * @param fraction_raw_observations The fraction of all raw observations that the duplex UMI accounts for.
    * @param fraction_unique_observations The fraction of all unique observations that the duplex UMI accounts for.
    * @param fraction_unique_observations_expected The fraction of all unique observations that are expected to be
    *                                              attributed to the duplex UMI based on the `fraction_unique_observations`
    *                                              of the two individual UMIs.
    */
  case class DuplexUmiMetric(umi: String,
                             var raw_observations: Count = 0,
                             var raw_observations_with_errors: Count = 0,
                             var unique_observations: Count = 0,
                             var fraction_raw_observations: Proportion = 0,
                             var fraction_unique_observations: Proportion = 0,
                             var fraction_unique_observations_expected: Proportion = 0
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
    |The following output files are produced:
    |
    |1. **<output>.family_sizes.txt**: metrics on the frequency of different types of families of different sizes
    |2. **<output>.duplex_family_sizes.txt**: metrics on the frequency of duplex tag families by the number of
    |                                        observations from each strand
    |3. **<output>.duplex_yield_metrics.txt**: summary QC metrics produced using 5%, 10%, 15%...100% of the data
    |4. **<output>.umi_counts.txt**: metrics on the frequency of observations of UMIs within reads and tag families
    |5. **<output>.duplex_qc.pdf**: a series of plots generated from the preceding metrics files for visualization
    |6. **<output>.duplex_umi_counts.txt**: (optional) metrics on the frequency of observations of duplex UMIs within
    |   reads and tag families. This file is only produced _if_ the `--duplex-umi-counts` option is used as it
    |   requires significantly more memory to track all pairs of UMIs seen when a large number of UMI sequences are present.
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
( @arg(flag='i', doc="Input BAM file generated by `GroupReadByUmi`.") val input: PathToBam,
  @arg(flag='o', doc="Prefix of output files to write.") val output: PathPrefix,
  @arg(flag='l', doc="Optional set of intervals over which to restrict analysis.") val intervals: Option[PathToIntervals] = None,
  @arg(flag='d', doc="Description of data set used to label plots. Defaults to sample/library.") val description: Option[String] = None,
  @arg(flag='u', doc="If true, produce the .duplex_umi_counts.txt file with counts of duplex UMI observations.") val duplexUmiCounts: Boolean = false,
  @arg(flag='a', doc="Minimum AB reads to call a tag family a 'duplex'.") val minAbReads: Int = 1,
  @arg(flag='b', doc="Minimum BA reads to call a tag family a 'duplex'.") val minBaReads: Int = 1,
  private val generatePlots: Boolean = true // not a CLP arg - here to allow disabling of plots to speed up testing
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
  private val duplexUmiMetricsMap    = mutable.Map[String,DuplexUmiMetric]()

  // A consensus caller used to generate consensus UMI sequences
  private val consensusBuilder = new SimpleConsensusCaller()

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
    val description = this.description.getOrElse { plotDescription(in, input) }

    // Do a bunch of metrics collection
    collect(iterator)
    in.safelyClose()

    // Write the output files
    write(description)
  }

  /** Consumes from the iterator and collects information internally from which to generate metrics. */
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

  /**
    * Updates all the various metrics for the UMIs contained on the given records.
    *
    * @param ssGroups a Seq of length one or two containing the one or two single-stranded tag families
    *                 that comprise an individual double-stranded tag family.  If there are two single-stranded
    *                 tag families, there is no guarantee of their ordering within the Seq.
    */
  private[umi] def updateUmiMetrics(ssGroups: Seq[Seq[SamRecord]]): Unit = {
    // ab and ba are just names here and _don't_ imply top/bottom strand, just different strands
    val (ab, ba) = ssGroups match {
      case Seq(a, b) => (a, b)
      case Seq(a)    => (a, Seq.empty)
      case _         => unreachable(s"Found ${ssGroups.size} single strand families in a double strand family!")
    }

    val umi1s = new ArrayBuffer[String] // ab UMI 1 (and ba UMI 2) sequences
    val umi2s = new ArrayBuffer[String] // ab UMI 2 (and ba UMI 1) sequences

    ab.iterator.map(r => r[String](RX).split('-')).foreach { case Array(u1, u2) => umi1s += u1; umi2s += u2 }
    ba.iterator.map(r => r[String](RX).split('-')).foreach { case Array(u1, u2) => umi1s += u2; umi2s += u1 }

    val Seq(abConsensusUmi, baConsensusUmi) = Seq(umi1s, umi2s).map{ umis =>
      val consensus = this.consensusBuilder.callConsensus(umis)
      val metric    = this.umiMetricsMap.getOrElseUpdate(consensus, UmiMetric(umi=consensus))
      metric.raw_observations    += umis.size
      metric.unique_observations += 1
      metric.raw_observations_with_errors += umis.filterNot(_ == consensus).size
      consensus
    }

    if (this.duplexUmiCounts) {
      // Make both possible consensus duplex UMIs.  We want to normalize to the pairing/orientation that we'd
      // see on an F1R2 read-pair. We can get this either by direct observation in the `ab` set of reads,
      // or by seeing an F2R1 in the `ba` set of reads.  This logic will pick more or less arbitrarily if the
      // group of reads doesn't consist of the expected all F1R2s in one set and all F2R1s in the other set.
      val duplexUmis = Seq(s"$abConsensusUmi-$baConsensusUmi", s"$baConsensusUmi-$abConsensusUmi")
      val duplexUmi  = if (ab.exists(r => r.firstOfPair && r.positiveStrand) || ba.exists(r => r.secondOfPair && r.positiveStrand)) duplexUmis(0) else duplexUmis(1)

      val metric = this.duplexUmiMetricsMap.getOrElseUpdate(duplexUmi, DuplexUmiMetric(umi=duplexUmi))
      metric.raw_observations             += ab.size + ba.size
      metric.unique_observations          += 1
      metric.raw_observations_with_errors += (ab.iterator ++ ba.iterator).map(r => r[String](RX)).count(umi => !duplexUmis.contains(umi))
    }
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

      DuplexYieldMetric(
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

  /** Generates the duplex UMI metrics from the current observations. */
  def duplexUmiMetrics(umiMetrics: Seq[UmiMetric]): Seq[DuplexUmiMetric] = {
    val singleUmiMetrics = umiMetrics.map(m => m.umi -> m).toMap
    val metrics          = this.duplexUmiMetricsMap.values.toIndexedSeq.sortBy(x => - x.unique_observations)
    val rawTotal         = metrics.map(_.raw_observations).sum.toDouble
    val uniqueTotal      = metrics.map(_.unique_observations).sum.toDouble

    metrics.foreach { m =>
      val Array(umi1, umi2)                   = m.umi.split('-')
      m.fraction_raw_observations             = m.raw_observations / rawTotal
      m.fraction_unique_observations          = m.unique_observations / uniqueTotal
      m.fraction_unique_observations_expected = singleUmiMetrics(umi1).fraction_unique_observations * singleUmiMetrics(umi2).fraction_unique_observations
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
    val Seq(fsPath, dfsPath, umiPath, duplexUmiPath, yieldPath, pdfPath) =
      Seq(FamilySizeMetricsExt, DuplexFamilySizeMetricsExt, UmiMetricsExt, DuplexUmiMetricsExt, YieldMetricsExt, PlotsExt).map { ext =>
        output.getParent.resolve(output.getFileName + ext)
      }

    Metric.write(fsPath, familySizeMetrics)
    Metric.write(dfsPath, duplexFamilySizeMetrics)
    Metric.write(yieldPath, yieldMetrics)

    val umiMetricValues = umiMetrics
    Metric.write(umiPath, umiMetricValues )
    if (this.duplexUmiCounts) Metric.write(duplexUmiPath, duplexUmiMetrics(umiMetricValues))

    if (generatePlots) {
      Rscript.execIfAvailable(PlottingScript, fsPath.toString, dfsPath.toString, yieldPath.toString, umiPath.toString, pdfPath.toString, description) match {
        case Failure(e) => logger.warning(s"Generation of PDF plots failed: ${e.getMessage}")
        case _ => Unit
      }
    }
  }
}
