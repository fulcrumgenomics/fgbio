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

import java.lang.Math.{max, min}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker

@clp(
  group = ClpGroups.Umi, description =
    """
      |Filters consensus reads generated by _CallMolecularConsensusReads_ or _CallDuplexConsensusReads_.
      |Two kinds of filtering are performed:
      |
      |  1. Masking/filtering of individual bases in reads
      |  2. Filtering out of reads (i.e. not writing them to the output file)
      |
      |Base-level filtering/masking is only applied if per-base tags are present (see _CallDuplexConsensusReads_ and
      |_CallMolecularConsensusReads_ for descriptions of these tags). Read-level filtering is always applied.  When
      |filtering reads, secondary alignments and supplementary records may be removed independently if they fail
      |one or more filters; if either R1 or R2 primary alignments fail a filter then all records for the template
      |will be filtered out.
      |
      |The filters applied are as follows:
      |
      |  1. Reads with fewer than min-reads contributing reads are filtered out
      |  2. Reads with an average consensus error rate higher than max-read-error-rate are filtered out
      |  3. Bases with quality scores below min-base-quality are masked to Ns
      |  4. Bases with fewer than min-reads contributing raw reads are masked to Ns
      |  5. Bases with a consensus error rate (defined as the fraction of contributing reads that
      |     voted for a different base than the consensus call) higher than max-base-error-rate are masked to Ns
      |  6. Reads with a proportion of Ns higher than max-no-call-fraction *after* per-base filtering are filtered out.
      |
      |When filtering _single-umi consensus_ reads generated by _CallMolecularConsensusReads_ a single value each
      |should be supplied for `--min-reads`, `--max-read-error-rate`, and `--max-base-error-rate`.
      |
      |When filtering duplex consensus reads generated by _CallDuplexConsensusReads_ each of the three parameters
      |may independently take 1-3 values. For example:
      |
      |```
      |FilterConsensusReads ... --min-reads 10 5 3 --max-base-error-rate 0.1
      |```
      |
      |In each case if fewer than three values are supplied, the last value is repeated (i.e. `80 40` -> `80 40 40`
      |and `0.1` -> `0.1 0.1 0.1`.  The first value applies to the final consensus read, the second value to one
      |single-strand consensus, and the last value to the other single-strand consensus. It is required that if
      |values two and three differ, the _more stringent value comes earlier_.
      |
      |In order to correctly filter reads in or out by template, if the input BAM is not `queryname` sorted or
      |grouped it will be sorted into queryname order.  The resulting records are `coordinate` sorted to efficiently
      |correct the `NM`, `UQ` and `MD` tags, and the output BAM is always written in coordinate order.
      |
      |The `--reverse-tags-per-base` option controls whether per-base tags should be reversed before being used on reads
      |marked as being mapped to the negative strand.  This is necessary if the reads have been mapped and the
      |bases/quals reversed but the consensus tags have not.  If true, the tags written to the output BAM will be
      |reversed where necessary in order to line up with the bases and quals.
    """
)
class FilterConsensusReads
( @arg(flag='i', doc="The input SAM or BAM file of consensus reads.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='r', doc="Reference fasta file.") val ref: PathToFasta,
  @arg(flag='R', doc="Reverse [complement] per base tags on reverse strand reads.") val reversePerBaseTags: Boolean = false,
  @arg(flag='M', minElements=1, maxElements=3, doc="The minimum number of reads supporting a consensus base/read.")
  val minReads: Seq[Int],
  @arg(flag='E', minElements=1, maxElements=3, doc="The maximum raw-read error rate across the entire consensus read.")
  val maxReadErrorRate: Seq[Double] = Seq(0.025),
  @arg(flag='e', minElements=1, maxElements=3, doc="The maximum error rate for a single consensus base.")
  val maxBaseErrorRate: Seq[Double] = Seq(0.1),
  @arg(flag='N', doc="Mask (make `N`) consensus bases with quality less than this threshold.")
  val minBaseQuality: PhredScore,
  @arg(flag='n', doc="Maximum fraction of no-calls in the read after filtering.")
  val maxNoCallFraction: Double = 0.2

) extends FgBioTool with LazyLogging {
  // Baseline input validation
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)
  if (maxReadErrorRate.exists(e => e < 0 || e > 1)) fail("max-read-error-rate must be between 0 and 1.")
  if (maxBaseErrorRate.exists(e => e < 0 || e > 1)) fail("max-base-error-rate must be between 0 and 1.")
  if (maxNoCallFraction < 0 || maxNoCallFraction > 1) fail("max-no-call-fraction must be between 0 and 1.")

  private val NoCall     = 'N'.toByte
  private val NoCallQual = PhredScore.MinValue
  private val MaxRecordsInMemoryWhenSorting = 250e3.toInt

  // Filtering objects that are used for the main read values, and if it's a duplex, for the two
  // single strand consensus reads.  Rules are as follows:
  //   - The first value is always used for the main read values (Vanilla or Duplex)
  //   - The second value is used for AB filtering if provided, otherwise the first value
  //   - The third value is used for BA filtering if provided, otherwise the second value, otherwise the first
  private[umi] case class Filters(minReads: Int, maxReadErrorRate: Double, maxBaseErrorRate: Double)
  private[umi] val ccFilters = Filters( // For the final Consensus Calls
    minReads         = minReads.head,
    maxReadErrorRate = maxReadErrorRate.head,
    maxBaseErrorRate = maxBaseErrorRate.head
  )
  private[umi] val abFilters = Filters(
    minReads         = minReads.take(2).last,
    maxReadErrorRate = maxReadErrorRate.take(2).last,
    maxBaseErrorRate = maxBaseErrorRate.take(2).last
  )
  private[umi] val baFilters = Filters(
    minReads         = minReads.last,
    maxReadErrorRate = maxReadErrorRate.last,
    maxBaseErrorRate = maxBaseErrorRate.last
  )

  // For depth thresholds it's required that ba <= ab <= cc
  validate(abFilters.minReads <= ccFilters.minReads, "min-reads values must be specified high to low.")
  validate(baFilters.minReads <= abFilters.minReads, "min-reads values must be specified high to low.")

  // For error rates only enforce that ab is more stringent than ba, because you might conceivably
  // want to say '--max-base-error-rate 0.1 0.05 0.2' to allow up to 10% but only if at least one read is < 5%
  validate(baFilters.maxReadErrorRate >= abFilters.maxReadErrorRate, "max-read-error-rate for AB must be <= than for BA")
  validate(baFilters.maxBaseErrorRate >= abFilters.maxBaseErrorRate, "max-read-error-rate for AB must be <= than for BA")

  // Variables for tracking how many reads meet which fate
  private var totalReads:  Long = 0
  private var keptReads:   Long = 0
  private var maskedBases: Long = 0
  private var totalBases:  Long = 0

  // Case class to track results of filtering a single read
  case class FilterResult(keepRead: Boolean, maskedBases: Int)
  private val EmptyFilterResult = FilterResult(keepRead=true, maskedBases=0)

  override def execute(): Unit = {
    val in        = SamReaderFactory.make().open(input)
    val header    = in.getFileHeader.clone()
    header.setSortOrder(SortOrder.coordinate)
    val sorter    = Bams.sortingCollection(SortOrder.coordinate, header, maxInMemory=MaxRecordsInMemoryWhenSorting)
    val out       = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(header, true, output.toFile, null)
    val progress1 = new ProgressLogger(logger, verb="Filtered and masked")

    // Go through the reads by template and do the filtering
    val templateIterator = Bams.templateIterator(in, maxInMemory=MaxRecordsInMemoryWhenSorting)
    logger.info("Filtering reads.")
    templateIterator.foreach { template =>
      val r1 = template.r1.getOrElse(throw new IllegalStateException(s"${template.name} had no R1."))
      if (r1.getReadPairedFlag) require(template.r2.isDefined, s"Paired read missing R2: ${template.name}")
      val primaryReadCount = if (template.r2.isDefined) 2 else 1
      totalReads += primaryReadCount

      // Reverse the tags on the reads if needs be
      if (reversePerBaseTags) template.allReads.foreach(reverseConsensusTagsIfNeeded)

      // Filter and mask
      val r1Result = filterRecord(r1)
      val r2Result = if (r1Result.keepRead && template.r2.isDefined) filterRecord(template.r2.get) else EmptyFilterResult

      if (r1Result.keepRead && r2Result.keepRead) {
        keptReads   += primaryReadCount
        totalBases  += r1.getReadLength + template.r2.map(_.getReadLength).getOrElse(0)
        maskedBases += r1Result.maskedBases + r2Result.maskedBases
        sorter.add(r1)
        progress1.record(r1)
        template.r2.foreach { r => sorter.add(r); progress1.record(r) }

        template.allSupplementaryAndSecondary.foreach { r =>
          val result = filterRecord(r)
          if (result.keepRead) {
            sorter.add(r)
            progress1.record(r)
          }
        }
      }
    }

    // Then iterate the reads in coordinate order and re-calculate key tags
    logger.info("Filtering complete; fixing tags and writing coordinate sorted reads.")
    val progress2 = new ProgressLogger(logger, verb="Wrote")
    val walker = new ReferenceSequenceFileWalker(ref.toFile)
    sorter.foreach { rec =>
      Bams.regenerateNmUqMdTags(rec, walker)
      out.addAlignment(rec)
      progress2.record(rec)
    }

    in.safelyClose()
    out.close()
    logger.info(f"Output ${keptReads}%,d of ${totalReads}%,d primary consensus reads.")
    logger.info(f"Masked ${maskedBases}%,d of ${totalBases}%,d bases in retained primary consensus reads.")
  }

  /**
    * Perform filtering on an individual read.
    *
    * @param  rec the read to filter
    * @return a tuple of (Boolean, Int) where the first value is true if the read should be kept
    *         and false if it should be filtered out, and the second value is the count of bases
    *         in the read that were masked to Ns as a result of per-base filtering
    */
  private[umi] def filterRecord(rec: SAMRecord): FilterResult = {
    // Checks that are common to vanilla and duplex reads
    val maxDepth = rec.getIntegerAttribute(ConsensusTags.PerRead.RawReadCount)
    val errRate  = rec.getFloatAttribute(ConsensusTags.PerRead.RawReadErrorRate)
    if (maxDepth == null || errRate == null) fail(s"Read ${rec.getReadName} does not appear to have consensus calling tags present.")

    // Only bother looking at reads where per-read criteria are met
    if (maxDepth < this.ccFilters.minReads || errRate > this.ccFilters.maxReadErrorRate) {
      FilterResult(false, 0)
    }
    else {
      val result = if (isDuplexRecord(rec)) filterDuplexConsensusRead(rec) else filterVanillaConsensusRead(rec)
      result.copy(keepRead = result.keepRead && fractionNoCalls(rec) <= this.maxNoCallFraction)
    }
  }

  /** Computes the fraction of the read that are no-calls. */
  private def fractionNoCalls(rec: SAMRecord): Double = {
    val bases = rec.getReadBases
    var ns = 0
    forloop(from = 0, until = rec.getReadLength) { i => if (bases(i) == NoCall) ns += 1 }
    ns / bases.length.toDouble
  }

  /** Performs filtering that is specific to Vanilla (single-umi/non-duplex) consensus reads. */
  private def filterVanillaConsensusRead(rec: SAMRecord): FilterResult = {
    val bases  = rec.getReadBases
    val quals  = rec.getBaseQualities
    val depths = rec.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadCount)
    val errors = rec.getSignedShortArrayAttribute(ConsensusTags.PerBase.RawReadErrors)
    val pb = depths != null && errors != null
    val filters = this.ccFilters

    // Do the per-base masking
    var maskedBasesThisRead = 0
    forloop(from = 0, until = rec.getReadLength) { i =>
      if ((quals(i) < this.minBaseQuality) ||
        (pb && depths(i) < filters.minReads) ||
        (pb && errors(i) / depths(i).toDouble > filters.maxBaseErrorRate)) {
        bases(i) = NoCall
        quals(i) = NoCallQual
        maskedBasesThisRead += 1
      }
    }

    rec.setReadBases(bases)
    rec.setBaseQualities(quals)
    FilterResult(keepRead=true, maskedBases=maskedBasesThisRead)
  }

  /** Performs filtering that is specific to Duplex consensus reads. */
  private def filterDuplexConsensusRead(rec: SAMRecord): FilterResult = {
    val failsReadLevelChecks = {
      import ConsensusTags.PerRead._
      val Seq(baDepth, abDepth) = Seq(AbRawReadCount,     BaRawReadCount    ).map(rec.getIntegerAttribute).sorted
      val Seq(abError, baError) = Seq(AbRawReadErrorRate, BaRawReadErrorRate).map(rec.getFloatAttribute).sorted
      abDepth < abFilters.minReads || abError > abFilters.maxReadErrorRate ||
        baDepth < baFilters.minReads || baError > baFilters.maxReadErrorRate
    }

    if (failsReadLevelChecks) {
      FilterResult(keepRead=false, maskedBases=0)
    }
    else {
      import ConsensusTags.PerBase._
      val bases   = rec.getReadBases
      val quals   = rec.getBaseQualities
      val depths1 = rec.getSignedShortArrayAttribute(AbRawReadCount)
      val depths2 = rec.getSignedShortArrayAttribute(BaRawReadCount)
      val errors1 = rec.getSignedShortArrayAttribute(AbRawReadErrors)
      val errors2 = rec.getSignedShortArrayAttribute(BaRawReadErrors)

      var maskedBases = 0
      forloop(from=0, until=bases.length) { i =>
        if (bases(i) != NoCall) {
          val abDepth = max(depths1(i), depths2(i))
          val baDepth = min(depths1(i), depths2(i))
          val abError = min(errors1(i)/depths1(i).toDouble, errors2(i) / depths2(i).toDouble)
          val baError = max(errors1(i)/depths1(i).toDouble, errors2(i) / depths2(i).toDouble)
          val totalDepth = abDepth + baDepth
          val totalError = (errors1(i) + errors2(i)) / totalDepth.toDouble

          if (totalDepth < ccFilters.minReads || totalError > ccFilters.maxBaseErrorRate ||
              abDepth    < abFilters.minReads || abError    > abFilters.maxBaseErrorRate ||
              baDepth    < baFilters.minReads || baError    > baFilters.maxBaseErrorRate ||
              quals(i) < this.minBaseQuality) {
            bases(i) = NoCall
            quals(i) = PhredScore.MinValue
            maskedBases += 1
          }
        }
      }

      rec.setReadBases(bases)
      rec.setBaseQualities(quals)
      FilterResult(keepRead=true, maskedBases=maskedBases)
    }
  }

  /** Quick check to see if a record is a duplex record. */
  private def isDuplexRecord(rec: SAMRecord): Boolean = rec.getAttribute(ConsensusTags.PerRead.AbRawReadCount) != null

  /** Reverses all the consensus tags if the right conditions are set. */
  private def reverseConsensusTagsIfNeeded(rec: SAMRecord): Unit = if (reversePerBaseTags && rec.getReadNegativeStrandFlag) {
    ConsensusTags.PerBase.AllPerBaseTags.foreach { tag =>
      val value = rec.getSignedShortArrayAttribute(tag)
      if (value != null) {
        reverse(value)
        rec.setAttribute(tag, value)
      }
    }
  }

  /** Reverses a short[]. */
  private def reverse(arr: Array[Short]): Unit = if (arr != null) {
    val lastIndex = arr.length - 1
    forloop (from=0, until=lastIndex/2) { i =>
      val tmp = arr(i)
      arr(i) = arr(lastIndex - i)
      arr(lastIndex - i) = tmp
    }
  }
}
