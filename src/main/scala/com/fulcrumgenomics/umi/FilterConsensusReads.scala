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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.reference.ReferenceSequence
import htsjdk.samtools.util.SequenceUtil

import java.io.Closeable
import java.lang.Math.{max, min}

/** Filter values for filtering consensus reads */
private[umi] case class ConsensusReadFilter(minReads: Int, maxReadErrorRate: Double, maxBaseErrorRate: Double)


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
      |  3. Reads with mean base quality of the consensus read, prior to any masking, less than min-mean-base-quality
      |     are filtered out (if specified)
      |  4. Bases with quality scores below min-base-quality are masked to Ns
      |  5. Bases with fewer than min-reads contributing raw reads are masked to Ns
      |  6. Bases with a consensus error rate (defined as the fraction of contributing reads that voted for a different
      |     base than the consensus call) higher than max-base-error-rate are masked to Ns
      |  7. For duplex reads, if require-single-strand-agreement is provided, masks to Ns any bases where the base was
      |     observed in both single-strand consensus reads and the two reads did not agree
      |  8. Reads with a proportion of Ns higher than max-no-call-fraction *after* per-base filtering are filtered out
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
  val maxNoCallFraction: Double = 0.2,
  @arg(flag='q', doc="The minimum mean base quality across the consensus read.")
  val minMeanBaseQuality: Option[PhredScore] = None,
  @arg(flag='s', doc="Mask (make `N`) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only).")
  val requireSingleStrandAgreement: Boolean = false,
  @arg(flag='S', doc="The sort order of the output. If not given, output will be in the same order as input if the input is query grouped, otherwise queryname order.")
  val sortOrder: Option[SamOrder] = None
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
  private[umi] val ccFilters = ConsensusReadFilter( // For the final Consensus Calls
    minReads         = minReads.head,
    maxReadErrorRate = maxReadErrorRate.head,
    maxBaseErrorRate = maxBaseErrorRate.head
  )
  private[umi] val abFilters = ConsensusReadFilter(
    minReads         = minReads.take(2).last,
    maxReadErrorRate = maxReadErrorRate.take(2).last,
    maxBaseErrorRate = maxBaseErrorRate.take(2).last
  )
  private[umi] val baFilters = ConsensusReadFilter(
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
    logger.info("Reading the reference fasta into memory")
    val refMap = ReferenceSequenceIterator(ref, stripComments=true).map { ref => ref.getContigIndex -> ref}.toMap
    logger.info(f"Read ${refMap.size}%,d contigs.")

    val in  = SamSource(input)
    val out = buildOutputWriter(in.header, refMap)

    // Go through the reads by template and do the filtering
    val templateIterator = Bams.templateIterator(in, maxInMemory=MaxRecordsInMemoryWhenSorting)
    logger.info("Filtering reads.")
    templateIterator.foreach { template =>
      val r1 = template.r1.getOrElse(throw new IllegalStateException(s"${template.name} had no R1."))
      if (r1.paired) require(template.r2.isDefined, s"Paired read missing R2: ${template.name}")
      val primaryReadCount = if (template.r2.isDefined) 2 else 1
      totalReads += primaryReadCount

      // Reverse the tags on the reads if needs be
      if (reversePerBaseTags) {
        template.allReads.foreach { rec =>
          reverseConsensusTagsIfNeeded(rec)
          reverseComplementConsensusTagsIfNeeded(rec)
        }
      }

      // Filter and mask
      val r1Result = filterRecord(r1)
      val r2Result = if (r1Result.keepRead && template.r2.isDefined) filterRecord(template.r2.get) else EmptyFilterResult

      if (r1Result.keepRead && r2Result.keepRead) {
        keptReads   += primaryReadCount
        totalBases  += r1.length + template.r2.map(_.length).getOrElse(0)
        maskedBases += r1Result.maskedBases + r2Result.maskedBases
        out += r1
        template.r2.foreach { r => out += r }

        template.allSupplementaryAndSecondary.foreach { r =>
          val result = filterRecord(r)
          if (result.keepRead) {
            out += r
          }
        }
      }
    }

    logger.info("Finalizing the output")
    in.safelyClose()
    out.close()
    logger.info(f"Output $keptReads%,d of $totalReads%,d primary consensus reads.")
    logger.info(f"Masked $maskedBases%,d of $totalBases%,d bases in retained primary consensus reads.")
  }

  /** Builds the writer to which filtered records should be written.
    *
    *  If the input order is [[SamOrder.Queryname]] or query grouped, then the filtered records will also be in the same
    *  order.  So if the output order is specified AND does not match the the input order, sorting will occur.
    *
    *  If the input order is not [[SamOrder.Queryname]] or query grouped, then the input records will be resorted into
    *  [[SamOrder.Queryname]].  So if the output order is specified AND is not [[SamOrder.Queryname]], sorting will occur.
    *
    *  Otherwise, we can skip sorting!
    *
    * */
  private def buildOutputWriter(inHeader: SAMFileHeader, refMap: Map[Int, ReferenceSequence]): Writer[SamRecord] with Closeable = {
    val outHeader    = inHeader.clone()

    val inSortOrder  = inHeader.getSortOrder
    val inGroupOrder = inHeader.getGroupOrder
    val inSubSort    = Option(inHeader.getAttribute("SS"))

    // Get the order after filtering
    val (afterFilteringSortOrder, afterFilteringGroupOrder, afterFilteringSubSort) = {
      if (inSortOrder == SortOrder.queryname || inGroupOrder == GroupOrder.query) { // no sorting occurred, so same as input
        (inSortOrder, inGroupOrder, inSubSort)
      }
      else { // sorting occurred, so it's queryname
        val order = SamOrder.Queryname
        (order.sortOrder, order.groupOrder, order.subSort)
      }
    }

    // Get the desired output order
    val (outputSortOrder, outputGroupOrder, outputSubSort) = this.sortOrder match {
      case None        => (inSortOrder, inGroupOrder, inSubSort) // same as input
      case Some(order) => (order.sortOrder, order.groupOrder, order.subSort) // specific output
    }

    val sort: Option[SamOrder] = {
      // if the order after filtering and the output order match, no need to re-sort the output
      if (afterFilteringSortOrder == outputSortOrder && afterFilteringGroupOrder == outputGroupOrder && afterFilteringSubSort == outputSubSort) {
        None
      } else { // output order and order after filtering do not match, we need to re-sort the output
        SamOrder.values.find { order =>
          order.sortOrder == outputSortOrder && order.groupOrder == outputGroupOrder && order.subSort == outputSubSort
        }.orElse {
          // this can only happen if the input order is unrecognized
          throw new IllegalArgumentException(
            s"The input BAM had an unrecognized sort order (SO:$inSortOrder GO:$inGroupOrder SS: $inSubSort) in $input"
          )
        }
      }
    }
    val writer = SamWriter(output, outHeader, sort=sort, maxRecordsInRam=MaxRecordsInMemoryWhenSorting)
    sort.foreach(o => logger.info(f"Output will be sorted into $o order"))

    // Create the final writer based on if the full reference has been loaded, or not
    new Writer[SamRecord] with Closeable {
      override def write(rec: SamRecord): Unit = {
        Bams.regenerateNmUqMdTags(rec, refMap(rec.refIndex))
        writer += rec
      }
      def close(): Unit = writer.close()
    }
  }

  /**
    * Perform filtering on an individual read.
    *
    * @param  rec the read to filter
    * @return a tuple of (Boolean, Int) where the first value is true if the read should be kept
    *         and false if it should be filtered out, and the second value is the count of bases
    *         in the read that were masked to Ns as a result of per-base filtering
    */
  private[umi] def filterRecord(rec: SamRecord): FilterResult = {
    // Checks that are common to vanilla and duplex reads
    val maxDepth = rec.get[Int](ConsensusTags.PerRead.RawReadCount)
    val errRate  = rec.get[Float](ConsensusTags.PerRead.RawReadErrorRate)
    if (maxDepth.isEmpty || errRate.isEmpty) fail(s"Read ${rec.name} does not appear to have consensus calling tags present.")

    // Only bother looking at reads where per-read criteria are met
    if (maxDepth.exists(_ < this.ccFilters.minReads) ||
      errRate.exists(_ > this.ccFilters.maxReadErrorRate) ||
      this.minMeanBaseQuality.exists(_ > rec.quals.foldLeft(0)(_ + _).toDouble / rec.length)) {
      FilterResult(keepRead=false, 0)
    }
    else {
      val result = if (isDuplexRecord(rec)) filterDuplexConsensusRead(rec) else filterVanillaConsensusRead(rec)
      result.copy(keepRead = result.keepRead && fractionNoCalls(rec) <= this.maxNoCallFraction)
    }
  }

  /** Computes the fraction of the read that are no-calls. */
  private def fractionNoCalls(rec: SamRecord): Double = {
    val bases = rec.bases
    var ns = 0
    forloop(from = 0, until = bases.length) { i => if (bases(i) == NoCall) ns += 1 }
    ns / bases.length.toDouble
  }

  /** Performs filtering that is specific to Vanilla (single-umi/non-duplex) consensus reads. */
  private def filterVanillaConsensusRead(rec: SamRecord): FilterResult = {
    val bases  = rec.bases
    val quals  = rec.quals
    val depths = rec[Array[Short]](ConsensusTags.PerBase.RawReadCount)
    val errors = rec[Array[Short]](ConsensusTags.PerBase.RawReadErrors)
    val pb = depths != null && errors != null
    val filters = this.ccFilters

    // Do the per-base masking
    var maskedBasesThisRead = 0
    forloop(from = 0, until = rec.length) { i =>
      if ((quals(i) < this.minBaseQuality) ||
        (pb && depths(i) < filters.minReads) ||
        (pb && errors(i) / depths(i).toDouble > filters.maxBaseErrorRate)) {
        bases(i) = NoCall
        quals(i) = NoCallQual
        maskedBasesThisRead += 1
      }
    }

    rec.bases = bases
    rec.quals = quals
    FilterResult(keepRead=true, maskedBases=maskedBasesThisRead)
  }

  /** Performs filtering that is specific to Duplex consensus reads. */
  private def filterDuplexConsensusRead(rec: SamRecord): FilterResult = {
    val failsReadLevelChecks = {
      import ConsensusTags.PerRead._
      val Seq(baMaxDepth, abMaxDepth) = Seq(AbRawReadCount,     BaRawReadCount    ).map(rec.apply[Int]).sorted
      val Seq(abError, baError)       = Seq(AbRawReadErrorRate, BaRawReadErrorRate).map(rec.apply[Float]).sorted

      abMaxDepth < abFilters.minReads || abError > abFilters.maxReadErrorRate ||
        baMaxDepth < baFilters.minReads || baError > baFilters.maxReadErrorRate
    }

    if (failsReadLevelChecks) {
      FilterResult(keepRead=false, maskedBases=0)
    }
    else {
      val bases   = rec.bases
      val quals   = rec.quals

      // NB: per-base values may not exist for the BA-strand if the consensus was generated from a single-strand.
      val (abValues, baValues) = DuplexConsensusPerBaseValues.getPerBaseValues(rec)

      var maskedBases = 0
      forloop(from=0, until=bases.length) { i =>
        if (bases(i) != NoCall) {
          val maskBase = DuplexConsensusPerBaseValues.maskBaseAt(abValues, baValues, idx=i, ccFilters=this.ccFilters, abFilters=this.abFilters, baFilters=this.baFilters)
          val singleStrandDisagree = this.requireSingleStrandAgreement && !DuplexConsensusPerBaseValues.ssAgrees(idx=i, abValues=abValues, baValues=baValues)
          if (maskBase || singleStrandDisagree || quals(i) < this.minBaseQuality) {
            bases(i) = NoCall
            quals(i) = PhredScore.MinValue
            maskedBases += 1
          }
        }
      }

      rec.bases = bases
      rec.quals = quals
      FilterResult(keepRead=true, maskedBases=maskedBases)
    }
  }

  /** Quick check to see if a record is a duplex record. */
  private def isDuplexRecord(rec: SamRecord): Boolean = rec.get(ConsensusTags.PerRead.AbRawReadCount).isDefined

  /** Reverses all the consensus tags if the right conditions are set. */
  private def reverseConsensusTagsIfNeeded(rec: SamRecord): Unit = if (reversePerBaseTags && rec.negativeStrand) {
    ConsensusTags.PerBase.TagsToReverse.foreach { tag =>
      rec.get[Any](tag).foreach {
        case short: Array[Short] => rec(tag) = reverse(short)
        case string: String      => rec(tag) = string.reverse
        case _ => throw new IllegalArgumentException(s"Tag '$tag' must be of type Array[Short] or String, was ${tag.getClass.getSimpleName}")
      }
    }
  }

  /** Reverse complements all the consensus tags if the right conditions are set. */
  private def reverseComplementConsensusTagsIfNeeded(rec: SamRecord): Unit = if (reversePerBaseTags && rec.negativeStrand) {
    ConsensusTags.PerBase.TagsToReverseComplement.foreach { tag =>
      rec.get[String](tag).foreach { value =>
        rec(tag) = SequenceUtil.reverseComplement(value)
      }
    }
  }

  /** Reverses an short[]. */
  private def reverse(arr: Array[Short]): Array[Short] = if (arr != null) {
    val lastIndex = arr.length - 1
    forloop (from=0, until=lastIndex/2) { i =>
      val tmp = arr(i)
      arr(i) = arr(lastIndex - i)
      arr(lastIndex - i) = tmp
    }
    arr
  } else { arr }
}

object DuplexConsensusPerBaseValues {
  import ConsensusTags.{PerBase, PerRead}

  /** Gets the per-base values for the ab and ba strand respectively.  If the duplex consensus was created from a single
    * strand, the record should have values for only the ab-strand.  In this case, the values for the ba strand are set
    * as follows: the bases are set to "N", the depths and errors are set to 0. */
  def getPerBaseValues(rec: SamRecord): (DuplexConsensusPerBaseValues, DuplexConsensusPerBaseValues) = {
    /** Gets the value of the given tag, ensuring that it exists. */
    def getAndRequire[A](tag: String): A = rec.get[A](tag).getOrElse {
      throw new IllegalStateException(s"Expected tag '$tag' to exists for record: $rec")
    }

    // The AB-strand should always have values
    val abValues = {
      require(rec[Int](PerRead.AbRawReadCount) > 0, s"Expected the AB-strand to have depth > 0 for read: $rec")
      DuplexConsensusPerBaseValues(
        abStrand = true,
        bases    = getAndRequire[String](PerBase.AbConsensusBases),
        depths   = getAndRequire[Array[Short]](PerBase.AbRawReadCount),
        errors   = getAndRequire[Array[Short]](PerBase.AbRawReadErrors)
      )
    }

    // The BA-strand should have values if the depth was > 0 for some position in the read (i.e. a single-strand
    // consensus was generated for the BA-strand).
    val baValues = if (rec[Int](PerRead.BaRawReadCount) > 0) {
      DuplexConsensusPerBaseValues(
        abStrand = false,
        bases    = getAndRequire[String](PerBase.BaConsensusBases),
        depths   = getAndRequire[Array[Short]](PerBase.BaRawReadCount),
        errors   = getAndRequire[Array[Short]](PerBase.BaRawReadErrors)
      )
    }
    else {
      DuplexConsensusPerBaseValues(
        abStrand = !abValues.abStrand,
        bases    = "N" * abValues.bases.length,
        depths   = Array.range(start=0, end=abValues.bases.length).map(_ => 0.toShort),
        errors   = Array.range(start=0, end=abValues.bases.length).map(_ => 0.toShort)
      )
    }
    (abValues, baValues)
  }

  /** Returns false if both the ab and ba strand had depth > 0 at the given position and their base call disagreed, otherwise true */
  def ssAgrees(idx: Int, abValues: DuplexConsensusPerBaseValues, baValues: DuplexConsensusPerBaseValues): Boolean = {
    abValues.depths(idx) == 0 || baValues.depths(idx) == 0 || abValues.bases(idx) == baValues.bases(idx)
  }

  /** Returns true if the base at the given index should be masked. */
  def maskBaseAt(abValues: DuplexConsensusPerBaseValues, baValues: DuplexConsensusPerBaseValues, idx: Int,
                 ccFilters: ConsensusReadFilter, abFilters: ConsensusReadFilter, baFilters: ConsensusReadFilter): Boolean = {
    val abDepth    = max(abValues.depths(idx), baValues.depths(idx))
    val baDepth    = min(abValues.depths(idx), baValues.depths(idx))
    val abError    = min(abValues.errors(idx)/abValues.depths(idx).toDouble, baValues.errors(idx) / baValues.depths(idx).toDouble)
    val baError    = max(abValues.errors(idx)/abValues.depths(idx).toDouble, baValues.errors(idx) / baValues.depths(idx).toDouble)
    val totalDepth = abDepth + baDepth
    val totalError = (abValues.errors(idx) + baValues.errors(idx)) / totalDepth.toDouble

    totalDepth < ccFilters.minReads || totalError > ccFilters.maxBaseErrorRate ||
      abDepth  < abFilters.minReads || abError    > abFilters.maxBaseErrorRate ||
      baDepth  < baFilters.minReads || baError    > baFilters.maxBaseErrorRate
  }
}

/** A little class to store the per-base values for ab and ba-strand*/
case class DuplexConsensusPerBaseValues(abStrand: Boolean, bases: String, depths: Array[Short], errors: Array[Short])
