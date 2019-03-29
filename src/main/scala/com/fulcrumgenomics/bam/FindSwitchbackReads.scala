/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

import java.lang.Math.{abs, max, min}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter}
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util._
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools.reference.ReferenceSequence
import htsjdk.samtools.util.{CoordMath, SequenceUtil}
import htsjdk.samtools.{CigarOperator, SAMFileHeader, SAMReadGroupRecord}


object FindSwitchbackReads {
  /** Resource path for R script for plotting. */
  private val PlottingScript = "com/fulcrumgenomics/bam/FindSwitchbackReads.R"

  /** The name of the BAM extended attribute tag used to store switchback information. */
  val SwitchTag: String = "sb"

  /** Parent trait for the different kinds of hits we can find. */
  sealed trait SwitchHit { def code: Char }

  /**
    * Information about a switchback found by looking at soft-clipped complementary sequence within
    * a single read.
    *
    * @param offset the offset between the reference position at the last base read on the original strand and
    *               the first base read on the opposite strand.
    * @param length the length of the post-switch sequence within the read
    * @param read optionally the read in which the hit was found (generally filled in later)
    */
  case class ReadBasedHit(var offset: Int, length: Int, var read: Option[SamRecord] = None) extends SwitchHit {
    override val code: Char = 'r'
  }

  /**
    * Information about a switchback found by examining a template with a tandem (i.e. FF or RR) pair orientation.
    * The offset in this case is the distance between the end of one read and the beginning of the next.  E.g.
    *
    * 01234567890123456789
    * ----->  ----->       gap = 2
    *  <----- <-----       gap = 1
    * ----->
    *     ------>          gap = -2
    * @param gap the unsequenced gap between the the 3' end of the "earlier" read and the 5' end of the "later" read.
    */
  case class TandemBasedHit(gap: Int) extends SwitchHit {
    override val code: Char = 't'
  }

  /**
    * Summary metrics regarding switchback reads found.
    *
    * @param sample the name of the sample sequenced.
    * @param library the name of the library sequenced.
    * @param templates the total number of templates (i.e. inserts, unique read names) seen in the input.
    * @param aligned_templates the number of templates that had at least one aligned read.
    * @param switchback_templates the number of templates identified as having a switchback event.
    * @param fraction_switchbacks the fraction of all templates that appear to have switchbacks in them.
    * @param read_based_switchbacks the count of switchback_templates that were identified by looking for soft-clipped
    *                               reverse complementary sequence at the ends of reads.
    * @param mean_length the mean length of the reverse complementary sequence in the `read_based_switchbacks`.
    * @param mean_offset the mean offset of the `read_based_switchbacks`.
    * @param tandem_based_switchbacks the count of switchback_templates that were identified because they had paired
    *                                 reads mapped in FF or RR orientations within the designated gap size range.
    * @param mean_gap the mean gap size for `tandem_based_switchbacks`
    */
  case class SwitchMetric(sample: String,
                          library: String,
                          templates: Count,
                          aligned_templates: Count,
                          switchback_templates: Count,
                          fraction_switchbacks: Proportion,
                          read_based_switchbacks: Count,
                          mean_length: Double,
                          mean_offset: Double,
                          tandem_based_switchbacks: Count,
                          mean_gap: Double
                         ) extends Metric

  /**
    * Attempts to find a switchback event based on soft-clipped sequence within the primary reads of the template.
    *
    * @param refs a map of reference sequence name to [[ReferenceSequence]] object containing the sequence.
    * @param template the template to be examined.
    * @param minLength the minimum length of switchback to be detected.
    * @param maxOffset the maximum offset on the reference between the end of the pre-switch sequence and the start
    *                  of the post-switch sequence. A value of `0` will cause this method to return `None`.
    * @param maxErrorRate the maximum rate of mismatches between the switched-back sequence and the reference
    * @return `Some(ReadBasedHit)` if a hit is found else `None`.
    */
  private[bam] def findReadBasedSwitchback(refs: Map[String, ReferenceSequence],
                                           template: Template,
                                           minLength: Int,
                                           maxOffset: Int,
                                           maxErrorRate: Double): Option[ReadBasedHit] = {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementation notes:
    //
    // Read based switchbacks are detected by looking for soft-clipped sequence at the 5' end of the
    // read - i.e. the start of F reads and the end of R reads.  To trigger further checks we require
    // that soft-clipping exist and be at least as long as `minLength`.
    //
    // The second part of the check is that the soft-clipped sequence matches the reverse-complement
    // of the reference, and that the match _starts_ reading on the opposite strand within `maxOffset`
    // of the 5'-most mapped base of the non-clipped sequence.  How we do this differs by strand.
    //
    // A key value we wish to determine is the `offset` - the gap between where we leave off reading
    // one strand, and where we start reading on the other strand.  This value is positive if
    // post-switch the enzyme starts outside the sequence read on the original strand and moves back
    // towards the last base sequenced on the original strand, and negative if the enzyme starts
    // reading the second strand recessed within the already-read sequence of the first strand.
    // The `find` method returns the offset within the reference window given, which is then
    // translated into the offset from the 5' alignment position in a strand dependent manner.
    ////////////////////////////////////////////////////////////////////////////////////////////////

    if (maxOffset <= 0) None else template.primaryReads.filter(_.mapped).flatMap { rec =>
      if (rec.positiveStrand && rec.cigar.head.operator == CigarOperator.S && rec.cigar.head.length >= minLength) {
        val ref         = refs(rec.refName)
        val bases       = Sequences.revcomp(rec.basesString.take(rec.cigar.head.length))
        val windowStart = max(1, rec.start - maxOffset)
        val windowEnd   = min(ref.length(), CoordMath.getEnd(rec.start, bases.length + maxOffset))
        val refBases    = new String(ref.getBases, windowStart - 1, CoordMath.getLength(windowStart, windowEnd))
        find(bases, refBases, maxErrorRate, preferLaterMatches=false)
          .map { hit => hit.read = Some(rec); hit.offset = -(hit.offset + windowStart - rec.start); hit }
      }
      else if (rec.negativeStrand && rec.cigar.last.operator == CigarOperator.S && rec.cigar.last.length >= minLength) {
        val ref         = refs(rec.refName)
        val bases       = Sequences.revcomp(rec.basesString.takeRight(rec.cigar.last.length))
        val windowStart = max(1, rec.end - bases.length - maxOffset + 1)
        val windowEnd   = min(ref.length(), rec.end + maxOffset)
        val refBases    = new String(ref.getBases, windowStart - 1, windowEnd - windowStart + 1)
        find(bases, refBases, maxErrorRate, preferLaterMatches=true)
          .map { hit => hit.read = Some(rec); hit.offset = -(rec.end - (hit.offset + windowStart + bases.length - 1)); hit }
      }
      else {
        None
      }
    }.find(_ => true)
  }

  /**
    * Attempts to find a switchback based on the orientation of the primary reads (FF or RR) and the distance between
    * the end of one read (in sequencing order) and the start of the next read (in sequencing order), i.e. the gap size.
    *
    * @param template the template to be examined.
    * @param maxGap the maximum allowable gap size
    * @return either `Some(TandemBasedHit)` if a hit can be found, else `None`
    */
  private[bam] def findTandemSwitchback(template: Template, maxGap: Int) : Option[TandemBasedHit] =
    if (maxGap <=0 ) None else template.pairOrientation match {
      case Some(PairOrientation.TANDEM) =>
        (template.r1, template.r2) match {
          case (Some(r1), Some(r2)) =>
            val (earlier, later) = if (r1.start <= r2.start) (r1, r2) else (r2, r1)
            val gap = later.start - earlier.end - 1
            if (abs(gap) <= maxGap) Some(TandemBasedHit(later.start - earlier.end - 1)) else None
        }
      case _ =>
        None
  }

  /**
    * Attempts to find the first ungapped match of the query sequence within the target sequence. "First" is either
    * from the start of the target sequence, or if `preferLaterMatches` is true, from the end of the target sequence.
    *
    * @param query the query sequence to be found
    * @param target the target sequence in which the query sequence is to be located
    * @param maxErrorRate the maximum rate of mismatches allowed when matching the sequences
    * @param preferLaterMatches if true prefer matches closer to the end of targetSequence, else prefer matches closer
    *                           to the start
    * @return a [[ReadBasedHit]] if a match is found, otherwise None.
    */
  def find(query: String, target: String, maxErrorRate: Double, preferLaterMatches: Boolean): Option[ReadBasedHit] = {
    val maxMismatches = Math.floor(query.length * maxErrorRate)
    val range = {
      if (preferLaterMatches) Range.inclusive(target.length-query.length, 0, step = -1)
      else Range.inclusive(0, target.length-query.length-1, step=1)
    }

    range.iterator.flatMap { targetOffset =>
      var mismatches = 0
      forloop (from=0, until=query.length) { i =>
        if (!SequenceUtil.basesEqual(query(i).toByte, target(targetOffset+i).toByte)) mismatches += 1
      }

      if (mismatches <= maxMismatches) Some(ReadBasedHit(offset=targetOffset, length=query.length))
      else None
    }.find(_ => true)
  }

  /** Returns the read number, either 1 or 2. */
  private def readNum(rec: SamRecord): Int = if (!rec.paired || rec.firstOfPair) 1 else 2

  /** Returns the concatenation of the unique values per read group, or the default if there are no values. */
  private def summarizeRg(header: SAMFileHeader, f: SAMReadGroupRecord => String, default: String = "n/a"): String = {
    val ss = header.getReadGroups.map(f).filter(_ != null).map(_.trim).filter(_.nonEmpty).toSeq.distinct.sorted
    if (ss.nonEmpty) ss.mkString(",") else default
  }
}

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Finds reads where a template switch occurred during library construction.  Some library construction methods,
    |notably ultra-low-input shotgun methods, are prone to template switching events that create molecules
    |(templates, inserts) that instead of being a linear copy of a segment of genomic DNA, instead are chimeras formed
    |by starting on one strand of DNA and then, switching to the opposite strand.  Frequently when the switch occurs
    |there may be a small offset between where the first strand was departed and the opposite strand joined.
    |
    |## Algorithm
    |
    |Templates that contain strand switch events (switch-backs) are found by this tool in two different ways:
    |
    |1. By looking at reads that contain soft-clipped bases at their 5' end that, when reverse complemented, matches the
    |   genome proximal to the 5'-most mapped base of the read.  We call these matches
    |   "read based switchbacks".  Finding read based switchbacks is based on several parameters:
    |
    |   1. `max-offset` controls how far away to search for the reverse-complemented sequence.  The default value of
    |      `35` allows matches to be found when the soft-clipped sequence matches the genome _starting_ at most 35bp
    |      from the 5' mapped position of the read, and reading in the opposite direction.
    |   2. `min-length` controls the minimum number of soft-clipped bases that must exist to trigger the search.
    |      Given that the search looks at `2 * max-offset` locations (default=70) it is important that `min-length`
    |      is set such that `4^min-length >> 2 * `max-offset` in order to avoid false positives.
    |   3. `max-error-rate` allows for some mismatches to exist between the soft-clipped sequence and the genome when matching.
    |
    |2. By identifying templates with `FF` or `RR` (aka tandem) orientation where it is surmised that the template
    |   switch occurred in the un-sequenced region of the template between R1 and R2.  We call these `tandem based
    |   switchbacks`.  This is controlled by a single parameter, `max-gap`, which causes the tool to only identify a
    |   tandem read pair as a switch-back _if_ the gap between the end of the first read and the start of the second
    |   read is `+/- max-gap`.
    |
    |By default, when a switch-back template is identified, the primary reads are made unmapped (and the original
    |alignment stored in the OA tag) and all secondary and supplementary alignments are discarded.  This can be
    |disabled with the `--dont-unmap` or `-d` option.
    |
    |All reads from a switch-back template are also tagged with an `sb` tag that describes the nature of the
    |switchback.  If the template was identified base on soft-clipped sequence within a read the format is:
    |
    |```sb:Z:r,[read|mate],{offset},{length}```
    |
    |If the template is identified due to it's tandem pair orientation then the format is:
    |
    |```sb:Z:t,{gap}```
    |
    |## Inputs and Outputs
    |
    |The tool takes as input a SAM or BAM file, and by default consumes from `stdin`.  The primary output is also a
    |SAM or BAM file, and defaults to compressed BAM on `stdout`.  This allows the tool to be run immediately after
    |an aligner in a pipe, e.g. `bwa mem ref.fa r1.fq r2.fq | fgbio -Xmx8g --ref=ref.fa | ...`.
    |
    |If the input is neither `queryname` sorted nor `queryname` grouped (i.e. all reads with the same name grouped
    |together) it will be sorted into `queryname` order by the tool before processing.
    |
    |By default the output BAM is produced in the order the reads were processed (i.e. the input ordering _or_
    |queryname sorted if sorting was required).  This can be overridden with the `--sort-order` option.
    |
    |A number of text files are also produced if the `--metrics` option is specified.  E.g. when specifying
    |`--metrics=s1.switchback` the following files are produced:
    |
    |1. `s1.switchback.summary.txt`: A table of summary metrics describing the number of reads, switchbacks, etc.
    |2. `s1.switchback.lengths.txt`: A table of the distribution of observed switchback lengths in read-based switchbacks.
    |3. `s1.switchback.offsets.txt`: A table of the distribution of observed offsets in read-based switchbacks.
    |4. `s1.switchback.gaps.txt`: A table of the distribution of gap lengths in tampl
    |5. `s1.switchback.plots.pdf`: A PDF containing plots of the distributions from 2-4.
    |
    |Note: because this tool accesses the reference genome in a random manner it pre-loads the entire reference fasta
    |into memory.  As a result the tool is best run with `-Xmx8g` to give it sufficient memory.
  """)
class FindSwitchbackReads
( @arg(flag='i', doc="Input BAM file.") val input: PathToBam = Io.StdIn,
  @arg(flag='o', doc="Output BAM file.") val output: PathToBam = Io.StdOut,
  @arg(flag='m', doc="Metrics output file.") val metrics: Option[PathPrefix] = None,
  @arg(flag='s', doc="Output sort order.") val sortOrder: Option[SamOrder] = None,
  @arg(flag='r', doc="Reference genome fasta file.") val ref: PathToFasta,
  @arg(flag='O', doc="Maximum offset between end the two segments of the read on the reference. Set to 0 to disable read-based checks.") val maxOffset: Int = 35,
  @arg(flag='g', doc="Maximum gap between R1 and R2 of tandem reads to call a template a switchback. Set to 0 to disable tandem-based checks.") val maxGap: Int = 500,
  @arg(flag='l', doc="Minimum match length of the switched back segment.") val minLength: Int = 6,
  @arg(flag='e', doc="Maximum mismatch error rate of switchback match to genome.") val maxErrorRate: Double = 0.1,
  @arg(flag='d', doc="IF true, do NOT unmap reads from switchback templates.") val dontUnmap: Boolean = false
) extends FgBioTool with LazyLogging {
  import FindSwitchbackReads._

  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)
  metrics.foreach(m => Io.assertCanWriteFile(m))

  validate(maxOffset    >= 0,   "max-offset must be >= 0")
  validate(maxGap       >= 0,   "max-gap must be >= 0")
  validate(minLength    >= 0,   "min-length must be >= 0")
  validate(maxErrorRate >= 0.0, "max-error-rate must be >= 0")

  override def execute(): Unit = {
    logger.info("Loading reference sequences.")
    val refs = ReferenceSequenceIterator(ref, stripComments=true).map(r => r.getName -> r).toMap

    val progress = ProgressLogger(logger)
    var total   = 0L
    var aligned = 0L
    val switchbackGaps    = new NumericCounter[Int]()
    val switchbackLengths = new NumericCounter[Int]()
    val switchbackOffsets = new NumericCounter[Int]()

    val in = SamSource(input)
    val iterator = Bams.templateIterator(in)
    val out = {
      val outHeader = in.header.clone()
      sortOrder match {
        case Some(so) => so.applyTo(outHeader)
        case None     =>
          if (in.header.getSortOrder != SortOrder.queryname && in.header.getGroupOrder != GroupOrder.query) {
            SamOrder.Queryname.applyTo(outHeader)
          }
      }

      SamWriter(output, header=outHeader, sort=sortOrder)
    }

    logger.info("Examining reads for switchbacks.")
    iterator.map { template =>
      total += 1
      if (template.primaryReads.forall(_.unmapped)) template else {
        aligned += 1

        val result: Option[SwitchHit] = findReadBasedSwitchback(refs, template, minLength, maxOffset, maxErrorRate)
          .orElse(findTandemSwitchback(template, maxGap))

        result.foreach {
          case hit: TandemBasedHit =>
            switchbackGaps.count(hit.gap)
            val tagValue = s"${hit.code},${hit.gap}"
            template.allReads.foreach(r => r(SwitchTag) = tagValue)
          case hit: ReadBasedHit =>
            switchbackLengths.count(hit.length)
            switchbackOffsets.count(hit.offset)
            val foundInReadNum = readNum(hit.read.getOrElse(unreachable("read should be populated!")))
            template.allReads.foreach { r =>
              val num = readNum(r)
              r(SwitchTag) = s"${hit.code},${if (num == foundInReadNum) "read" else "mate"},${hit.offset},${hit.length}"
            }
        }

        if (result.isDefined && !this.dontUnmap) template.unmapped else template
      }
    }
    .foreach { template =>
      template.allReads.foreach { rec =>
        out += rec
        progress.record(rec)
      }
    }

    in.safelyClose()
    out.close()

    // Output the metrics
    this.metrics.foreach { prefix =>
      val m = SwitchMetric(
        sample                   = summarizeRg(header=in.header, f = _.getSample),
        library                  = summarizeRg(header=in.header, f = _.getLibrary),
        templates                = total,
        aligned_templates        = aligned,
        switchback_templates     = switchbackLengths.total + switchbackGaps.total,
        fraction_switchbacks     = (switchbackLengths.total + switchbackGaps.total) / aligned.toDouble,
        read_based_switchbacks   = switchbackLengths.total,
        mean_length              = switchbackLengths.mean(),
        mean_offset              = switchbackOffsets.mean(),
        tandem_based_switchbacks = switchbackGaps.total,
        mean_gap                 = switchbackGaps.mean()
      )

      Metric.write(PathUtil.pathTo(prefix + ".summary.txt"), m)

      Seq(
        ("lengths", switchbackLengths, "length"),
        ("offsets", switchbackOffsets, "offset"),
        ("gaps",    switchbackGaps,    "gap_length")
      ).foreach { case (suffix, counter, header) =>
        val writer = Io.toWriter(PathUtil.pathTo(s"$prefix.$suffix.txt"))
        writer.write(s"$header\tcount\n")
        counter.foreach { case (x, count) => writer.write(s"$x\t$count\n")}
        writer.close()
      }

      Rscript.execIfAvailable(PlottingScript, input.getFileName.toString, s"$prefix.lengths.txt", s"$prefix.offsets.txt", s"$prefix.gaps.txt", s"$prefix.plots.pdf")
    }
  }
}

