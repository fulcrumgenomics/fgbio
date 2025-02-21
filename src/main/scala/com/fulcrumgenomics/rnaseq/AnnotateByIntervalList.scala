/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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
package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{FilePath, PathToBam, PathToIntervals}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.rnaseq.OverlapAssociationStrategy.{Duplicate, LeastSpecific, MostSpecific}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.util.{CoordMath, Interval, IntervalList, OverlapDetector}

@clp(group = ClpGroups.Fastq, description =
  """
    |Collects base and read coverage statistics across intervals by name.
    |A regular expression may used to classify the name of the intervals to be counted and how they should be grouped.  Otherwise, the name is used as-is.
    |
    |For example for the following interval list:
    |
    |```
    |chr1    981971  982021  +       utr
    |chr1    982022  982065  +       intron
    |chr1    982066  982117  -       CDS;1
    |chr1    998963  999059  +       utr
    |chr1    999060  999432  -       CDS;2
    |chr1    999433  999526  +       intron
    |chr1    999527  999613  -       CDS;3
    |chr1    999614  999692  +       intron
    |chr1    999693  999973  -       CDS;4
    |```
    |
    |Using the the regex `([^;]*).*` would produce 3 separate interval counts CDS, utr and intron.
    |
    |Template bases (paired-end) or read bases (single end) are only counted if they fully enclosed within an interval.
    |By default if a read is enclosed by multiple intervals it will be counted for both intervals. This behavior
    |can be changed using a different `--overlap-association-strategy`.
  """)
class AnnotateByIntervalList
(
  @arg(flag = 'i', doc = "Path to a mapped BAM file.") val input: PathToBam,
  @arg(flag = 'l', doc = "Path to intervals") val intervals: PathToIntervals,
  @arg(flag = 'o', doc = "Path to the output metrics file") val output: FilePath,
  @arg(flag = 'R', doc = "The regular expression to use to classify the name of the interval") val nameRegex: Option[String] = None,
  @arg(flag = 'u', doc = "The name of the intervals not in the given set of intervals") val unmatched: String = "unmatched",
  @arg(flag = 'm', doc = "The minimum mapping quality for a read to be included") val minMappingQuality: Int = 20,
  @arg(flag = 'q', doc = "The minimum base quality for a base to be included") val minBaseQuality: Int = 0,
  @arg(flag = 's', doc = "The maximum insert size for an un-fragmented template") val maxInsertSize: Int = 1000,
  @arg(flag = 'a', doc = "The strategy to determine which interval overlapping templates should be associated with")
  val overlapAssociationStrategy: OverlapAssociationStrategy = OverlapAssociationStrategy.Duplicate
) extends FgBioTool with LazyLogging {

  Io.assertReadable(this.input)
  Io.assertReadable(this.intervals)
  Io.assertCanWriteFile(this.output)

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, verb = "Processed", noun = "loci", unit = 1000000)
    val reader = SamSource(this.input)
    val dict = reader.dict

    val (counters, overlapDetector) = {
      val baseCounter            = new SimpleCounter[String]()
      val uniqueBaseCounter      = new SimpleCounter[String]()
      val templateCounter        = new SimpleCounter[String]()
      val uniqueTemplateCounter  = new SimpleCounter[String]()
      val intervals              = buildIntervals(dict.asSam)
      intervals.foreach { interval => baseCounter.count(interval.getName, 0) }

      val _overlapDetector = OverlapDetector.create[Interval](intervals.iterator.toJavaList)
      val _counters: Map[CounterType, SimpleCounter[String]] = Map(
        CounterType.Base            -> baseCounter,
        CounterType.UniqueBase      -> uniqueBaseCounter,
        CounterType.Template        -> templateCounter,
        CounterType.UniqueTemplate  -> uniqueTemplateCounter
      )

      (_counters, _overlapDetector)
    }

    // Templates are assigned to an interval together
    val iterator = reader.iterator

    val records = iterator.filterNot { rec =>
      rec.duplicate || rec.secondary || rec.supplementary || !rec.pf || (rec.paired && !rec.isFrPair) || rec.mapq < minMappingQuality || pairStartTieBreaker(rec) || getInsertSize(rec) > maxInsertSize
    }

    val enclosingIntervalToRecords: Map[Interval, SamRecord] = records.map(record =>
      // Interval should be end of read one to the beginning of read two
      new Interval(
        record.refName,
        // Single end read mateStart will be set to 0
        Math.min(record.start, record.mateStart),
        // Set mate end to math min if it isn't defined (single end read)
        Math.max(record.end, record.mateEnd.getOrElse(Int.MinValue))
      ) -> record
    ).toMap

    enclosingIntervalToRecords.foreach { intervalToRecord =>
      // Overlaps should only be flagged when they enclose the template bases
      val recordInterval = intervalToRecord._1
      val record = intervalToRecord._2

      val overlaps = overlapDetector.getOverlaps(recordInterval).toSeq.filter { interval =>
        CoordMath.encloses(
          interval.getStart,
          interval.getEnd,
          recordInterval.getStart,
          recordInterval.getEnd)
      }

      if (overlaps.isEmpty) {
        counters(CounterType.Base).count(this.unmatched, getCountOfSequencedBases(record))
        counters(CounterType.Template).count(this.unmatched, 1)
      } else {
        val names: Seq[String] = overlaps.map(_.getName).distinct

        if (names.length > 1) {
          // If there are multiple enclosing intervals by name use an associate strategy
          overlapAssociationStrategy match {
            case Duplicate      => countForNames(namesToCount = names, record = record, counters = counters)
            case MostSpecific   => countForNames(namesToCount = Seq(overlaps.minBy(_.length).getName), record = record, counters = counters)
            case LeastSpecific  => countForNames(namesToCount = Seq(overlaps.maxBy(_.length).getName), record = record, counters = counters)
          }
        } else {
          // When there is only one interval matched we associate all intersected bases with the single interval
          countForNames(namesToCount = names, record = record, counters = counters)
        }
      }
      progress.record(recordInterval.getContig, recordInterval.getStart)
    }
    progress.logLast()
    reader.close()

    val numLociMultiMatched = counters(CounterType.Template).total - counters(CounterType.UniqueTemplate).total
    logger.info(f"Found $numLociMultiMatched%,d templates with multiple classifications.")

    // Get the total # of bases covered by each set of intervals
    val lengthByIntervalName = new SimpleCounter[String]()
    overlapDetector.getAll.toSeq.groupBy(_.getName).foreach { case (name, intvs) =>
      val intervalList = new IntervalList(dict.asSam)
      intervalList.addall(intvs.iterator.toJavaList)

      val length = intervalList.uniqued().map(_.length()).sum
      lengthByIntervalName.count(name, length)
    }
    lengthByIntervalName.count(this.unmatched, dict.length - lengthByIntervalName.total)

    val totalBases  = counters(CounterType.Base).total.toDouble
    val metrics     = counters(CounterType.Base).map { case (name, bases) =>
      // The total # of genomic bases covered by the intervals with the given name
      val span = lengthByIntervalName.countOf(name)

      CoverageByIntervalMetrics(
        name              = name,
        bases             = bases,
        unique_bases      = counters(CounterType.UniqueBase).countOf(name),
        span              = span,
        frac_genome       = span / dict.length.toDouble,
        frac_total_bases  = if (totalBases > 0) bases / totalBases else 0f,
        mean_coverage     = if (span > 0) bases / span.toDouble else 0f,
        templates         = counters(CounterType.Template).countOf(name),
        unique_templates  = counters(CounterType.UniqueTemplate).countOf(name)
      )
    }.toSeq.sortBy(_.name)
    Metric.write(this.output, metrics)
  }

  private def getInsertSize(rec: SamRecord): Int = {
    if (rec.paired && rec.mateMapped && (rec.mateRefName == rec.refName)) {
      Math.max(
        Math.abs(rec.start - rec.mateEnd.getOrElse(rec.start - rec.end)),
        Math.abs(rec.mateStart - rec.end)
      )
    } else 0
  }

  /**
    * If a record is paired, both reads are mapped and on the same contig and their start positions are identical,
    * then return true if for the second of the pair. This method is meant to be used in a `filterNot` such that
    * the first of the pair will be selected.
    * @param rec The record to check both reads for identical starts on
    * @return Returns true if both reads are mapped on the same contig and have the same start location.
    */
  private def pairStartTieBreaker(rec: SamRecord): Boolean = {
    rec.paired && rec.mapped && rec.mateMapped && rec.refName == rec.mateRefName && rec.start == rec.mateStart && rec.secondOfPair
  }

  /**
    * Gets the count of sequenced bases for a record. This does not include soft clipped bases.
    * @param record The record to count bases for.
    * @return The count of sequenced bases for this record.
    */
  private def getCountOfSequencedBases(record: SamRecord): Int = {
    val mateSequencedBases = if (record.mateMapped) {
      // If the mate is mapped and there is no mateEnd it should fail with the require, but just to be safe set
      // length to 0
      record.mateEnd.getOrElse(record.mateStart) - record.mateStart + 1
    } else {
      0
    }
    (record.end - record.start + 1) + mateSequencedBases
  }

  private def countForNames(
    namesToCount: Seq[String],
    record: SamRecord,
    counters: Map[CounterType, SimpleCounter[String]],
  ): Unit = {
    namesToCount.foreach { name =>
      val bases = getCountOfSequencedBases(record)
      counters(CounterType.Base).count(name, bases)
      counters(CounterType.Template).count(name, 1)
      // If we only have one enclosing interval count unique bases and templates as well
      if (namesToCount.length == 1) {
        counters(CounterType.UniqueBase).count(name, bases)
        counters(CounterType.UniqueTemplate).count(name, 1)
      }
    }
  }

  private def buildIntervals(dict: SAMSequenceDictionary): List[Interval] = {
    val inputIntervalList = IntervalList.fromPath(this.intervals)
    this.nameRegex.map(_.r) match {
      case None        => inputIntervalList.uniqued(true).toList
      case Some(regex) =>
        val intervals = inputIntervalList.map { i =>
          val name = i.getName match {
            case regex(n) => n
            case _        => throw new IllegalArgumentException(s"Could not match name '${i.getName}' with regex '${this.nameRegex}'")
          }
          new Interval(i.getContig, i.getStart, i.getEnd, i.isNegativeStrand, name)
        }.toJavaList
        val intervalList = new IntervalList(dict)
        intervalList.addall(intervals)
        intervalList.toList
    }
  }
}

/** Metric for [[AnnotateByIntervalList]].
  *
  * @param name             the name of the intervals
  * @param bases            the number of bases mapping to intervals with the given name.
  * @param unique_bases     the number of bases mapping exclusively to intervals with the given name,
  * @param span             the number of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_genome      the fraction of unique genomic loci/bases coverage by intervals with the given name.
  * @param frac_total_bases the fraction of mapped bases covered by intervals with the given name.
  * @param mean_coverage    the mean coverage of intervals with the given name.
  * @param templates        the number of templates mapping to the intervals with the given name
  * @param unique_templates the number of templates mapping exclusively to the intervals with the given name
  */
case class CoverageByIntervalMetrics(
  name: String,
  unique_bases: Long,
  bases: Long,
  span: Long,
  frac_genome: Double,
  frac_total_bases: Double,
  mean_coverage: Double,
  templates: Long,
  unique_templates: Long
) extends Metric


sealed trait CounterType extends EnumEntry
object CounterType extends FgBioEnum[CounterType] {
  def values: IndexedSeq[CounterType] = findValues
  case object Base extends CounterType
  case object UniqueBase extends CounterType
  case object Template extends CounterType
  case object UniqueTemplate extends CounterType
}
/**
  * Enum for determining how to associate a read that is enclosed in more than a single interval.
  */
sealed trait OverlapAssociationStrategy extends EnumEntry

object OverlapAssociationStrategy extends FgBioEnum[OverlapAssociationStrategy] {
  def values: IndexedSeq[OverlapAssociationStrategy] = findValues

  /**
    * Duplicate counting associates the template base with all enclosing intervals.
    */
  case object Duplicate extends OverlapAssociationStrategy

  /**
    * MostSpecific counting associates the template bases with the smallest enclosing interval.
    */
  case object MostSpecific extends OverlapAssociationStrategy

  /**
    * LeastSpecific counting associates the template bases with the largest enclosing interval.
    */
  case object LeastSpecific extends OverlapAssociationStrategy
}
