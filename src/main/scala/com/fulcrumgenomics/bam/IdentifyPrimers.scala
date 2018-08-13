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

import java.io.File
import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathToBam, forloop}
import com.fulcrumgenomics.alignment.{Alignable, Aligner, Alignment, AlignmentTask, InteractiveAlignmentScorer, Mode => AlignmentMode}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{javaIteratorAsScalaIterator, unreachable, _}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{DelimitedDataParser, LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.util.{Interval, Locatable, OverlapDetector, SequenceUtil}

import scala.collection.immutable
import scala.reflect.runtime.{universe => ru}

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Identifies primers that generate the reads in the input BAM file.
    |
    |Note: when comparing a primer sequence to a read sequence, it "matches" if the two sequences are within at most
    |`--max-mismatches` to the assigned primers with at most `--max-mismatch-delta` difference between between number
    |of mismatches in the best and second best primer matches.
    |
    |## Primers
    |
    |The input primers file should contain a header, with one row per primer, as follows:
    |
    |  * pair_id - the unique identifier for the primer pair
    |  * primer_id - the unique identifier for the primer in a primer pair
    |  * sequence - the DNA sequence, including degenerate bases, in the sequencing order (i.e. if forward is False, the
    |               sequence should be the reverse complement of the forward genomic strand).
    |  * chr - the refNameosome/refName
    |  * start - the start position (1-based)
    |  * end - the end position (1-based, inclusive)
    |  * forward - True if the forward primer, false otherwise
    |
    |The primer file must contain headers, e.g:
    |
    |```
    |pair_id	primer_id	sequence	chr	start	end	forward
    |1	1_1	GATTACA	chr1	1010873	1010894	True
    |1	1_2	ACATTAG	chr1	1011118	1011137	True
    |```
    |
    |Note: a canonical primer pair is the pair of primers with the same `pair_id`.  Other possible pairs are called
    |"non-canonical".
    |
    |## Paired-End Matching (5' only)
    |
    |If the input are paired-end reads, then the reads are matched as follows for each read pair:
    |
    |1. Assign to a canonical primer pair if the 5' end of the mates agree with the mapping
    |location of their respective primers.
    |2. Assign to a non-canonical primer pair if the 5' end of the mates agree with the mapping location of their
    |respective primers, and the primer and read sequences "match" respectively for each primer.
    |3. Assign each mate of a read pair to a primer based on fully-gapped alignment.
    |
    |## Paired-End Matching (5' and 3')
    |
    |The `--three-prime` option can be used to also search the 3' end of every read for a primer sequence as follows:
    |
    |1. If a primer pair is identified using the 5' end, that pair is searched on the 3' end.
    |2. Otherwise, assign a primer based on fully-gapped alignment.
    |
    |The primer matched on the 5' end of read one must match the primer matched on the 3' end of read two, and the
    |primer matched on the 3' end of read one must match the primer matched on the 5' end of read two.  This will only
    |be enforced if both ends of a read match a primer.
    |
    |## Fragment Matching (5' only)
    |
    |If the input are fragment reads, then only the 5' end of a read is used to match against the input list of primers:
    |
    |1. Assign to a primer if the 5' end agrees with the mapping location of the respective primer, and the primer
    |and read sequence "match".
    |2. Assign to a primer based on fully-gapped alignment.
    |
    |## Fragment Matching (5' and 3')
    |
    |The 3' end of the read can also be searched to identify primer pairs using the `--three-prime` option.  In this
    |case the same algorithm will be employed as paired-end matching for the 5' only case, where each end of the fragment
    |read is equivalent to the 5' end of the mates of the paired end rad.
    |
    |## Unmapped Data
    |
    |If no reference is given or mapping information is not available, matching using mapping locations is skipped.
    |
    |## Speeding Up Alignment
    |
    |Install [the ksw executable](https://github.com/nh13/ksw) manually or via conda (`conda install -c bioconda ksw`)
    |and supply the path to the `ksw` executable via `--ksw`.
  """)
class IdentifyPrimers
(@arg(flag='i', doc="Input BAM file.")  val input: PathToBam,
 @arg(flag='p', doc="File containing information about the primers.") val primers: FilePath,
 @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Path prefix for the metrics files.") val metrics: PathPrefix,
 @arg(flag='S', doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
 @arg(          doc="Search for primers on the 3' end of each read as well (5' only by default).") val threePrime: Boolean = false,
 @arg(          doc="Maximum mismatches for a primer to be considered a match.") val maxMismatches: Int = 1,
 @arg(          doc="The minimum alignment score for fully gapped alignment for a primer to be considered a match.") val minAlignmentScore: Int = 10,
 @arg(          doc="The match score to use for aligning to primer sequences.") val matchScore: Int = 1,
 @arg(          doc="The mismatch score to use for aligning to primer sequences.") val mismatchScore: Int = -4,
 @arg(          doc="The gap open score to use for aligning to primer sequences.") val gapOpen: Int = -6,
 @arg(          doc="The gap extension score to use for aligning to primer sequences.") val gapExtend: Int = -1,
 @arg(          doc="The minimum difference in mismatches between the best match and second best match to accept a primer match.  Not used for full alignment.") val minMismatchDelta: Int = 1,
 @arg(flag='t', doc="The number of threads to use.") val threads: Int = 1,
 @arg(flag='n', doc="Examine the first N templates.") val numTemplates: Option[Int] = None,
 @arg(          doc="Skip full-alignment matching") val skipFullAlignment: Boolean = false,
 @arg(          doc="Path to the ksw aligner.") val ksw: Option[String] = None,
 @arg(          doc="The maximum number of templates in memory.") val maxTemplatesInRam: Option[Int] = None,
 @arg(          doc="The number of templates to process at a time per thread.") val templatesPerThread: Int = 1000
) extends FgBioTool with LazyLogging {

  import IdentifyPrimers._

  private val kswExecutable: Option[FilePath] = this.ksw.map(Paths.get(_)).flatMap {
    case p if Files.exists(p) => Some(p)
    case name =>
      val path = System.getenv("PATH")
      validate(path != null, "PATH environment variable was not found; required for --ksw")
      val systemPath = path.split(File.pathSeparatorChar).view.map(PathUtil.pathTo(_))
      systemPath.map(p => p.resolve(name)).find(ex => Files.exists(ex))
        .orElse {
          throw new ValidationException(s"Could not find ksw executable ${ksw} in PATH: $path")
        }
  }

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  Io.assertCanWriteFile(metrics)
  kswExecutable.foreach(Io.assertReadable)

  validate(maxMismatches >= 0,   "--max-mismatches must be >= 0")
  validate(minMismatchDelta > 0, "--min-mismatch-delta must be > 0")
  validate(matchScore > 0,       "--match-score must be > 0")
  validate(mismatchScore < 0,    "--mismatch-score must be < 0")
  validate(gapOpen < 0,          "--gap-open must be < 0")
  validate(gapExtend < 0,        "--gap-extend must be < 0")
  maxTemplatesInRam.foreach { m => validate(templatesPerThread < m, "--max-templates-in-ram must be greater than or equal to --templates-per-thread")}

  private val aligner: Aligner = Aligner(matchScore = matchScore, mismatchScore = mismatchScore, gapOpen = gapOpen, gapExtend = gapExtend, mode = AlignmentMode.Glocal)

  private val (fwdMatcher, revMatcher) = PrimerMatcher.build(
    path              = this.primers,
    aligner           = aligner,
    slop              = slop,
    maxMismatches     = maxMismatches,
    minMismatchDelta  = minMismatchDelta,
    minAlignmentScore = minAlignmentScore,
    withFullAlignment = !skipFullAlignment,
    revComp           = true
  )

  var numAlignments: AtomicLong = new AtomicLong(0)

  private val PrimerPairMatchTypeTag: String = "pp"
  private val PrimerInfoForwardTag: String   = "pf"
  private val PrimerInfoReverseTag: String   = "pr"
  private val NoPrimerMatchInfo: String      = "none"
  private val comments: Seq[String]          = {
    Seq(
      "The pf/pr tags store the primer match metadata for the forward and reverse strand respectively.",
      "The pf/pr tags are formatted as follow: <pair_id>,<primer_id>,<refName>:<start>-<end>,<strand>,<5-prime-offset>,<match-type>,<match-type-info>.",
      s"The match-type is 'location', 'mismatch', or 'aligmnent' based on if the match was found using the location, mismatch-based alignment, or full-alignment.",
      s"The pp tag is either 'canonical', 'non-canonical', 'single', or '$NoPrimerMatchInfo' based on if the matches for the primer pair agree, disagree, or are missing."
    )
  }

  override def execute(): Unit = {
    this.kswExecutable match {
      case Some(p) => logger.info(s"Using ksw aligner: $p")
      case None    => logger.info("Using the scala aligner.")
    }
    val metricCounter = new SimpleCounter[TemplateTypeMetric]()

    // Group the reads by template
    val in: SamSource = SamSource(this.input)
    val iterator: Iterator[Template] = {
      val iter = Bams.templateIterator(in)
      numTemplates.map { n => iter.take(n) }.getOrElse { iter }
    }

    // NB: Add comments explaining tags in the output writer
    val out: SamWriter = {
      val header = {
        val h = in.header.clone()
        comments.foreach { comment => h.addComment(comment) }
        h
      }
      SamWriter(path = output, header = header, sort = Some(SamOrder.Coordinate))
    }

    // Get the number of templates to process in a single thread
    val batchSize = numTemplates.map(_ / threads) match {
      case None    => templatesPerThread
      case Some(n) => math.max(1, math.min(n, templatesPerThread))
    }
    // The number of templates to batch across threads
    val majorBatchSize  = math.min(batchSize * threads * 4, maxTemplatesInRam.getOrElse(Int.MaxValue))
    val readingProgress = ProgressLogger(this.logger, verb="read", unit=batchSize*threads*2)

    // Batch templates to process them in individual threads.
    val outputIterator: Iterator[Seq[Template]] = if (threads > 1) {
      logger.info(f"Batching $majorBatchSize%,d templates with $batchSize%,d templates per thread.")

      import com.fulcrumgenomics.commons.CommonsDef.ParSupport
      // Developer Note: Iterator does not support parallel operations, so we need to group together records into a
      // [[List]] or [[Seq]].  A fixed number of records are grouped to reduce memory overhead.
      iterator
        .grouped(majorBatchSize)
        .flatMap { templates =>
          templates
            .grouped(batchSize)
            .toStream
            .parWith(threads, fifo = false)
            .map { templates => processBatch(templates, metricCounter, readingProgress) }
            .seq // ensures that the only parallelism is processBatch
            .toIterator
        }
    }
    else {
      logger.info(f"Batching $batchSize%,d templates.")
      iterator
        .grouped(batchSize)
        .map { templates => processBatch(templates, metricCounter, readingProgress) }
    }

    // Write the results
    val writingProgress = ProgressLogger(this.logger, "written", unit=1000000)
    outputIterator
      .flatMap(_.flatMap(_.allReads))
      .foreach { rec =>
        writingProgress.record(rec)
        out += rec
      }

    val rate = numAlignments.get() / readingProgress.getElapsedSeconds.toDouble
    logger.info(f"Performed ${numAlignments.get()}%,df full alignments in total ($rate%,.2f alignments/second).")
    logger.info(f"Wrote ${readingProgress.getCount}%,d records.")

    in.safelyClose()
    out.close()


    // Write metrics
    // Detailed metrics
    {
      val total   = metricCounter.total.toDouble
      val metrics = metricCounter.map { case (metric, count) =>
        metric.copy(count = count, frac = count / total)
      }.toSeq.sortBy(-_.count)
      Metric.write(PathUtil.pathTo(this.metrics + ".detailed.txt"), metrics)
    }
    // Summary metrics
    Metric.write(PathUtil.pathTo(this.metrics + ".summary.txt"), IdentifyPrimersMetric(metricCounter))
  }

  /** Creates a new [[InteractiveAlignmentScorer]]. */
  private def newAligner: InteractiveAlignmentScorer[Primer, SamRecordAlignable] = {
    InteractiveAlignmentScorer(matchScore, mismatchScore, gapOpen, gapExtend, AlignmentMode.Glocal, this.kswExecutable)
  }

  // NB: batching alignment inputs to the aligner (i.e. ksw) is empirically faster than given them all or just one at a time.
  private def completePrimerMatch(matchOption: Option[Either[PrimerMatch, Seq[AlignmentTask[Primer, SamRecordAlignable]]]],
                                  aligner: InteractiveAlignmentScorer[Primer, SamRecordAlignable],
                                  alignmentInputBatchSize: Int = 100): Option[PrimerMatch] = {
    case class PrimerAndScore(primer: Primer, score: Int)
    import IdentifyPrimers.MaxTwoBy
    matchOption.flatMap {
      case Left(primerMatch)      => Some(primerMatch)
      case Right(alignmentInputs) =>
        // Add the alignment inputs to the aligner
        alignmentInputs.grouped(alignmentInputBatchSize).foreach { batch =>
          aligner.append(batch: _*)
        }
        // Get the alignment results
        val alignmentResults = alignmentInputs.toIterator.zip(aligner.iterator)
          // Keep results that meet the minimum score and align with a minimum # of query bases
          .filter { case (alignmentInput, alignmentResult) =>
            val minQueryEnd = alignmentInput.queryLength - slop
            alignmentResult.score >= minAlignmentScore && alignmentResult.queryEnd >= minQueryEnd
          }
          // now just keep the primer and score
          .map { case (alignmentInput, alignmentResult) => PrimerAndScore(alignmentInput.query, alignmentResult.score) }
          // get the best two alignment results by score
          .maxTwoBy(_.score)

        // Create a primer match if alignments were found (i.e. a Option[PrimerMatch])
        alignmentResults match {
          case (None, _)                    => None
          case (Some(best), None)           => Some(FullAlignmentPrimerMatch(primer = best.primer, score = best.score, secondBestScore = minAlignmentScore))
          case (Some(best), Some(nextBest)) => Some(FullAlignmentPrimerMatch(primer = best.primer, score = best.score, secondBestScore = nextBest.score))
        }
    }
  }

  /** Processes a single batch of templates. */
  def processBatch(templates: Seq[Template],
                   metricCounter: SimpleCounter[TemplateTypeMetric],
                   progress: ProgressLogger): Seq[Template] = {

    val counter = new SimpleCounter[TemplateTypeMetric]()
    val aligner  = newAligner

    templates.foreach { template =>
      // match up to but not including full alignment
      val r1MatchOption = template.r1.map { r1 => fwdMatcher.matchUntilAlignment(r1) }
      val r2MatchOption = template.r2.map { r1 => revMatcher.matchUntilAlignment(r1) }

      // full in the alignment
      val r1Match = completePrimerMatch(r1MatchOption, aligner)
      val r2Match = completePrimerMatch(r2MatchOption, aligner)

      // add the three prime matching
      val r1ThreePrimeMatch = template.r1.flatMap { r1 => matchThreePrime(r1, r1Match, fwdMatch=true) }
      val r2ThreePrimeMatch = template.r2.flatMap { r2 => matchThreePrime(r2, r2Match, fwdMatch=false) }

      // Get information about the matches
      val templateTypes = (template.r1, template.r2) match {
        case (Some(r1), Some(r2)) =>
          val rType = (r1.mapped, r2.mapped) match {
            case (true, true)   => TemplateType.MappedPair
            case (true, false)  => TemplateType.Unpaired
            case (false, true)  => TemplateType.Unpaired
            case (false, false) => TemplateType.UnmappedPair
          }
          val mType = formatPair(r1, r2, r1Match, r2Match, r1ThreePrimeMatch, r2ThreePrimeMatch)
          TemplateTypeMetric(rType, mType, r1Match, r2Match)
        case (Some(r1), None) =>
          require(!r1.paired, s"Found paired read but missing R2 for ${r1.name}")
          val rType = if (r1.mapped) TemplateType.MappedFragment else TemplateType.UnmappedFragment
          val mType = formatFragment(r1, r1Match, r1ThreePrimeMatch)
          TemplateTypeMetric(rType, mType, r1Match, r2Match)
        case _ =>
          unreachable(s"Template did not have an R1: ${template.name}")
      }

      counter.count(templateTypes)
    }

    require(aligner.numAdded == aligner.numAligned, s"added: ${aligner.numAdded} alignments: ${aligner.numAligned}")

    aligner.close()

    numAlignments.addAndGet(aligner.numAligned)
    metricCounter.synchronized {
      addToSimpleCounter(src = counter, dest = metricCounter)
      templates.flatMap(_.allReads).foreach { rec =>
        if (progress.record(rec)) {
          val rate = numAlignments.get() / progress.getElapsedSeconds.toDouble
          logger.info(f"Performed ${numAlignments.get()}%,d full alignments so far ($rate%,.2f alignments/second).")
        }
      }
    }

    templates
  }

  // FIXME: add to SimpleCounter in commons
  private def addToSimpleCounter[T](src: SimpleCounter[T], dest: SimpleCounter[T]): Unit = {
    src.foreach { case (item, count) => dest.count(item, count) }
  }

  private def matchThreePrime(rec: SamRecord, fivePrimeMatch: Option[PrimerMatch], fwdMatch: Boolean): Option[PrimerMatch] = if (!threePrime) None else {
    // FIXME
    throw new NotImplementedError("Three-prime matching currently not supported, but is planned.")
  }

  /** Adds tags to the records based on the primer matching results. */
  private def formatPair(r1: SamRecord,
                         r2: SamRecord,
                         r1FivePrimeMatch: Option[PrimerMatch],
                         r2FivePrimeMatch: Option[PrimerMatch],
                         r1ThreePrimeMatch: Option[PrimerMatch],
                         r2ThreePrimeMatch: Option[PrimerMatch]): PrimerPairMatchType = {
    import PrimerPairMatchType._

    // Set the primer pair match type
    val matchType: PrimerPairMatchType = (r1FivePrimeMatch, r2FivePrimeMatch) match {
      case (Some(fwd), Some(rev))            => if (fwd.primer.pair_id == rev.primer.pair_id) Canonical else NonCanonical
      case (Some(_), None) | (None, Some(_)) => Single
      case _                                 => NoMatch
    }
    r1(PrimerPairMatchTypeTag) = matchType.toString

    // Set the per-primer/per-read match metadata
    val forwardInfo = r1FivePrimeMatch.map(_.info(r1, forward = true)).getOrElse(NoPrimerMatchInfo)
    val reverseInfo = r2FivePrimeMatch.map(_.info(r2, forward = false)).getOrElse(NoPrimerMatchInfo)
    Seq(r1, r2).foreach { rec =>
      rec(PrimerInfoForwardTag) = forwardInfo
      rec(PrimerInfoReverseTag) = reverseInfo
    }

    // TODO: three prime

    matchType
  }

  /** Adds tags to the record based on the primer matching results. */
  private def formatFragment(frag: SamRecord, fragFivePrimeMatch: Option[PrimerMatch], fragThreePrimeMatch: Option[PrimerMatch]): PrimerPairMatchType = {
    import PrimerPairMatchType._
    val forwardInfo = fragFivePrimeMatch.map(_.info(frag, forward = true)).getOrElse(NoPrimerMatchInfo)

    val matchType: PrimerPairMatchType = fragFivePrimeMatch match {
      case Some(pm) => Single
      case None     => NoMatch
    }

    frag(PrimerPairMatchTypeTag) = matchType.toString
    frag(PrimerInfoForwardTag)   = forwardInfo
    frag(PrimerInfoReverseTag)   = NoPrimerMatchInfo

    // TODO: three prime

    matchType
  }
}



object IdentifyPrimers {
  /** An implicit class that adds support for getting the maximum two values, if they exist. */
  implicit class MaxTwoBy[A](values: Iterator[A]) {
    def maxTwoBy[B](f: A => B)(implicit cmp: Ordering[B]): (Option[A], Option[A]) = {
      var best: Option[A]          = None
      var bestValue: Option[B]     = None
      var nextBest: Option[A]      = None
      var nextBestValue: Option[B] = None

      values.foreach { cur =>
        val curValue = f(cur)

        if (bestValue.forall ( v => cmp.lt(v, curValue))) {
          best = Some(cur)
          bestValue = Some(curValue)
        }
        else if (nextBestValue.forall ( v => cmp.lt(v, curValue))) {
          nextBest = Some(cur)
          nextBestValue = Some(curValue)
        }
      }

      (best, nextBest)

      // NB: this is slightly slower than the implementation above.
      /*
      values.foldLeft((Option.empty[A], Option.empty[A])) { case ((maybeBest, maybeNextBest), value) =>
        (maybeBest, maybeNextBest) match {
          case (None, _)                    => (Some(value), None)
          case (Some(best), None)           => if (cmp.lt(f(best), f(value))) (maybeBest, None) else (maybeBest, None)
          case (Some(best), Some(nextBest)) =>
            if (cmp.lt(f(best), f(value))) (Some(value), maybeBest)
            else if (cmp.lt(f(nextBest), f(value))) (maybeBest, Some(value))
            else (maybeBest, maybeNextBest)
        }
      }
      */
    }
  }



  object IdentifyPrimersMetric {

    def apply(templateTypesCounter: SimpleCounter[TemplateTypeMetric]): IdentifyPrimersMetric = {
      import PrimerPairMatchType._
      import TemplateType._

      val readTypeCounter    = new SimpleCounter[TemplateType]()
      val matchTypeCounter   = new SimpleCounter[PrimerPairMatchType]()
      val primerMatchCounter = new SimpleCounter[Class[_ <: PrimerMatch]]()

      templateTypesCounter.foreach { case (templateTypes, count) =>
        readTypeCounter.count(templateTypes.template_type, count)
        matchTypeCounter.count(templateTypes.primer_pair_match_type, count)
        Seq(templateTypes.r1_primer_match_type, templateTypes.r2_primer_match_type).flatten.foreach { primerMatchType =>
          primerMatchCounter.count(primerMatchType, count)
        }
      }

      val mapped_pairs       = readTypeCounter.countOf(MappedPair)
      val unpaired           = readTypeCounter.countOf(Unpaired)
      val unmapped_pairs     = readTypeCounter.countOf(UnmappedPair)
      val mapped_fragments   = readTypeCounter.countOf(MappedFragment)
      val unmapped_fragments = readTypeCounter.countOf(UnmappedFragment)

      new IdentifyPrimersMetric(
        templates          = readTypeCounter.total,
        pairs              = mapped_pairs + unpaired + unmapped_pairs,
        fragments          = mapped_fragments + unmapped_fragments,
        mapped_pairs       = mapped_pairs,
        unpaired           = unpaired,
        unmapped_pairs     = unmapped_pairs,
        mapped_fragments   = mapped_fragments,
        unmapped_fragments = unmapped_fragments,
        total_types        = matchTypeCounter.total,
        canonical          = matchTypeCounter.countOf(Canonical),
        non_canonical      = matchTypeCounter.countOf(NonCanonical),
        single             = matchTypeCounter.countOf(Single),
        no_match           = matchTypeCounter.countOf(NoMatch),
        total_matches      = primerMatchCounter.total,
        location           = primerMatchCounter.countOf(classOf[LocationBasedPrimerMatch]),
        mismatch           = primerMatchCounter.countOf(classOf[MismatchAlignmentPrimerMatch]),
        full_alignment     = primerMatchCounter.countOf(classOf[FullAlignmentPrimerMatch])
      )
    }
  }

  // TODO: document
  /**
    *
    * @param templates
    * @param pairs
    * @param fragments
    * @param mapped_pairs
    * @param unpaired
    * @param unmapped_pairs
    * @param mapped_fragments
    * @param unmapped_fragments
    * @param total_types
    * @param canonical
    * @param non_canonical
    * @param single
    * @param no_match
    * @param total_matches
    * @param location
    * @param mismatch
    * @param full_alignment
    */
  case class IdentifyPrimersMetric
  ( templates: Long = 0,
    pairs: Long = 0,
    fragments: Long = 0,
    mapped_pairs: Long = 0,
    unpaired: Long = 0,
    unmapped_pairs: Long = 0,
    mapped_fragments: Long = 0,
    unmapped_fragments: Long = 0,
    total_types: Long = 0,
    canonical: Long = 0,
    non_canonical: Long = 0,
    single: Long = 0,
    no_match: Long = 0,
    total_matches: Long = 0,
    location: Long = 0,
    mismatch: Long = 0,
    full_alignment: Long = 0
  ) extends Metric

  object TemplateTypeMetric {
    def apply(readType: TemplateType,
              primerPairMatchType: PrimerPairMatchType,
              r1PrimerMatchType: Option[PrimerMatch],
              r2PrimerMatchType: Option[PrimerMatch]
             ): TemplateTypeMetric = {
      TemplateTypeMetric(
        template_type              = readType,
        primer_pair_match_type = primerPairMatchType,
        r1_primer_match_type   = r1PrimerMatchType.map(_.getClass),
        r2_primer_match_type   = r2PrimerMatchType.map(_.getClass)
      )
    }
  }


  // TODO: document
  /**
    *
    * @param template_type
    * @param primer_pair_match_type
    * @param r1_primer_match_type
    * @param r2_primer_match_type
    * @param count
    * @param frac
    */
  case class TemplateTypeMetric
  (template_type: TemplateType,
   primer_pair_match_type: PrimerPairMatchType,
   r1_primer_match_type: Option[Class[_ <: PrimerMatch]],
   r2_primer_match_type: Option[Class[_ <: PrimerMatch]],
   count: Long = 0,
   frac: Double = 0d
  ) extends Metric

  sealed trait TemplateType extends EnumEntry

  // TODO: document
  object TemplateType extends FgBioEnum[TemplateType] {
    case object MappedPair extends TemplateType
    case object Unpaired extends TemplateType
    case object UnmappedPair extends TemplateType
    case object MappedFragment extends TemplateType
    case object UnmappedFragment extends TemplateType

    override def values: immutable.IndexedSeq[TemplateType] = findValues
  }

  sealed trait PrimerPairMatchType extends EnumEntry

  // TODO: document
  object PrimerPairMatchType extends FgBioEnum[PrimerPairMatchType] {
    case object Canonical extends PrimerPairMatchType
    case object NonCanonical extends PrimerPairMatchType
    case object Single extends PrimerPairMatchType
    case object NoMatch extends PrimerPairMatchType

    override def values: immutable.IndexedSeq[PrimerPairMatchType] = findValues
  }

  object PrimerMatch {
    val InfoDelimiter: String = ","
  }

  // TODO: document
  sealed trait PrimerMatch {
    def primer: Primer

    final def info(rec: SamRecord, forward: Boolean): String = {
      val baseInfo = Seq(
        primer.pair_id,
        primer.primer_id,
        primer.ref_name + ":" + primer.start + "-" + primer.end,
        if (forward) "+" else "-",
        0, // TODO: offset from the 5' end,
        this.toName
      ).map(_.toString)
      (baseInfo ++ this._info(rec, forward)).mkString(PrimerMatch.InfoDelimiter)
    }

    protected def _info(rec: SamRecord, forward: Boolean): Seq[Any]

    private def toName: String = this match {
      case _: LocationBasedPrimerMatch     => "location"
      case _: MismatchAlignmentPrimerMatch => "mismatch"
      case _: FullAlignmentPrimerMatch     => "alignment"
      case _                               => unreachable(s"Unknown primer match type: ${this.getClass.getSimpleName}.")
    }
  }
  // TODO: document
  case class LocationBasedPrimerMatch(primer: Primer, numMismatches: Int) extends PrimerMatch {
    protected def _info(rec: SamRecord, forward: Boolean): Seq[Any] = Seq(numMismatches)
  }

  // TODO: document
  case class MismatchAlignmentPrimerMatch(primer: Primer, numMismatches: Int, nextNumMismatches: Int) extends PrimerMatch {
    protected def _info(rec: SamRecord, forward: Boolean): Seq[Any] = {
      val nextOrNa = if (nextNumMismatches == Int.MaxValue) "na" else nextNumMismatches
      Seq(numMismatches, nextOrNa)
    }
  }

  // TODO: document
  case class FullAlignmentPrimerMatch(primer: Primer, score: Int, secondBestScore: Int) extends PrimerMatch {
    protected def _info(rec: SamRecord, forward: Boolean): Seq[Any] = Seq(score, secondBestScore)
  }

  object Primer {
    /** Writes the given primers to file.  Reverse strand primers are written with their sequence reverse complemented. */
    def write(path: FilePath, primers: TraversableOnce[Primer]): Unit = {
      val newPrimers = primers.map { primer =>
        if (primer.forward) primer
        else primer.copy(sequence = SequenceUtil.reverseComplement(primer.sequence))
      }
      Metric.write(path, newPrimers)
    }
  }

  /** A locatable class for a single primer of a primer pair.
    *
    * The bases should be 5'->3' on the genomic forward strand, which facilitates matching against [[SamRecord]]s.
    *
    * @param pair_id the canonical primer pair identifier, unique across all primer pairs.
    * @param primer_id the canonical primer identifier, unique across all primers.
    * @param sequence the primer sequence, in 5'->3' on the genomic forward strand.
    * @param ref_name the reference name to which this primer targets.
    * @param start the reference start position at which this primer starts.
    * @param end the reference end position at which this primer starts.
    * @param forward true if the primer maps to the forward genomic strand, false otherwise.
    */
  case class Primer(pair_id: String,
                    primer_id: String,
                    sequence: String,
                    ref_name: String,
                    start: Int,
                    end: Int,
                    forward: Boolean) extends Locatable with Alignable with Metric {
    override def getContig: String = ref_name
    override def getStart: Int = start
    override def getEnd: Int = end
    override def length: Int = end - start + 1
    def reverse: Boolean = !this.forward

    val bases: Array[Byte] = sequence.getBytes

    this.sequence.zipWithIndex.foreach { case (base, index) =>
      if (!SequenceUtil.isValidBase(base.toByte) && !SequenceUtil.isIUPAC(base.toByte)) {
        val prefix = sequence.substring(0, index)
        val base   = sequence.charAt(index)
        val suffix = sequence.substring(index+1)
        throw new IllegalArgumentException(s"Found an invalid base for primer $pair_id/$primer_id: $prefix[$base]$suffix")
      }
    }
  }

  /** Companion class to [[PrimerMatcher]] class that provides factory methods. */
  object PrimerMatcher {

    /** Builds a forward and reverse matcher from the primers at the given path. */
    def build(path: FilePath,
              aligner: Aligner,
              slop: Int,
              maxMismatches: Int,
              minMismatchDelta: Int,
              minAlignmentScore: Int,
              withFullAlignment: Boolean,
              revComp: Boolean = true): (PrimerMatcher, PrimerMatcher) = {
      val parser  = DelimitedDataParser(path, '\t')
      val primers = parser.map { row =>
        val forward  = row.apply[Boolean]("forward")
        val sequence = row.apply[String]("sequence")
        Primer(
          pair_id   = row.apply[String]("pair_id"),
          primer_id = row.apply[String]("primer_id"),
          sequence  = if (forward) sequence else SequenceUtil.reverseComplement(sequence),
          ref_name  = row.apply[String]("ref_name"),
          start     = row.apply[Int]("start"),
          end       = row.apply[Int]("end"),
          forward   = forward
        )
      }.toSeq
      build(primers, aligner, slop, maxMismatches, minMismatchDelta, minAlignmentScore, withFullAlignment)
    }

    /** Builds a forward and reverse matcher from the given primers. */
    def build(primers: Seq[Primer],
              aligner: Aligner,
              slop: Int,
              maxMismatches: Int,
              minMismatchDelta: Int,
              minAlignmentScore: Int,
              withFullAlignment: Boolean): (PrimerMatcher, PrimerMatcher) = {
      val (fwdPrimers: Seq[Primer], revPrimers: Seq[Primer]) = primers.partition(_.forward)
      val fwdMatcher = new PrimerMatcher(fwdPrimers, aligner, slop, maxMismatches, minMismatchDelta, minAlignmentScore, isForward=true, withFullAlignment=withFullAlignment)
      val revMatcher = new PrimerMatcher(revPrimers, aligner, slop, maxMismatches, minMismatchDelta, minAlignmentScore, isForward=false, withFullAlignment=withFullAlignment)

      // Validate a bunch of things
      validate(primers, fwdMatcher, revMatcher)

      (fwdMatcher, revMatcher)
    }

    /** Validates a forward and reverse matcher given a set of primers. */
    private def validate(primers: Seq[Primer], fowdMatcher: PrimerMatcher, revMatcher: PrimerMatcher): Unit = {
      val fwdDetector = fowdMatcher.detector
      val revDetector = revMatcher.detector

      // Validate that two primers exist for each pair, and that one is forward and the other is reverse
      primers.groupBy(_.pair_id).foreach { case (pairId, primersForPair) =>
        primersForPair match {
          case Seq(first, second) =>
            if (first.forward == second.forward) {
              val tpe = if (first.forward) "forward" else "reverse"
              throw new IllegalArgumentException(s"Found two $tpe primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
            }
          case _ =>
            throw new IllegalArgumentException(s"Found ${primersForPair.length} primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
        }
      }

      // Validate we do not have forward primers that have the same refName/start, or reverse primers with the same refName/end
      primers.foreach { primer =>
        val overlaps = (fowdMatcher.detector.getOverlaps(primer).iterator() ++ revDetector.getOverlaps(primer).iterator())
          .filter { overlap => primer != overlap && overlap.forward == primer.forward && overlap.ref_name == primer.ref_name }
          .filter { overlap => if (primer.forward) overlap.start == primer.start else overlap.end == primer.end }
        if (overlaps.nonEmpty) {
          throw new IllegalArgumentException(s"Found primers had the same location:\n$primer\n" + overlaps.map(_.toString).mkString("\n"))
        }
      }

      // Validate we do not have forward primers that have the same refName/start or reverse primers with the same refName/end
      primers.foreach { primer =>
        val overlaps = (fwdDetector.getOverlaps(primer).iterator() ++ revDetector.getOverlaps(primer).iterator())
          .filter { overlap => primer != overlap && overlap.forward == primer.forward && overlap.ref_name == primer.ref_name }
          .filter { overlap => if (primer.forward) overlap.start == primer.start else overlap.end == primer.end }
        if (overlaps.nonEmpty) {
          throw new IllegalArgumentException(s"Found primers had the same location:\n$primer\n" + overlaps.map(_.toString).mkString("\n"))
        }
      }
    }

    /** Counts the # of mismatches, allowing the primer to have IUPAC bases. */
    private[bam] def numMismatches(bases: Array[Byte], primer: Array[Byte]): Int = {
      // NB: the primer may be longer than bases!
      val length = math.min(bases.length, primer.length)
      var count = 0
      forloop(from = 0, until = length) { index =>
        if (!SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(bases(index), primer(index))) count += 1
      }
      count
      //bases.zip(primer).count { case (l, r) => !SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(l.toByte, r.toByte) }
    }
  }

  /** A little class that wraps a [[SamRecord]] to make it [[Alignable]] and also to defer extracting the bases to align.
    * NB: end is exclusive, and start and end are zero-based. */
  case class SamRecordAlignable(rec: SamRecord, start: Int, end: Int) extends Alignable {
    def bases: Array[Byte]   = rec.bases.slice(start, end)
    override def length: Int = end - start
  }

  /** A primer matcher.  Attempts to match based on location, then based on a mismatch-only alignment, then based on a full-alignment.
    *
    * @param primers the primers against which this matcher matches.
    * @param aligner the aligner for full-alignment-based matching.
    * @param slop the slop allowed by this matcher when matching based on location.
    * @param maxMismatches the maximum number of mismatches allowed for a location or mismatch-alignment based match.
    * @param minMismatchDelta the minimum difference between the best and next-best mismatch-alignment for a mismatch-alignment based match.
    * @param minAlignmentScore the minimum difference between the best and next-best mismatch-alignment for a mismatch-alignment based match.
    * @param isForward true if we are to align the left side of the read, otherwise the right side of the read
    * @param withFullAlignment true to perform full alignment, otherwise skip full alignment.
    */
  class PrimerMatcher(val primers: Seq[Primer],
                      val aligner: Aligner,
                      val slop: Int,
                      val maxMismatches: Int,
                      val minMismatchDelta: Int,
                      val minAlignmentScore: Int,
                      val isForward: Boolean,
                      val withFullAlignment: Boolean) {
    import PrimerMatcher._

    require(minMismatchDelta > 0)

    /** The number of primers. */
    def length: Int = primers.length

    /** Attempts to match the 5' end of the read with the given set of primers. This stops short of performing a full
      * alignment, and instead, produces an alignment task that can be performed by the caller.  This may be useful if
      * alignment is expensive or is required many times. */
    def matchUntilAlignment(rec: SamRecord): Either[PrimerMatch, Seq[AlignmentTask[Primer, SamRecordAlignable]]] = {
      matchOnLocation(rec).orElse { matchWithMismatchAlignment(rec) } match {
        case Some(primerMatch) => Left(primerMatch)
        case None              =>
          val inputs = if (isForward) {
            this.primers.iterator.map { primer =>
              val targetEnd = math.min(primer.sequence.length + slop, rec.length)
              AlignmentTask(query=primer, target=SamRecordAlignable(rec, 0, targetEnd))
            }
          }
          else {
            this.primers.iterator.map { primer =>
              val targetStart = math.max(0, rec.length - primer.length - slop)
              AlignmentTask(query=primer, target=SamRecordAlignable(rec, targetStart, rec.length))
            }
          }
          Right(inputs.toSeq)
      }
    }

    /** Attempts to match the 5' end of the read with the given set of primers. */
    def matchIt(rec: SamRecord): Option[PrimerMatch] = {
      matchOnLocation(rec)
        .orElse { matchWithMismatchAlignment(rec) }
        .orElse { if (withFullAlignment) matchWithFullAlignment(rec) else None }
    }

    /** The [[OverlapDetector]] for the primers. */
    private val detector: OverlapDetector[Primer] = {
      val d = new OverlapDetector[Primer](0, 0)
      this.primers.foreach { primer => d.addLhs(primer, primer) }
      d
    }

    /** Returns a location-based match, otherwise None */
    private[bam] def matchOnLocation(rec: SamRecord): Option[LocationBasedPrimerMatch] = if (rec.unmapped) None else {
      val pos = if (isForward) rec.unclippedStart else rec.unclippedEnd
      matchOnLocation(rec.refName, pos)
        .flatMap { p => toLocationBasedPrimerMatch(mismatchAlign(rec, p)) }
    }

    /** Returns a mismatch-based alignment match, otherwise None */
    private[bam] def mismatchAlign(rec: SamRecord, primer: Primer): Option[MismatchAlignmentPrimerMatch] = {
      // TODO: should slop be considered?
      val mm = {
        if (isForward) numMismatches(rec.bases, primer.bases)
        else numMismatches(rec.bases.drop(math.max(0, rec.length - primer.length)), primer.bases)
      }
      if (mm <= maxMismatches) Some(MismatchAlignmentPrimerMatch(primer, mm, Int.MaxValue)) else None
    }

    /** Converts a [[MismatchAlignmentPrimerMatch]] to a [[LocationBasedPrimerMatch]]. */
    private def toLocationBasedPrimerMatch(primerMatch: Option[MismatchAlignmentPrimerMatch]): Option[LocationBasedPrimerMatch] = {
      primerMatch.map { m => LocationBasedPrimerMatch(primer = m.primer, numMismatches = m.numMismatches) }
    }

    /** Determines if the primer is close enoguh to the given read defined by start/end. */
    private[bam] def locationCloseEnough(primer: Primer, pos: Int): Boolean = {
      if (isForward) {
        math.abs(primer.start - pos) <= slop
      }
      else {
        math.abs(primer.end - pos) <= slop
      }
    }

    /** Matches the read defined by refName/start/end based on location. */
    private[bam] def matchOnLocation(refName: String, pos: Int): Option[Primer] = {
      val interval = new Interval(refName, math.max(1, pos-slop), pos+slop) // TODO: off the end of the refNameosome?
      detector.getOverlaps(interval)
        .iterator()
        .filter { primer => locationCloseEnough(primer, pos)}
        .toSeq match {
        case Seq()       => None
        case Seq(primer) => Some(primer)
        case _primers    => unreachable(s"Found multiple primers for $refName:$pos\n${_primers.mkString("\n")}")
      }
    }

    /** Examines all primers to find the best mismatch-based alignment match. */
    private[bam] def matchWithMismatchAlignment(rec: SamRecord): Option[MismatchAlignmentPrimerMatch] = {
      val alignments = primers.flatMap { primer: Primer => mismatchAlign(rec, primer) }
      getBestMismatchAlignment(alignments)
    }

    /** Examines all primers to find the best full-alignment-based match. */
    private[bam] def getBestMismatchAlignment(alignments: Seq[MismatchAlignmentPrimerMatch]): Option[MismatchAlignmentPrimerMatch] = {
      if (alignments.isEmpty) None
      else {
        val best: MismatchAlignmentPrimerMatch = alignments.minBy(_.numMismatches)
        alignments.filter(_.numMismatches != best.numMismatches) match {
          case Seq() => Some(best.copy(nextNumMismatches = Int.MaxValue))
          case alns  =>
            val nextBest = alns.minBy(_.numMismatches)
            if (nextBest.numMismatches - best.numMismatches >= minMismatchDelta) {
              Some(best.copy(nextNumMismatches = nextBest.numMismatches))
            }
            else {
              None
            }
        }
      }
    }

    // A little clss to store an alignment to a given primer
    private case class AlignmentAndPrimer(alignment: Alignment, primer: Primer)

    /** A little helper method for matchWithFullAlignment to create a [[FullAlignmentPrimerMatch]]. */
    private def toFullAlignmentPrimerMatch(best: AlignmentAndPrimer, nextBestScore: Option[Int]): Option[FullAlignmentPrimerMatch] = {
      val primerMatch = FullAlignmentPrimerMatch(
        primer          = best.primer,
        score           = best.alignment.score,
        secondBestScore = nextBestScore.map(math.max(_, minAlignmentScore)).getOrElse(minAlignmentScore)
      )
      Some(primerMatch)
    }

    /** Returns a full-alignment-based match, otherwise None */
    private[bam] def matchWithFullAlignment(rec: SamRecord): Option[FullAlignmentPrimerMatch] = {
      // TODO: cache alignments?
      val alignments: IndexedSeq[AlignmentAndPrimer] = if (isForward) {
        this.primers.iterator.map { primer =>
          val target = rec.bases.slice(0, math.min(primer.sequence.length + slop, rec.length))
          AlignmentAndPrimer(this.aligner.align(query = primer.bases, target = target), primer)
        }.filter(_.alignment.score >= minAlignmentScore).toIndexedSeq
      }
      else {
        this.primers.iterator.map { primer =>
          val startIndex = math.max(0, rec.length - primer.bases.length - slop)
          val target = rec.bases.drop(startIndex)
          AlignmentAndPrimer(this.aligner.align(query = primer.bases, target = target), primer)
        }.filter(_.alignment.score >= minAlignmentScore).toIndexedSeq
      }

      if (alignments.isEmpty) None else {
        val best: AlignmentAndPrimer = alignments.maxBy(_.alignment.score)
        alignments.filter(_.alignment.score != best.alignment.score) match {
          case Seq() => toFullAlignmentPrimerMatch(best, None)
          case alns  =>
            val nextBest = alns.maxBy(_.alignment.score)
            if (best.alignment.score - nextBest.alignment.score >= minMismatchDelta) {
              toFullAlignmentPrimerMatch(best, Some(nextBest.alignment.score))
            }
            else {
              None
            }
        }
      }
    }
  }
}