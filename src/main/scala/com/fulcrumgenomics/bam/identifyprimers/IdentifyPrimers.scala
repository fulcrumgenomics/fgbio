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

package com.fulcrumgenomics.bam.identifyprimers

import java.io.File
import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam}
import com.fulcrumgenomics.alignment.{Alignable, Aligner, AlignmentTask, BatchAligner, Mode => AlignmentMode}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{unreachable, _}
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger}
import htsjdk.samtools.util.SequenceUtil

import scala.collection.immutable


@clp(group=ClpGroups.SamOrBam, description=
  """
    |Identifies primers that generate the reads in the input BAM file. For sequence data that is generated by targeted
    |PCR amplification, will identify the most likely primer for both the R1 and R2 ends of each template. Will detect
    |canonical (i.e. expected) primer pairings, unexpected primer pairings and primer-dimer products.
    |
    |## Primers
    |
    |The input primers file must be tab-separated, and contain a header, and then one row per primer:
    |
    |  * `pair_id` - the unique identifier for the primer pair
    |  * `primer_id` - the unique identifier for the primer in a primer pair
    |  * `sequence` - the DNA sequence of the primer as ordered, including degenerate bases.
    |  * `ref_name` - the chromosome/contig
    |  * `start` - the start position (1-based)
    |  * `end` - the end position (1-based, inclusive)
    |  * `strand` - F/f/+ if the forward primer, R/r/- otherwise
    |
    |If the primer is on the reverse strand, the sequence should be the reverse complement of the forward genomic strand.
    |Typically, the forward primers will match the positive genomic strand and reverse primers will match the reverse
    |complement of the positive genomic strand.
    |
    |The primer file must contain headers, e.g:
    |
    |```
    |pair_id primer_id sequence chr  start   end     strand
    |1       1_1       GATTACA  chr1 1010873 1010894 +
    |1       1_2       ACATTAG  chr1 1011118 1011137 -
    |```
    |
    |Primers without mapping coordinates should have an empty `ref_name`.  The strand should still be given.
    |
    |Use the same `pair_id` to have `N` forward and `M` primers as part of a primer pair set (where either `N > 1` or
    |`M > 1).  This is important for determing if both reads of a pair match the expected primer pair set (i.e. called
    |"canonical", otherwise "non-canonical", "self-dimer", or "dimer"; see below).  The `--multi-primer-pairs` options
    |must be set for this case.
    |
    |## Primer Matching
    |
    |Generally, the alignment position of a read is first used to find primers that overlap.  Next, the primer and read
    |sequence are compared based on mismatch-only (ungapped) alignment, then gapped alignment if necessary.
    |
    |When matching based on location or based on mismatch-only (ungapped) alignment, the primer sequence compared
    |to the read sequence must have at most the maximum mismatch rate to be accepted as a match
    |(see `--max-mismatch-rate`).  For gapped alignment, the alignment score rate (i.e. alignment score divided by the
    |primer length) must be at least `--min-alignment-score-rate` to be accepted as a match.
    |
    |The `--ignore-primer-strand` can be used if to which strand the primers belong is unknown.  This may the case for
    |primers designed against an assembly with no canonical top or bottom strand.
    |
    |### Paired-End Matching (5' only)
    |
    |If the input are paired-end reads, then the reads are further classified as follows for each read pair:
    |
    |1. If both R1 and R2 have primer matches on the 5' ends, then classify as:
    |  a. Self-dimer if the primers are the same.
    |  b. Canonical if the primers from the same primer pair set (i.e. same `pair_id`) but different primers.
    |  c. Cross-dimer if the primers are from different primer pair sets (i.e. different `pair_id`) and either the p
    |     primers map to different chromosomes or the inferred insert length (based on the span of the primer pair) is
    |     too small or too large (see the `--min-insert-length` and `--max-insert-length` options).
    |  d. Non-canonical otherwise.
    |2. Otherwise, if only one read (R1 or R2) has a primer match, classify as a single primer match.
    |3. Otherwise, classify the reads as matching no primers.
    |
    |## Matching Primers on the 3' End
    |
    |The `--three-prime` option can be used to also search the 3' end of every read for a primer sequence as follows:
    |
    |1. If a primer pair is identified using the 5' end, the expected paired primer is searched for at the 3' end.
    |2. Otherwise, assign a primer based on gapped alignment.
    |
    |When a primer match is found at both 5' and 3' ends of reads, the tool enforces that the 5' end of R1 and the 3'
    |end of R2 match the same primer, and vice versa.
    |
    |### Unmapped Data or Primers without Mappings
    |
    |If no reference is given or mapping information is not available, matching using mapping locations is skipped.
    |
    |### Speeding Up Gapped Alignment
    |
    |Install [the ksw executable](https://github.com/nh13/ksw) manually or via conda (`conda install -c bioconda ksw`)
    |and supply the path to the `ksw` executable via `--ksw`.  This may help when many gapped alignments need to be
    |performed.
    |
    |## SAM Tags
    |
    |Each read will be annotated with SAM tags based on results of primer matching.
    |
    |The `pp` tag stores the type of match found for the read pair.  Values include: 'Canonical', 'SelfDimer',
    |'CrossDimer', 'NonCanonical', 'Single', or 'NoMatch'.  See the Paired-End Matching section for a description of
    |these values.
    |
    |The `f5`, `r5` describe the primer match for the forward and reverse strand primers matching the 5' end of the read
    |respectively.  The `f3`, `r3` describe the primer match for the forward and reverse strand primers matching the 3'
    |end of the read respectively (when `--three-prime` is used).  If no match was found, then the tag is set "none",
    |otherwise, a comma-delimited set of values are given as follows
    |
    |  * `pair_id` - the unique identifier for the primer pair
    |  * `primer_id` - the unique identifier for the primer in a primer pair
    |  * `<ref_name>:<start>-<end>` - the primer chromosome/contig, start (1-based), and end (1-based, inclusive)
    |  * `strand` - '+' if the forward primer, '-' otherwise
    |  * `read_num` - the read number of the match (1 or 2).
    |  * `match_offset` - the offset in the read where the primer matches
    |  * `match_length`- the length of the primer match in the read
    |  * `match_type` - how the match was found; valid values are 'location', 'gapped', or 'ungapped'
    |  * `match_type_info` - additional information based on the `match_type`, containing one or more values
    |
    |The `match_type_info` for each `match_type` is as follow:
    |  * location:
    |    * `num_mismatches` - the number of mismatches between the primer and read sequence
    |  * ungapped:
    |    * `num_mismatches` - the number of mismatches between the primer and read sequence
    |    * 'next_best` - the number of mismatches in the next best ungapped alignment (or "na" if none was found)
    |  * gapped:
    |    * `score` - the alignment score
    |    * 'next_best` - the alignment score of then next best alignment
  """)
class IdentifyPrimers
(@arg(flag='i', doc="Input BAM file.")  val input: PathToBam,
 @arg(flag='p', doc="File containing information about the primer pairs.") val primerPairs: FilePath,
 @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Path prefix for the metrics files.") val metrics: PathPrefix,
 @arg(flag='S', doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
 @arg(flag='3', doc="Search for primers on the 3' end of each read as well (5' only by default).") val threePrime: Boolean = false,
 @arg(          doc="Maximum per-base mismatch rate for a primer to be considered a match.") val maxMismatchRate: Double = 0.05,
 @arg(          doc="The minimum per-base alignment score rate for gapped alignment for a primer to be considered a match. " +
                    "This is the total score divided by the primer length.") val minAlignmentScoreRate: Double = 0.5,
 @arg(          doc="The match score to use for aligning to primer sequences (must be >= 0).") val matchScore: Int = 1,
 @arg(          doc="The mismatch score to use for aligning to primer sequences (must be <= 0).") val mismatchScore: Int = -4,
 @arg(          doc="The gap open score to use for aligning to primer sequences (must be <= 0).") val gapOpen: Int = -6,
 @arg(          doc="The gap extension score to use for aligning to primer sequences (must be <= 0).") val gapExtend: Int = -1,
 @arg(flag='t', doc="The number of threads to use.") val threads: Int = 1,
 @arg(          doc="Skip gapped-alignment matching") val skipGappedAlignment: Boolean = false,
 @arg(          doc="Path to the ksw aligner.") val ksw: Option[String] = None,
 @arg(flag='k', doc="Skip ungapped alignment if no kmer of this length is common between any primer and a given read.") val ungappedKmerLength: Option[Int] = None,
 @arg(flag='K', doc="Skip gapped alignment if no kmer of this length is common between any primer and a given read.") val gappedKmerLength: Option[Int] = None,
 @arg(          doc="Allow multiple primers on the same strand to have the same `pair_id`.") val multiPrimerPairs: Boolean = false,
 @arg(          doc="Ignore the strand of the primer when matching.") val ignorePrimerStrand: Boolean = false,
 @arg(          doc="The minimum insert length for a non-canonical primer pair match, otherwise a dimer.") val minInsertLength: Int = 50,
 @arg(          doc="The maximum insert length for a non-canonical primer pair match, otherwise a dimer") val maxInsertLength: Int = 250
) extends FgBioTool with LazyLogging {

  /** The maximum number of templates in memory. */
  private val maxTemplatesInRam: Option[Int] = None

  /** The number of templates to process at a time per thread. */
  private val templatesPerThread: Int = 5000

  private val kswExecutable: Option[FilePath] = this.ksw.map(Paths.get(_)).flatMap {
    case p if Files.exists(p)                  => Some(p)
    case name if name.contains(File.separator) =>
      throw new ValidationException(s"Path to the ksw executable does not exist: $ksw")
    case name                                  =>
      // try the system path
      val path = System.getenv("PATH")
      validate(path != null, s"Searching for the ksw executable '$ksw' on the PATH, but the PATH environment variable was not defined.")
      path.split(File.pathSeparatorChar)
        .view
        .map(PathUtil.pathTo(_))
        .map(p => p.resolve(name)).find(ex => Files.exists(ex))
        .orElse {
          throw new ValidationException(s"Is the path to the ksw executable mis-typed? Could not find ksw executable ${ksw} in PATH: $path")
        }
  }

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  Io.assertCanWriteFile(metrics)
  kswExecutable.foreach(Io.assertReadable)

  validate(maxMismatchRate >= 0, "--max-mismatches must be >= 0")
  validate(matchScore >= 0,       "--match-score must be >= 0")
  validate(mismatchScore <= 0,    "--mismatch-score must be <= 0")
  validate(gapOpen <= 0,          "--gap-open must be <= 0")
  validate(gapExtend <= 0,        "--gap-extend must be <= 0")
  validate(minAlignmentScoreRate >= 0, "--min-alignment-score-rate must be >= 0")
  maxTemplatesInRam.foreach { m => validate(templatesPerThread < m, "--max-templates-in-ram must be greater than or equal to --templates-per-thread")}

  private val (locationBasedMatcher, ungappedBasedMatcher, gappedBasedMatcher, pairIdToPrimers) = {
    val primers  = Primer.read(this.primerPairs, multiPrimerPairs = multiPrimerPairs).toIndexedSeq
    val location = new LocationBasedPrimerMatcher(primers, slop, ignorePrimerStrand, maxMismatchRate)
    val ungapped = new UngappedAlignmentBasedPrimerMatcher(primers, slop, ignorePrimerStrand, maxMismatchRate, ungappedKmerLength)
    val gapped   = {
      val aligner: Aligner = Aligner(matchScore = matchScore, mismatchScore = mismatchScore, gapOpen = gapOpen, gapExtend = gapExtend, mode = AlignmentMode.Glocal)
      new GappedAlignmentBasedPrimerMatcher(primers, slop, ignorePrimerStrand, aligner, minAlignmentScoreRate, gappedKmerLength)
    }
    (location, ungapped, gapped, primers.groupBy(_.pair_id))
  }

  val numAlignments: AtomicLong = new AtomicLong(0)

  private val PrimerPairMatchTypeTag: String       = "pp"
  private val PrimerInfoForward5PrimeTag: String   = "f5"
  private val PrimerInfoReverse5PrimeTag: String   = "r5"
  private val PrimerInfoForward3PrimeTag: String   = "f3"
  private val PrimerInfoReverse3PrimeTag: String   = "r3"
  private val NoPrimerMatchInfo: String            = "none"
  private val comments: Seq[String]          = {
    val allTags = Seq(PrimerInfoForward5PrimeTag, PrimerInfoReverse5PrimeTag, PrimerInfoForward3PrimeTag, PrimerInfoReverse3PrimeTag).mkString("/")
    Seq(
      s"The $allTags tags store the primer match metadata for the forward and reverse strand respectively.",
      s"The $allTags tags are formatted as follow: <pair_id>,<primer_id>,<ref_name>:<start>-<end>,<strand>,<read-num>,<match-offset>,<match-length>,<match-type>,<match-type-info>.",
      s"The match-type is 'location', 'gapped', or 'ungapped' based on if the match was found using the location, mismatch-based (ungapped) alignment, or gapped-alignment.",
      s"The ${PrimerPairMatchTypeTag} tag is either 'canonical', 'self-dimer', 'cross-dimer', 'non-canonical', 'single', or '$NoPrimerMatchInfo', based on how the primers for pairs match."
    )
  }

  override def execute(): Unit = {
    this.kswExecutable match {
      case Some(p) => logger.info(s"Using ksw aligner: $p")
      case None    => logger.info("Using the internal aligner; install ksw for faster operation.")
    }
    val metricCounter = new SimpleCounter[TemplateTypeMetric]()

    // Group the reads by template
    val in: SamSource = SamSource(this.input)
    val iterator: Iterator[Template] = Bams.templateIterator(in)

    // NB: Add comments explaining tags in the output writer
    val (out: SamWriter, programGroupId: Option[String]) = {
      val header = {
        val h = in.header.clone()
        comments.foreach { comment => h.addComment(comment) }
        h
      }
      val pgId   = this.toolInfo.map { info => info.applyTo(header) }
      val writer = SamWriter(path = output, header = header, sort = SamOrder(in.header))
      (writer, pgId)
    }

    // We take the input records, and batch them into `majorBatchSize` number of records.  For each major-batch, we
    // split those into `templatesPerThread` sub-batches.  Each sub-batch is processed by a single thread.  We process
    // each major-batch serially, meaning we wait for one major-batch to complete before moving onto the next one.  This
    // is so we don't have to read all records into memory to parallize.  We set the `majorBatchSize` to have more
    // sub-batches than just `templatesPerThread * threads` so that we can more efficiently utlize the available threads.
    val majorBatchSize  = templatesPerThread * threads * 4
    val readingProgress = ProgressLogger(this.logger, verb="read", unit=5e4.toInt)

    // Batch templates to process them in individual threads.
    val outputIterator: Iterator[Seq[Template]] = if (threads > 1) {
      logger.info(f"Batching $majorBatchSize%,d templates with $templatesPerThread%,d templates per thread.")

      import com.fulcrumgenomics.commons.CommonsDef.ParSupport
      // Developer Note: Iterator does not support parallel operations, so we need to group together records into a
      // [[List]] or [[Seq]].  A fixed number of records are grouped to reduce memory overhead.
      iterator
        .grouped(majorBatchSize)
        .flatMap { templates =>
          templates
            .grouped(templatesPerThread)
            .toStream
            .parWith(threads, fifo = false)
            .map { templates => processBatch(templates, metricCounter, readingProgress) }
            .seq // ensures that the only parallelism is processBatch
            .toIterator
        }
    }
    else {
      val aligner = newAligner
      logger.info(f"Batching $templatesPerThread%,d templates.")
      val iter = iterator
        .grouped(templatesPerThread)
        .map { templates => processBatch(templates, metricCounter, readingProgress, Some(aligner)) }
      new SelfClosingIterator[Seq[Template]](iter, () => aligner.close())
    }

    // Write the results
    val writingProgress = ProgressLogger(this.logger, "written", unit=1000000)
    outputIterator
      .flatMap(_.flatMap(_.allReads))
      .foreach { rec =>
        programGroupId.foreach { pgId => rec("PG") = pgId }
        writingProgress.record(rec)
        out += rec
      }

    val rate = numAlignments.get() / readingProgress.getElapsedSeconds.toDouble
    logger.info(f"Performed ${numAlignments.get()}%,d gapped alignments in total ($rate%,.2f alignments/second).")
    logger.info(f"Wrote ${readingProgress.getCount}%,d records.")

    in.safelyClose()
    out.close()

    // Create the [[TemplateTypeMetric]]
    val templateTypeMetrics = TemplateTypeMetric.metricsFrom(metricCounter)

    // Write metrics
    // Detailed metrics
    Metric.write(PathUtil.pathTo(this.metrics + ".detailed.txt"), templateTypeMetrics)
    // Summary metrics
    Metric.write(PathUtil.pathTo(this.metrics + ".summary.txt"), IdentifyPrimersMetric(templateTypeMetrics))
  }

  /** Creates a new [[BatchAligner]]. */
  private def newAligner: BatchAligner[Primer, SamRecordAlignable] = {
    BatchAligner(matchScore, mismatchScore, gapOpen, gapExtend, AlignmentMode.Glocal, this.kswExecutable)
  }

  // NB: batching alignment inputs to the aligner (i.e. ksw) is empirically faster than given them all or just one at a time.
  private def getAlignmentResults[T <: Alignable](alignmentTasks: Seq[AlignmentTask[Primer, T]],
                                  aligner: BatchAligner[Primer, T]): Option[PrimerMatch] = {
    case class PrimerAndScore(primer: Primer, score: Int)

    // Get the alignment results
    val alignmentResults: Seq[PrimerAndScore] = alignmentTasks.toIterator.zip(aligner.iterator)
      // Keep results that meet the minimum score and align with a minimum # of query bases
      .filter { case (alignmentInput, alignmentResult) =>
        val minQueryEnd = alignmentInput.queryLength - slop
        val alignmentScoreRate = alignmentResult.score * alignmentInput.query.length.toDouble
        alignmentScoreRate >= minAlignmentScoreRate && alignmentResult.queryEnd >= minQueryEnd
      }
      // now just keep the primer and score
      .map { case (alignmentInput, alignmentResult) => PrimerAndScore(alignmentInput.query, alignmentResult.score) }
      // get the best two alignment results by score
      .maxNBy(n = 2, _.score)

    // Create a primer match if alignments were found (i.e. a Option[PrimerMatch])
    alignmentResults match {
      case Seq()               => None
      case Seq(best)           =>
        val secondBestScore = (minAlignmentScoreRate * best.primer.length).toInt
        Some(GappedAlignmentPrimerMatch(primer = best.primer, score = best.score, secondBestScore = secondBestScore))
      case Seq(best, nextBest) => Some(GappedAlignmentPrimerMatch(primer = best.primer, score = best.score, secondBestScore = nextBest.score))
      case _                   => unreachable("Should have returned at most two items.")
    }
  }

  /** Matches the read based on location, then ungapped alignment.  If no match was found, return a list of alignment tasks. */
  private def toPrimerMatchOrAlignmentTasks(rec: SamRecord): Either[PrimerMatch, Seq[AlignmentTask[Primer, SamRecordAlignable]]] = {
    locationBasedMatcher.find(rec).orElse { ungappedBasedMatcher.find(rec) } match {
      case Some(pm) => Left(pm)
      case None     => Right(gappedBasedMatcher.toAlignmentTasks(rec))
    }
  }

  /** Processes a single batch of templates. */
  def processBatch(templates: Seq[Template],
                   metricCounter: SimpleCounter[TemplateTypeMetric],
                   progress: ProgressLogger,
                   batchAligner: Option[BatchAligner[Primer, SamRecordAlignable]] = None): Seq[Template] = {
    val counter = new SimpleCounter[TemplateTypeMetric]()
    val aligner = batchAligner.getOrElse(newAligner)
    val threePrimerAligner: BatchAligner[Primer, Alignable] = {
      // NB: we run this in local mode to get partial matches
      BatchAligner(matchScore, mismatchScore, gapOpen, gapExtend, AlignmentMode.Local, this.kswExecutable)
    }

    // reset the provided aligner
    batchAligner.foreach(_.reset())

    val templateMatchOptions = templates.toIterator.map { template =>
      // match based on location, then ungapped alignment, and if no match was found, then return the alignment tasks
      val r1MatchOrTasks = template.r1.map { r => toPrimerMatchOrAlignmentTasks(r) }
      val r2MatchOrTasks = template.r2.map { r => toPrimerMatchOrAlignmentTasks(r) }

      // NB: if we have a match for one end only of a pair, we could add that to the list of alignment tasks to perform
      // to ensure we detect self-dimers properly.  If the template is really short, we may not find the primer on the
      // 5' end side, as there may be sequence preceding it.

      // add any alignment tasks
      Seq(r1MatchOrTasks, r2MatchOrTasks).flatten.foreach {
        case Right(alignmentTasks) => alignmentTasks.foreach(t => aligner.append(t))
        case Left(_)               => Unit
      }

      (template, r1MatchOrTasks, r2MatchOrTasks)
    }

    // run through the match options, retrieving any alignment tasks that were performed
    templateMatchOptions.foreach { case (template, r1MatchOrTasks, r2MatchOrTasks) =>
      // get the final primer match, if any.  We must retrieve results from the aligner if we sent tasks
      val r1Match = r1MatchOrTasks.flatMap {
        case Left(pm)              => Some(pm)
        case Right(alignmentTasks) => getAlignmentResults(alignmentTasks, aligner)
      }
      val r2Match = r2MatchOrTasks.flatMap {
        case Left(pm)              => Some(pm)
        case Right(alignmentTasks) => getAlignmentResults(alignmentTasks, aligner)
      }

      // add the three prime matching
      val r1ThreePrimeMatch = template.r1.flatMap { r1 => matchThreePrime(r1, r1Match, r2Match, threePrimerAligner) }
      val r2ThreePrimeMatch = template.r2.flatMap { r2 => matchThreePrime(r2, r2Match, r1Match, threePrimerAligner) }

      // Get information about the matches
      val templateTypes = (template.r1, template.r2) match {
        case (Some(r1), Some(r2)) =>
          val rType = TemplateType(r1, Some(r2))
          val mType = formatPair(r1, r2, r1Match, r2Match, r1ThreePrimeMatch, r2ThreePrimeMatch)
          TemplateTypeMetric(rType, r1.isFrPair, mType, r1Match, r2Match)
        case (Some(r1), None) =>
          require(!r1.paired, s"Found paired read but missing R2 for ${r1.name}")
          val rType = TemplateType(r1, None)
          val mType = formatFragment(r1, r1Match, r1ThreePrimeMatch)
          TemplateTypeMetric(rType, r1.isFrPair, mType, r1Match, r2Match)
        case _ =>
          throw new IllegalStateException(s"Template did not have an R1: ${template.name}")
      }

      counter.count(templateTypes)
    }

    require(aligner.numAdded == aligner.numRetrieved, s"added: ${aligner.numAdded} alignments: ${aligner.numRetrieved}")
    require(aligner.numAvailable == 0, s"Found alignments available: ${aligner.numAvailable}")

    if (batchAligner.isEmpty) aligner.close()

    numAlignments.addAndGet(aligner.numRetrieved)
    metricCounter.synchronized {
      metricCounter += counter
      templates.flatMap(_.allReads).foreach { rec =>
        if (progress.record(rec)) {
          val rate = numAlignments.get() / progress.getElapsedSeconds.toDouble
          logger.info(f"Performed ${numAlignments.get()}%,d gapped alignments so far ($rate%,.2f alignments/second).")
        }
      }
    }

    templates
  }

  private def matchThreePrime(rec: SamRecord,
                              fivePrimeMatch: Option[PrimerMatch],
                              otherReadFivePrimerMatch: Option[PrimerMatch],
                              aligner: BatchAligner[Primer, Alignable]): Option[PrimerMatch] = if (!threePrime) None else {
    import com.fulcrumgenomics.alignment.Alignable

    // Get the primers to check.  If the read was paired, then use the primer found for the mate.  Otherwise, if no 5'
    // match, use all the primers, otherwise, get the primers from the same primer pair set and return the ones on the
    // opposite strand.
    val primersToCheck = (fivePrimeMatch, otherReadFivePrimerMatch) match {
      case (_, Some(pm)) => Seq(pm.primer)
      case (None,     _) => this.gappedBasedMatcher.primers
      case (Some(pm), _) => pairIdToPrimers(pm.primer.pair_id).filter(_.positive_strand == pm.primer.negativeStrand)
    }

    // create the alignment tasks
    val alignmentTasks = {
      val alignableCache = scala.collection.mutable.HashMap[(Int,Int), Alignable]()
      // for 3' matching, it is the same as 5' matching, just that we match the end of the reads for positive strand
      // primers (or unmapped) and start of the reads for negative strand primers
      primersToCheck.map { primer =>
        val (targetOffset, targetLength) = {
          if (primer.negativeStrand) (0, math.min(primer.sequence.length + slop, rec.length))
          else {
            val offset = math.max(0, rec.length - primer.length - slop)
            val length = rec.length - offset
            (offset, length)
          }
        }
        val target       = alignableCache.getOrElseUpdate((targetOffset, targetLength), SamRecordAlignable(rec, targetOffset, targetLength))
        AlignmentTask(query = primer, target = target)
      }
    }

    // submit them
    alignmentTasks.foreach { task => aligner.append(task)}

    // get the best match, if any
    getAlignmentResults(alignmentTasks, aligner)
  }

  private case class PrimerMatches(rec: SamRecord, fivePrimeMatch: Option[PrimerMatch], threePrimeMatch: Option[PrimerMatch], positiveStrand: Boolean)

  /** Adds tags to the records based on the primer matching results. */
  private def formatPair(r1: SamRecord,
                         r2: SamRecord,
                         r1FivePrimeMatch: Option[PrimerMatch],
                         r2FivePrimeMatch: Option[PrimerMatch],
                         r1ThreePrimeMatch: Option[PrimerMatch],
                         r2ThreePrimeMatch: Option[PrimerMatch]): PrimerPairMatchType = {
    // Get which rec is on the forward strand, and which is reverse.  First looks at the primer match for r1, then the
    // primer match for r2, otherwise r1 is forward.
    val (forwardPrimerMatch, reversePrimerMatch) = {
      val r1PrimerMatch = PrimerMatches(r1, r1FivePrimeMatch, r1ThreePrimeMatch, positiveStrand=true)
      val r2PrimerMatch = PrimerMatches(r2, r2FivePrimeMatch, r2ThreePrimeMatch, positiveStrand=false)
      (r1FivePrimeMatch, r2FivePrimeMatch) match {
        case (Some(m), _)    =>
          if (m.primer.positive_strand) (r1PrimerMatch, r2PrimerMatch) else  (r2PrimerMatch.copy(positiveStrand=true), r1PrimerMatch.copy(positiveStrand=false))
        case (None, Some(m)) =>
          if (m.primer.negativeStrand) (r1PrimerMatch, r2PrimerMatch) else  (r2PrimerMatch.copy(positiveStrand=true), r1PrimerMatch.copy(positiveStrand=false))
        case _               => (r1PrimerMatch, r2PrimerMatch)
      }
    }
    require(forwardPrimerMatch.positiveStrand && !reversePrimerMatch.positiveStrand)

    // Set the primer pair match type
    val matchType: PrimerPairMatchType = PrimerPairMatchType(forwardPrimerMatch.fivePrimeMatch, reversePrimerMatch.fivePrimeMatch, minInsertLength, maxInsertLength)
    Seq(r1, r2).foreach { rec => rec(PrimerPairMatchTypeTag) = matchType.toString }

    // Set the primer match info to place in the forward/reverse primer match tags
    def tagit(fwdTag: String, revTag: String, f: PrimerMatches => Option[PrimerMatch]): Unit = {
      val forwardInfo = f(forwardPrimerMatch).map(_.info(forwardPrimerMatch.rec)).getOrElse(NoPrimerMatchInfo)
      val reverseInfo = f(reversePrimerMatch).map(_.info(reversePrimerMatch.rec)).getOrElse(NoPrimerMatchInfo)
      Seq(r1, r2).foreach { rec =>
        rec(fwdTag) = forwardInfo
        rec(revTag) = reverseInfo
      }
    }
    tagit(PrimerInfoForward5PrimeTag, PrimerInfoReverse5PrimeTag, _.fivePrimeMatch)
    if (threePrime) tagit(PrimerInfoForward3PrimeTag, PrimerInfoReverse3PrimeTag, _.threePrimeMatch)

    matchType
  }

  /** Adds tags to the record based on the primer matching results. */
  private def formatFragment(frag: SamRecord, fragFivePrimeMatch: Option[PrimerMatch],
                             fragThreePrimeMatch: Option[PrimerMatch]): PrimerPairMatchType = {
    val forwardInfo = fragFivePrimeMatch.map(_.info(frag)).getOrElse(NoPrimerMatchInfo)
    val matchType: PrimerPairMatchType = PrimerPairMatchType(fragFivePrimeMatch, None, minInsertLength, maxInsertLength)

    frag(PrimerPairMatchTypeTag)     = matchType.toString
    frag(PrimerInfoForward5PrimeTag) = forwardInfo
    frag(PrimerInfoReverse5PrimeTag) = NoPrimerMatchInfo

    if (threePrime) {
      frag(PrimerInfoForward5PrimeTag) = fragThreePrimeMatch.map(_.info(frag)).getOrElse(NoPrimerMatchInfo)
      frag(PrimerInfoReverse5PrimeTag) = NoPrimerMatchInfo
    }

    matchType
  }
}
