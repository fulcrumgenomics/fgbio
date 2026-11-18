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
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.collection.LeastRecentlyUsedCache
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.CorrectUmis._
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util._

import scala.collection.mutable


object CorrectUmis {
  /**
    * Holds information about the match of a UMI to an expected UMI.
    * @param matched true if the match is acceptable, false otherwise
    * @param umi the fixed UMI sequence that was the closest match
    * @param mismatches the number of mismatches between the UMI and the reported best matching fixed UMI
    *  */
  private[umi] case class UmiMatch(matched: Boolean, umi: String, mismatches: Int)

  object UmiMatch {
    def apply(matched: Boolean, umi: Array[Byte], mismatches: Int): UmiMatch = {
      UmiMatch(matched=matched, umi=new String(umi), mismatches)
    }
  }

  /**
    * Metrics produced by `CorrectUmis` regarding the correction of UMI sequences to a fixed set of known UMIs.
    *
    * @param umi The corrected UMI sequence, or all `N`s for UMIs that could not be matched.  There is one row of `N`s
    *            for each length of the fixed UMIs, holding the unmatched UMIs of that length.
    * @param total_matches The number of UMI sequences that matched/were corrected to this UMI.
    * @param perfect_matches The number of UMI sequences that were perfect matches to this UMI.
    * @param one_mismatch_matches The number of UMI sequences that matched with a single mismatch.
    * @param two_mismatch_matches The number of UMI sequences that matched with two mismatches.
    * @param other_matches The number of UMI sequences that matched with three or more mismatches.
    * @param fraction_of_matches The fraction of all UMIs that matched or were corrected to this UMI, or zero if no UMIs
    *                            were examined.
    * @param representation The `total_matches` for this UMI divided by the _mean_ `total_matches` for all fixed UMIs
    *                       (excluding the rows of `N`s), or zero if no UMIs matched a fixed UMI.
    */
  case class UmiCorrectionMetrics(umi: String,
                                  var total_matches: Count = 0,
                                  var perfect_matches: Count = 0,
                                  var one_mismatch_matches: Count = 0,
                                  var two_mismatch_matches: Count = 0,
                                  var other_matches: Count = 0,
                                  var fraction_of_matches: Proportion = 0,
                                  var representation: Double = 0
                                 ) extends Metric

  /**
    * Finds pairs of UMIs within the given set that are within some edit distance of each other.
    *
    * @param umis the set of umis to check
    * @param distance the maximum edit distance at which to report pairs of UMIs
    * @return a Seq of three-tuples containing (umi1, umi2, edit_distance)
    */
  def findUmiPairsWithinDistance(umis: Seq[String], distance: Int): Seq[(String,String,Int)] = {
    val builder = Seq.newBuilder[(String, String, Int)]
    forloop (from=0, until=umis.length) { i =>
      forloop (from=i+1, until=umis.length) { j =>
        val d = Sequences.countMismatchesWithMax(umis(i).getBytes, umis(j).getBytes, distance+1)
        if (d <= distance) builder += ((umis(i), umis(j), d))
      }
    }
    builder.result()
  }
}

@clp(group=ClpGroups.Umi, description=
  """
    |Corrects UMIs stored in BAM files when a set of fixed UMIs is in use.  If the set of UMIs used in
    |an experiment is known and is a _subset_ of the possible randomers of the same length, it is possible
    |to error-correct UMIs prior to grouping reads by UMI.  This tool takes an input BAM with UMIs in a
    |tag (`RX` by default) and set of known UMIs (either on the command line or in a file) and produces:
    |
    |  1. A new BAM with corrected UMIs in the same tag the UMIs were found in
    |  2. Optionally a set of metrics about the representation of each UMI in the set
    |  3. Optionally a second BAM file of reads whose UMIs could not be corrected within the specific parameters
    |
    |The fixed UMIs may have different lengths.  Each UMI in the BAM file is compared only to the fixed UMIs of the
    |same length (e.g. a 6bp UMI is never corrected to a 7bp fixed UMI), and a record is rejected if any of its UMIs has
    |a length that matches no fixed UMI.  The metrics have one row of `N`s for each fixed UMI length, counting the UMIs
    |of that length that could not be matched.
    |
    |The UMIs in the BAM file must already be trimmed to their true lengths.  For example, if 7 bases were extracted
    |for every read but some of the fixed UMIs are 6bp, the 6bp UMIs (with an extra base) are compared only to the 7bp
    |fixed UMIs.  A warning is logged at the end for each fixed UMI length that no UMI in the BAM file had.
    |
    |Multiple UMIs that are concatenated with hyphens (e.g. `AACCAGT-AGGTAGA`) are split apart, corrected individually
    |and then re-assembled.  A read is accepted only if all the UMIs can be corrected.
    |
    |Correction is controlled by two parameters that are applied per-UMI:
    |
    |  1. _--max-mismatches_ controls how many mismatches (no-calls are counted as mismatches) are tolerated
    |         between a UMI as read and a fixed UMI.
    |  2. _--min-distance_ controls how many *more* mismatches the next best hit must have; it is not applied when
    |     only one fixed UMI has the same length as the UMI
    |
    |For example, with two fixed UMIs `AAAAA` and `CCCCC` and `--max-mismatches=3` and `--min-distance=2` the
    |following would happen:
    |
    |  - AAAAA would match to AAAAA
    |  - AAGTG would match to AAAAA with three mismatches because CCCCC has five mismatches and 5 >= 3 + 2
    |  - AACCA would be rejected because it is 2 mismatches to AAAAA and 3 to CCCCC and 3 < 2 + 2
    |
    |The set of fixed UMIs may be specified on the command line using `--umis umi1 umi2 ...` or via one or
    |more files of UMIs with a single sequence per line using `--umi-files umis.txt more_umis.txt`.  Each fixed UMI
    |may contain only the bases `A`, `C`, `G`, `T` and `N` (in either case), and may not be all `N`s.  If there
    |are multiple UMIs per template, leading to hyphenated UMI tags, the values for the fixed UMIs should
    |be single, non-hyphenated UMIs (e.g. if a record has `RX:Z:ACGT-GGCA`, you would use `--umis ACGT GGCA`).
    |
    |Records which have their UMIs corrected (i.e. the UMI is not identical to one of the expected UMIs but is
    |close enough to be corrected) will by default have their original UMI stored in the `OX` tag. This can be
    |disabled with the `--dont-store-original-umis` option.
    |
    |For a large number of input UMIs, the `--cache-size` option may used to speed up the tool.  To disable
    |using a cache, set the value to `0`.
    |
    |The reverse complement (using `--revcomp`) option will reverse complement the UMI in place.  In the case of
    |multiple UMIs concatenated together, the individual UMIs are reverse complemented and the order
    |reversed (eg. `AAGG-ACTG` is changed to `CAGT-CCTT`).
  """)
class CorrectUmis
(@arg(flag='i', doc="Input SAM or BAM file.")  val input: PathToBam,
 @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
 @arg(flag='r', doc="Reject BAM file to save unassigned reads.") val rejects: Option[PathToBam] = None,
 @arg(flag='M', doc="Metrics file to write.") val metrics: Option[FilePath] = None,
 @arg(flag='m', doc="Maximum number of mismatches between a UMI and an expected UMI.") val maxMismatches: Int,
 @arg(flag='d', name="min-distance", doc="Minimum difference (of mismatch distance) to next-best UMI.") val minDistanceDiff: Int,
 @arg(flag='u', doc="Expected UMI sequences.", minElements=0) val umis: Seq[String] = Seq.empty,
 @arg(flag='U', doc="File of UMI sequences, one per line.", minElements=0) val umiFiles: Seq[FilePath] = Seq.empty,
 @arg(flag='t', doc="Tag in which UMIs are stored.") val umiTag: String = ConsensusTags.UmiBases,
 @arg(flag='x', doc="Don't store original UMIs upon correction.") val dontStoreOriginalUmis: Boolean = false,
 @arg(doc="The number of uncorrected UMIs to cache; zero will disable the cache.") val cacheSize: Int = 100000,
 @arg(doc="The minimum ratio of kept UMIs to accept. A ratio below this will cause a failure (but all files will still be written).  Not checked if the input has no records.") val minCorrected: Option[Double] = None,
 @arg(doc="Reverse complement the UMIs in the BAM file prior to correcting.") val revcomp: Boolean = false
) extends FgBioTool with LazyLogging {

  validate(umis.nonEmpty || umiFiles.nonEmpty, "At least one UMI or UMI file must be provided.")
  Io.assertReadable(input)
  Io.assertReadable(umiFiles)
  Io.assertCanWriteFile(output)
  rejects.foreach(Io.assertCanWriteFile(_))
  minCorrected.foreach(m => validate(m >= 0 && m <= 1, "--min-corrected must be between 0 and 1."))

  // Construct the cache
  private lazy val cache = new LeastRecentlyUsedCache[String,UmiMatch](maxEntries = cacheSize)

  override def execute(): Unit = {
    // Construct the full set of UMI sequences to match against, grouped by length
    val umiSequencesByLength: Map[Int, Array[String]] = {
      val set = mutable.HashSet[String](umis.map(_.toUpperCase):_*)
      umiFiles.foreach(Io.readLines(_).map(_.trim.toUpperCase).filter(_.nonEmpty).foreach(set.add))
      validate(set.nonEmpty, s"At least one UMI sequence must be provided; none found in files ${umiFiles.mkString(", ")}")
      val invalid = set.filterNot(_.matches("[ACGTN]+")).toSeq.sorted
      validate(invalid.isEmpty, s"Fixed UMIs may contain only A, C, G, T and N (and may not be hyphenated): ${invalid.mkString(", ")}")
      // An all-N fixed UMI matches nothing meaningful, and would share its metrics row with the unmatched UMIs
      val allNs = set.filter(_.forall(_ == 'N')).toSeq.sorted
      validate(allNs.isEmpty, s"Fixed UMIs may not be all Ns: ${allNs.mkString(", ")}")
      set.toArray.groupBy(_.length)
    }
    val umiBytesByLength: Map[Int, Array[Array[Byte]]] = umiSequencesByLength.view.mapValues(_.map(_.getBytes)).toMap
    val expectedUmiLengths = umiSequencesByLength.keys.toSeq.sorted.mkString(", ")

    // Warn if any of the UMIs are too close together
    umiSequencesByLength.values.foreach { umiSequences =>
      CorrectUmis.findUmiPairsWithinDistance(umiSequences.toSeq, minDistanceDiff - 1).foreach { case (umi1, umi2, distance) =>
        logger.warning(s"Umis $umi1 and $umi2 are $distance edits apart which is less than the min distance diff: $minDistanceDiff")
      }
    }

    // Construct the UMI metrics objects, with one all-N "unmatched" UMI per fixed UMI length
    val unmatchedUmiByLength = umiSequencesByLength.keys.map(umiLength => umiLength -> "N" * umiLength).toMap
    val umiMetrics           = (umiSequencesByLength.values.flatten ++ unmatchedUmiByLength.values).map(umi => umi -> UmiCorrectionMetrics(umi)).toMap

    // Now go through and correct the UMIs in the BAM
    var (totalRecords, missingUmisRecords, wrongLengthRecords, mismatchedRecords) = (0L, 0L, 0L, 0L)
    val observedUmiLengths = mutable.HashSet[Int]()
    val in        = SamSource(input)
    val out       = SamWriter(output, in.header)
    val rejectOut = rejects.map(r => SamWriter(r, in.header))
    val progress  = ProgressLogger(logger)

    in.foreach { rec =>
      totalRecords += 1

      rec.get[String](umiTag) match {
        case None | Some("") =>
          if (missingUmisRecords == 0) logger.warning(s"Read (${rec.name}) detected without UMI in tag $umiTag")
          missingUmisRecords += 1
          rejectOut.foreach(w => w += rec)
        case Some(umi: String) =>
          val sequences = if (revcomp) umi.split("-", -1).map(Sequences.revcomp).reverse else umi.split("-", -1)
          sequences.foreach(seq => observedUmiLengths += seq.length)
          if (sequences.exists(seq => !umiSequencesByLength.contains(seq.length))) {
            if (wrongLengthRecords == 0) {
              logger.warning(s"Read (${rec.name}) detected with unexpected length UMI(s): ${sequences.mkString(" ")}.")
              logger.warning(s"Expected UMI length(s): $expectedUmiLengths")
            }
            wrongLengthRecords += 1
            rejectOut.foreach(w => w += rec)
          }
          else {
            // Find matches for all the UMIs
            val matches = sequences.map(seq => findBestMatch(seq.toUpperCase, umiBytesByLength(seq.length)))

            // Update the metrics
            matches.foreach { m =>
              if (m.matched) {
                val metric = umiMetrics(m.umi)
                metric.total_matches += 1
                m.mismatches match {
                  case 0 => metric.perfect_matches      += 1
                  case 1 => metric.one_mismatch_matches += 1
                  case 2 => metric.two_mismatch_matches += 1
                  case _ => metric.other_matches        += 1
                }
              }
              else {
                umiMetrics(unmatchedUmiByLength(m.umi.length)).total_matches += 1
              }
            }

            // Output the corrected read
            if (matches.forall(_.matched)) {
              // Store the original UMI if enabled and there are mismatches
              if (!dontStoreOriginalUmis && !matches.forall(_.mismatches == 0)) {
                rec(ConsensusTags.OriginalUmiBases) = rec[String](this.umiTag)
              }

              val correctedUmi = matches.map(_.umi).mkString("-")
              rec(this.umiTag) = correctedUmi
              out += rec
            }
            else {
              mismatchedRecords += 1
              rejectOut.foreach(r => r += rec)
            }
          }
      }

      progress.record(rec)
    }

    in.safelyClose()
    out.close()
    rejectOut.foreach(_.close())

    // Finalize the metrics.  The denominators are zero when no UMIs were examined or none matched a fixed UMI, in
    // which case the fractions are reported as zero rather than NaN or Infinity.
    val sortedMetrics        = umiMetrics.values.toSeq.sortBy(_.umi)
    val fixedUmiMetrics      = sortedMetrics.filterNot(m => unmatchedUmiByLength.get(m.umi.length).contains(m.umi))
    val totalWithUnmatched   = sortedMetrics.map(_.total_matches).sum.toDouble
    val meanWithoutUnmatched = fixedUmiMetrics.map(_.total_matches).sum / fixedUmiMetrics.length.toDouble
    sortedMetrics.foreach { m =>
      m.fraction_of_matches = if (totalWithUnmatched > 0)   m.total_matches / totalWithUnmatched   else 0
      m.representation      = if (meanWithoutUnmatched > 0) m.total_matches / meanWithoutUnmatched else 0
    }

    this.metrics.foreach(path => Metric.write(path, sortedMetrics))

    // Log some summary info & warnings
    val discarded = missingUmisRecords + wrongLengthRecords + mismatchedRecords
    val kept      = totalRecords - discarded
    logger.info(f"Read ${totalRecords}%,d; kept ${kept}%,d and rejected ${discarded}%,d")

    if (missingUmisRecords > 0 || wrongLengthRecords > 0) {
      logger.error("###################################################################")
      if (missingUmisRecords > 0) logger.error(s"# ${missingUmisRecords} were missing UMI attributes in the BAM file!")
      if (wrongLengthRecords > 0) logger.error(s"# ${wrongLengthRecords} had UMIs whose lengths matched no fixed UMI length in the BAM file!")
      logger.error("###################################################################")
    }

    // A fixed UMI length that no UMI had may be a malformed fixed UMI, or UMIs that were not trimmed to their true length
    umiSequencesByLength.keys.filterNot(observedUmiLengths.contains).toSeq.sorted.foreach { umiLength =>
      val fixed = umiSequencesByLength(umiLength).sorted.mkString(", ")
      logger.warning(s"No UMI in the input had length $umiLength, the length of fixed UMI(s): $fixed")
    }

    minCorrected.foreach { min =>
      if (totalRecords == 0) logger.warning("No records were read; skipping the --min-corrected check.")
      else {
        val ratioKept = 1.0 * kept / totalRecords
        assert(ratioKept >= min,
          f"# Final ratio of reads kept / total was ${ratioKept}%2.2f (user specified minimum was ${min}%2.2f) " +
            "This could indicate a mismatch between library preparation and the provided UMI file."
        )
      }
    }
  }

  /** Given a UMI sequence and a set of fixed UMIs, report the best match.
   *
   * Both bases and all UMIs in umis must be normalized to upper case prior to calling this method.
   *
   * @param bases the input UMI sequence to match
   * @param umis the set of fixed UMI sequences to compare against
   * @return a UmiMatch indicating whether a match was found, the best matching UMI, and the number of mismatches
   * @throws ArrayIndexOutOfBoundsException if `bases.length` differs from the length of any UMI in `umis`
   * @note '''Preconditions (caller must ensure):'''
   *       - all UMIs in `umis` must have the same length as `bases` (required - undefined behavior if violated)
   *       - `umis` must be non-empty (required - undefined behavior if violated)
   * */
  private[umi] def findBestMatch(bases: String, umis: Array[Array[Byte]]): UmiMatch = {
    val cachedResult = if (this.cacheSize > 0) cache.get(bases) else None
    cachedResult match {
      case Some(result) => result
      case None         =>
        var minMismatches      = bases.length + 1
        var minIndex           = -1
        var nextBestMismatches = bases.length + 1
        var umiIndex           = 0
        val bytes              = bases.getBytes
        while (umiIndex < umis.length) {
          val mismatches = Sequences.countMismatchesWithMax(bytes, umis(umiIndex), nextBestMismatches)
          if (mismatches < minMismatches) {
            nextBestMismatches = minMismatches
            minMismatches      = mismatches
            minIndex           = umiIndex
          } else if (mismatches < nextBestMismatches) {
            nextBestMismatches = mismatches
          }
          umiIndex += 1
        }
        // With a single candidate (e.g. the only fixed UMI of this length) there is no next-best UMI to be far from
        val distanceOk = umis.length == 1 || nextBestMismatches >= minMismatches + this.minDistanceDiff
        val matched    = (minMismatches <= maxMismatches) && distanceOk
        val umiMatch = UmiMatch(matched, umis(minIndex), minMismatches)
        if (cacheSize > 0) cache.put(bases, umiMatch)
        umiMatch
    }
  }
}
