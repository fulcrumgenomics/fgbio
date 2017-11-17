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
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.CorrectUmis._
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util.{Metric, _}

import scala.collection.mutable

object CorrectUmis {
  /**
    * Holds information about the match of a UMI to an expected UMI.
    * @param matched true if the match is acceptable, false otherwise
    * @param umi the fixed UMI sequence that was the closest match
    * @param mismatches the number of mismatches between the UMI and the reported best matching fixed UMI
    *  */
  private[umi] case class UmiMatch(matched: Boolean, umi: String, mismatches: Int)

  /**
    * Metrics produced by `CorrectUmis` regarding the correction of UMI sequences to a fixed set of known UMIs.
    *
    * @param umi The corrected UMI sequence (or all `N`s for unmatched).
    * @param total_matches The number of UMI sequences that matched/were corrected to this UMI.
    * @param perfect_matches The number of UMI sequences that were perfect matches to this UMI.
    * @param one_mismatch_matches The number of UMI sequences that matched with a single mismatch.
    * @param two_mismatch_matches The number of UMI sequences that matched with two mismatches.
    * @param other_matches The number of UMI sequences that matched with three or more mismatches.
    * @param fraction_of_matches The fraction of all UMIs that matched or were corrected to this UMI.
    * @param representation The `total_matches` for this UMI divided by the _mean_ `total_matches` for all UMIs.
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
    umis.tails.flatMap {
      case x +: ys => ys.map(y => (x, y, Sequences.countMismatches(x, y))).filter(_._3 <= distance)
      case _      => Seq.empty
    }.toList
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
    |All of the fixed UMIs must be of he same length, and all UMIs in the BAM file must also have the same
    |length.  Multiple UMIs that are concatenated with hyphens (e.g. `AACCAGT-AGGTAGA`) are split apart,
    |corrected individually and then re-assembled.  A read is accepted only if all the UMIs can be corrected.
    |
    |Correction is controlled by two parameters that are applied per-UMI:
    |
    |  1. _--max-mismatches_ controls how many mismatches (no-calls are counted as mismatches) are tolerated
    |         between a UMI as read and a fixed UMI.
    |  2. _--min-distance_ controls how many more mismatches the next best hit must have
    |
    |For example, with two fixed UMIs `AAAAA` and `CCCCC` and `--max-mismatches=3` and `--min-distance=2` the
    |following would happen:
    |
    |  - AAAAA would match to AAAAA
    |  - AAGTG would match to AAAAA with three mismatches because CCCCCC has six mismatches and 6 >= 3 + 2
    |  - AACCA would be rejected because it is 2 mismatches to AAAAA and 3 to CCCCCC and 3 <= 2 + 2
    |
    |The set of fixed UMIs may be specified on the command line using `--umis umi1 umi2 ...` or via one or
    |more files of UMIs with a single sequence per line using `--umi-files umis.txt more_umis.txt`.  If there
    |are multiple UMIs per template, leading to hyphenated UMI tags, the values for the fixed UMIs should
    |be single, non-hyphenated UMIs (e.g. if a record has `RX:Z:ACGT-GGCA`, you would use `--umis ACGT GGCA`).
    |
    |Records which have their UMIs corrected (i.e. the UMI is not identical to one of the expected UMIs but is
    |close enough to be corrected) will by default have their original UMI stored in the `OX` tag. This can be
    |disabled with the `--dont-store-original-umis` option.
  """)
class CorrectUmis
( @arg(flag='i', doc="Input SAM or BAM file.")  val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='r', doc="Reject BAM file to save unassigned reads.") val rejects: Option[PathToBam] = None,
  @arg(flag='M', doc="Metrics file to write.") val metrics: Option[FilePath] = None,
  @arg(flag='m', doc="Maximum number of mismatches between a UMI and an expected UMI.") val maxMismatches: Int,
  @arg(flag='d', doc="Minimum distance (in mismatches) to next best UMI.") val minDistance: Int,
  @arg(flag='u', doc="Expected UMI sequences.", minElements=0) val umis: Seq[String] = Seq.empty,
  @arg(flag='U', doc="File of UMI sequences, one per line.", minElements=0) val umiFiles: Seq[FilePath] = Seq.empty,
  @arg(flag='t', doc="Tag in which UMIs are stored.") val umiTag: String = ConsensusTags.UmiBases,
  @arg(flag='x', doc="Don't store original UMIs upon correction.") val dontStoreOriginalUmis: Boolean = false
) extends FgBioTool with LazyLogging {

  validate(umis.nonEmpty || umiFiles.nonEmpty, "At least one UMI or UMI file must be provided.")
  Io.assertReadable(input)
  Io.assertReadable(umiFiles)
  Io.assertCanWriteFile(output)
  rejects.foreach(Io.assertCanWriteFile(_))

  override def execute(): Unit = {
    // Construct the full set of UMI sequences to match again
    val (umiSequences, umiLength) = {
      val set = mutable.HashSet[String](umis:_*)
      umiFiles.foreach(Io.readLines(_).map(_.trim).filter(_.nonEmpty).foreach(set.add))
      validate(set.nonEmpty, s"At least one UMI sequence must be provided; none found in files ${umiFiles.mkString(", ")}")

      val lengths = set.map(_.length)
      validate(lengths.size == 1, s"UMIs of multiple lengths found. Lengths: ${lengths.mkString(", ")}")
      (set.toArray, lengths.head)
    }

    // Warn if any of the UMIs are too close together
    CorrectUmis.findUmiPairsWithinDistance(umiSequences, minDistance-1).foreach { case (umi1, umi2, distance) =>
        logger.warning(s"Umis $umi1 and $umi2 are $distance edits apart which is less than the min distance: $minDistance")
    }

    // Construct the UMI metrics objects
    val unmatchedUmi = "N" * umiLength
    val umiMetrics   = (umiSequences ++ Seq(unmatchedUmi)).map(umi => umi -> UmiCorrectionMetrics(umi)).toMap

    // Now go through and correct the UMIs in the BAM
    var (totalRecords, missingUmisRecords, wrongLengthRecords, mismatchedRecords) = (0L, 0L, 0L, 0L)
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
        case Some(umi) =>
          val sequences = umi.split('-')
          if (sequences.exists(_.length != umiLength)) {
            if (wrongLengthRecords == 0) logger.warning(s"Read (${rec.name}) detected with unexpected length UMI(s): ${sequences.mkString(" ")}")
            wrongLengthRecords += 1
            rejectOut.foreach(w => w += rec)
          }
          else {
            // Find matches for all the UMIs
            val matches = sequences.map(findBestMatch(_, umiSequences))

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
                umiMetrics(unmatchedUmi).total_matches += 1
              }
            }

            // Output the corrected read
            if (matches.forall(_.matched)) {
              // Store the original UMI if enabled and there are mismatches
              if (!dontStoreOriginalUmis && !matches.forall(_.mismatches == 0)) {
                rec(ConsensusTags.OriginalUmiBases) = rec(this.umiTag)
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

    // Finalize the metrics
    val sortedMetrics        = umiMetrics.values.toSeq.sortBy(_.umi)
    val totalWithUnmatched   = sortedMetrics.map(_.total_matches).sum.toDouble
    val meanWithoutUnmatched = sortedMetrics.filter(_.umi != unmatchedUmi).map(_.total_matches).sum / (sortedMetrics.size - 1d)
    sortedMetrics.foreach { m =>
      m.fraction_of_matches = m.total_matches / totalWithUnmatched
      m.representation      = m.total_matches / meanWithoutUnmatched
    }

    this.metrics.foreach(path => Metric.write(path, sortedMetrics))

    // Log some summary info & warnings
    val discarded = missingUmisRecords + wrongLengthRecords + mismatchedRecords
    val kept      = totalRecords - discarded
    logger.info(f"Read ${totalRecords}%,d; kept ${kept}%,d and rejected ${discarded}%,d")

    if (missingUmisRecords > 0 || wrongLengthRecords > 0) {
      logger.error("###################################################################")
      if (missingUmisRecords > 0) logger.error(s"# ${mismatchedRecords} were missing UMI attributes in the BAM file!")
      if (wrongLengthRecords > 0) logger.error(s"# ${wrongLengthRecords} had unexpected UMIs of differing lengths in the BAM file!")
      logger.error("###################################################################")
    }
  }

  /** Given a UMI sequence and a set of fixed UMIs, report the best match. */
  private[umi] def findBestMatch(bases: String, umis: Array[String]): UmiMatch = {
    val mismatchLimit = this.maxMismatches + minDistance
    val mismatches    = umis.map(umi => Sequences.countMismatches(bases, umi))
    val min           = mismatches.min
    val matched       = (min <= maxMismatches) && (mismatches.count(m => m < min + this.minDistance) == 1)

    UmiMatch(matched, umis(mismatches.indexOf(min)), min)
  }
}

