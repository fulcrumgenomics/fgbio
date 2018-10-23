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

import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.alignment.{Aligner, Mode}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.{SAMFileHeader, SAMUtils, SamPairUtil}
import org.scalatest.OptionValues

final class IdentifyPrimersTest extends UnitSpec with OptionValues {

  // Some basic settings for the tests
  private val primerLength = 11
  private val maxMismatches = 3
  private val minAlignmentScore = 5

  // Defaults for all tests
  private val slop                  = 1
  private val maxMismatcheRate      = maxMismatches.toDouble / primerLength
  private val minAlignmentScoreRate = minAlignmentScore.toDouble / primerLength

  // The default aligner
  private val matchScore    = 1
  private val mismatchScore = -3
  private val gapOpen       = -6
  private val gapExtend     = -1
  private val aligner       = Aligner(matchScore, mismatchScore, gapOpen, gapExtend, mode=Mode.Glocal)

  /** Companion object to [[AlignmentResult]] to help making objects with defaults. */
  private object AlignmentResult {
    def apply(mmScore: Int): AlignmentResult = new AlignmentResult(mmScore=Some(mmScore))
    def apply(mmScore: Int, mmNextScore: Int): AlignmentResult = new AlignmentResult(mmScore=Some(mmScore), mmNextScore = Some(mmNextScore))
    def apply(mmScore: Int, mmNextScore: Int, fullNextBest: Int): AlignmentResult = new AlignmentResult(mmScore = Some(mmScore), mmNextScore = Some(mmNextScore))
  }

  /** Stores the results of running [[UngappedAlignmentBasedPrimerMatcher]] and [[GappedAlignmentBasedPrimerMatcher]]*/
  private case class AlignmentResult
  (
    mmScore: Option[Int]     = None,
    mmNextScore: Option[Int] = None,
    fullNextBest: Int        = minAlignmentScore
  ) {
    if (mmNextScore.isDefined) require(mmScore.isDefined)
  }

  /** The primers to test. */
  private val primers = IndexedSeq(
    // this pair matches no other pair
    Primer("1", "1+", "AAAAAAAAAA", "chr1", 1,      10, true), Primer("1", "1-", "TTTTTTTTTT", "chr1", 101,   110, false),
    // the next two pairs are one mismatch apart, so will return a second best hit on full alignment, but are too close
    // for the mismatch alignment (see minMismatchDelta)
    Primer("2", "2+", "TTTTTTTTTT", "chr2", 1,      10, true), Primer("2", "2-", "AAAAAAAAAA", "chr2", 101,   110, false),
    Primer("3", "3+", "TTTTTATTTT", "chr2", 1001, 1010, true), Primer("3", "3-", "AAAAATAAAA", "chr2", 1001, 1010, false),
    // the next two pairs are three mismatch apart, so will not return a second best hit on full alignment (see the
    // scoring parameters make it fall below the minAlignmentScore), but are fare enough apart for the mismatch
    // alignment (see minMismatchDelta)
    Primer("4", "4+", "GGGGGGGGGG", "chr3", 2001, 2010, true), Primer("4", "4-", "CCCCCCCCCC", "chr3", 2001, 2010, false),
    Primer("5", "5+", "GGGGAAAGGG", "chr3", 3001, 3010, true), Primer("5", "5-", "CCCCTTTCCC", "chr3", 3001, 3010, false),
    // this pair matches no other pair
    Primer("6", "6+", "GATTACAGAT", "chr4", 1001, 1010, true), Primer("6", "6-", "GATTACAGAT", "chr4", 1001, 1010, false),
  )

  /** The alignment results for each primer, in the same order. */
  private val results = IndexedSeq(
    AlignmentResult(0),                AlignmentResult(0),
    AlignmentResult(fullNextBest = 6), AlignmentResult(fullNextBest = 6),
    AlignmentResult(fullNextBest = 6), AlignmentResult(fullNextBest = 6),
    AlignmentResult(0, 3),             AlignmentResult(0, 3),
    AlignmentResult(0, 3),             AlignmentResult(0, 3),
    AlignmentResult(0),                AlignmentResult(0),
  )
  require(results.length == primers.length)

  private def buildBuilder(primers: Seq[Primer]): SamBuilder = {
    val builder = new SamBuilder(readLength=10)
    // location match: use all the primers!
    primers.foreach { primer =>
      fromPrimer(builder, primer).foreach { r =>
        r("lc") = 1
        r("id") = primer.primer_id
      }
    }

    // mismatch match: change reference and add mismatches to the bases
    Range.inclusive(0, maxMismatches+1).foreach { numMismatches =>
      primers.foreach { primer =>
        fromPrimer(builder, primer, Some(Mismatch(numMismatches)))
      }
    }

    // full alignment: change reference and delete a base
    Range.inclusive(1, 2).foreach { indelLength =>
      primers.foreach { primer =>
        fromPrimer(builder, primer, Some(Deletion(indelLength)))
      }
      primers.foreach { primer =>
        fromPrimer(builder, primer, Some(Insertion(indelLength)))
      }
    }

    // no matches
    builder.addFrag(bases = "G"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Plus)
    builder.addFrag(bases = "C"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Minus)

    // unmapped and no matches
    builder.addFrag(bases = "G"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Plus, unmapped = true)
    builder.addFrag(bases = "C"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Minus, unmapped = true)
    builder
  }

  private def pairRecords(builder: SamBuilder): Seq[SamRecord] = {
    val pairs = builder.iterator.grouped(2)
      .map { case Seq(r1, r2) => Seq(r1.clone(), r2.clone()) }
      .flatMap { case Seq(r1, r2) =>
        r2.name         = r1.name
        r1.paired       = true
        r2.paired       = true
        r1.firstOfPair  = true
        r2.secondOfPair = true
        SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
        Seq(r1, r2)
      }.toList
    require(pairs.length == builder.toSeq.length, s"${pairs.length} == ${builder.toSeq.length}")
    pairs
  }

  /** Create paired end reads from the given primers. */
  private val (pairs: Seq[SamRecord], header: SAMFileHeader) = {
    // Paired end reads are created as fragments, with pairing information added later.
    val fragBuilder: SamBuilder = {
      val builder = buildBuilder(primers)

      // non-canonical when made into a pairs
      {
        fromPrimer(builder, primers(0))
        fromPrimer(builder, primers(3))
      }

      builder
    }

    (pairRecords(fragBuilder), fragBuilder.header)
  }

  // The type of edit to perform
  private sealed trait EditType
  private case class Mismatch(numMismatches: Int) extends EditType
  private case class Insertion(length: Int) extends EditType
  private case class Deletion(length: Int) extends EditType

  /** Adds `numNs` Ns into the middle of the read. */
  private def addNsInTheMiddle(seq: String, numNs: Int, fwd: Boolean): String = {
    val str = seq.length - numNs match {
      case n if n <= 0 => seq
      case 1 if fwd    => ("N" :+ seq.drop(1)).mkString
      case 1 if !fwd   => (seq.drop(1) +: "N").mkString
      case n           => seq.take((n+1)/2) ++ ("N" * numNs) ++ seq.takeRight(n/2)
    }
    require(str.count(_ == 'N') == numNs, s"n=$numNs $seq => $str")
    require(str.length == seq.length, s"n=$numNs $seq => $str")
    str
  }

  /** Deletes the given number of bases in the middle of the read. */
  private def deleteInTheMiddle(seq: String, indelLength: Int): String = {
    val remaining = seq.length - indelLength
    require(remaining > 0)
    val str = seq.take((remaining+1)/2) + seq.takeRight(remaining/2)
    require(str.length == remaining)
    str
  }

  /** Inserts Ns into the middle of the read. */
  private def insertInTheMiddle(seq: String, indelLength: Int): String = {
    val str = seq.take((seq.length+1)/2) + ("N"*indelLength) + seq.takeRight(seq.length/2)
    require(str.length == seq.length + indelLength)
    str
  }

  /** Creates a [[SamRecord]] for the given primer, modifying it based on the edit type. */
  private def fromPrimer(builder: SamBuilder,
                         primer: Primer,
                         edit: Option[EditType] = None): Option[SamRecord] = {

    val dict       = builder.dict
    val newRefName = dict.getSequence(dict.getSequences.size() - 1).getSequenceName
    val newPrimer: Primer  = edit match {
      case None                 => primer
      case Some(mm: Mismatch)   => primer.copy(ref_name=newRefName, sequence = addNsInTheMiddle(primer.sequence, mm.numMismatches, primer.positiveStrand))
      case Some(ins: Insertion) => primer.copy(ref_name=newRefName, sequence = insertInTheMiddle(primer.sequence, ins.length), end=primer.end+ins.length)
      case Some(del: Deletion)  => primer.copy(ref_name=newRefName, sequence = deleteInTheMiddle(primer.sequence, del.length), start=primer.start+del.length)
      case _                    => unreachable(s"Unknown edit: $edit")
    }

    val refIndex = dict.getSequenceIndex(newPrimer.ref_name)
    val strand   = if (newPrimer.positiveStrand) Plus else Minus
    val bases    = new String(newPrimer.basesInGenomicOrder)
    val quals    = (33 + builder.baseQuality).toChar.toString * newPrimer.length
    builder.addFrag(bases = bases, quals = quals, contig = refIndex, start = newPrimer.start, cigar = s"${newPrimer.length}M", strand = strand).map { rec =>
      val editTagValue = edit match {
        case None                 => "no_edits"
        case Some(mm: Mismatch)   => s"mis_${mm.numMismatches}"
        case Some(ins: Insertion) => s"ins_${ins.length}"
        case Some(del: Deletion)  => s"del_${del.length}"
        case _                    => unreachable(s"Unknown edit: $edit")
      }
      rec("et") = editTagValue
      rec
    }
  }

  "IdentifyPrimers" should s"run end to end on a two primer pairs" in {
    val inputPrimers = Seq(
      Primer("1", "1+", "AGCTAGCTAG", "chr1", 1,      10, true),
      Primer("1", "1-", "CTAGCTAGCT", "chr1", 101,   110, false),
      Primer("2", "2+", "ACTGACTGAC", "chr2", 1,      10, true),
      Primer("2", "2-", "GTCAGTCAGT", "chr3", 101,   110, false)
    )
    val (inputRecords, header) = {
      val builder = buildBuilder(inputPrimers)
      (pairRecords(builder), builder.header)
    }

    val input = {
      val path = makeTempFile("input.", ".bam")
      val writer = SamWriter(path, header)
      writer ++= inputRecords
      writer.close()
      path
    }

    val metrics = makeTempFile("metrics.", ".prefix")
    val output  = makeTempFile("output.", ".bam")
    val primers = makeTempFile("primers.", ".tab")

    Primer.write(primers, inputPrimers)

    val tool = new IdentifyPrimers(
      input                 = input,
      primerPairs           = primers,
      metrics               = metrics,
      output                = output,
      slop                  = slop,
      maxMismatchRate       = maxMismatcheRate,
      minAlignmentScoreRate = 0.1, // so we get at most a 1bp indel
      matchScore            = matchScore,
      mismatchScore         = mismatchScore,
      gapOpen               = gapOpen,
      gapExtend             = gapExtend
    )
    executeFgbioTool(tool)

    val outputRecords = readBamRecs(output).toIndexedSeq
    outputRecords.length shouldBe inputRecords.length

    outputRecords.foreach { rec =>
      val matchInfo = if (rec.positiveStrand) rec[String]("f5") else rec[String]("r5")
      val matchType = if (matchInfo.contains(",")) matchInfo.split(',')(7) else matchInfo
      rec.get[String]("et") match {
        case None     => matchType shouldBe "none"
        case Some(et) =>
          et.split('_').toSeq match {
            case Seq("no", "edits") => matchType shouldBe PrimerMatch.LocationName
            case Seq("mis", n)      => matchType shouldBe (if (n.toInt < 3) PrimerMatch.UngappedName else "none")
            case Seq("ins", n)      => matchType shouldBe (if (n.toInt < 2) PrimerMatch.GappedName else "none")
            case Seq("del", n)      => matchType shouldBe (if (n.toInt < 2) PrimerMatch.GappedName else "none")
          }
      }
    }
  }

  Seq(true, false).foreach { primersAreMapped =>
    val primerMessage = if (primersAreMapped) "mapped primers" else "unmapped primers"
    Seq(true, false).foreach { readsAreMapped =>
      val readsMessage = if (readsAreMapped) "mapped reads" else "unmapped reads"
      it should s"run end to end $primerMessage and $readsMessage" in {
        val input   = {
          val path = makeTempFile("input.", ".bam")
          val writer = SamWriter(path, header)
          if (readsAreMapped) {
            writer ++= this.pairs
          }
          else {
            writer ++= this.pairs.map { rec =>
              val r = rec.clone()
              SAMUtils.makeReadUnmapped(r.asSam)
              r
            }
          }
          writer.close()
          path
        }
        val metrics = makeTempFile("metrics.", ".prefix")
        val output  = makeTempFile("output.", ".bam")
        val primers = makeTempFile("primers.", ".tab")

        if (primersAreMapped) {
          Primer.write(primers, this.primers)
        }
        else {
          Primer.write(primers, this.primers.map(_.copy(ref_name="")))
        }

        val tool = new IdentifyPrimers(
          input                 = input,
          primerPairs           = primers,
          metrics               = metrics,
          output                = output,
          slop                  = slop,
          maxMismatchRate       = maxMismatcheRate,
          minAlignmentScoreRate = minAlignmentScoreRate,
          matchScore            = matchScore,
          mismatchScore         = mismatchScore,
          gapOpen               = gapOpen,
          gapExtend             = gapExtend
        )

        executeFgbioTool(tool)

        val actual = Metric.read[IdentifyPrimersMetric](PathUtil.pathTo(metrics + ".summary.txt")) match {
          case Seq(s)    => s
          case summaries => fail(s"Found ${summaries.length} summary metrics, should be only one.")
        }

        // relationships across metric groups
        (actual.templates * 2) shouldBe this.pairs.length
        (actual.pairs * 2) shouldBe this.pairs.length
        (actual.paired_matches + actual.unpaired_matches + actual.no_paired_matches) shouldBe (this.pairs.length / 2)
        actual.match_attempts  shouldBe (actual.location + actual.ungapped + actual.gapped + actual.no_match)

        // create the expected set of metrics
        val expected = {
          // Common metrics between mapped and unmapped test scenarios
          val commonMetrics = new IdentifyPrimersMetric(
            // read pair match counts
            paired_matches     = 48,
            unpaired_matches   = 0,
            no_paired_matches  = 15,
            // counts of template types
            templates          = 63,
            pairs              = 63,
            fragments          = 0,
            mapped_pairs       = if (readsAreMapped) 62 else 0,
            unpaired           = 0,
            unmapped_pairs     = if (readsAreMapped) 1 else 63,
            fr_pairs           = if (readsAreMapped) 61 else 0,
            mapped_fragments   = 0,
            unmapped_fragments = 0,
            // counts of types of individual primer matches
            match_attempts     = 126,
            no_match           = 30
          )
          if (primersAreMapped && readsAreMapped) {
            // location-based matching was performed
            commonMetrics.copy(
              // counts of types of individual primer matches
              location           = 14,
              ungapped           = 78,
              gapped             = 4
            )
          }
          else {
            // only ungapped and gapped alignment
            commonMetrics.copy(
              // counts of types of individual primer matches
              location           = 0,
              ungapped           = 92,
              gapped             = 4
            )
          }
        }

        actual.zip(expected).foreach { case (act, exp) => act shouldBe exp }  // NB: this helps show **which** metric is different
      }
    }
  }
}