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

import java.nio.file.Path

import com.fulcrumgenomics.alignment.{Aligner, Mode}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus, Strand}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.util.SequenceUtil
import org.scalatest.OptionValues
import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.bam.identifyprimers.IdentifyPrimersMetric
//import com.fulcrumgenomics.bam.IdentifyPrimers._
import com.fulcrumgenomics.bam.identifyprimers.IdentifyPrimers
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.SamPairUtil

/*
trait PrimerMatcherTestData {
  // Defaults for all tests
  protected val slop                  = 1
  protected val maxMismatcheRate      = 3/10d
  protected val minAlignmentScoreRate = 5/10d

  // The default aligner
  protected val matchScore    = 1
  protected val mismatchScore = -3
  protected val gapOpen       = -6
  protected val gapExtend     = -1
  protected val aligner       = Aligner(matchScore, mismatchScore, gapOpen, gapExtend, mode=Mode.Glocal)

  /** Companion object to [[AlignmentResult]] to help making objects with defaults. */
  protected object AlignmentResult {
    def apply(mmScore: Int): AlignmentResult = new AlignmentResult(mmScore=Some(mmScore))
    def apply(mmScore: Int, mmNextScore: Int): AlignmentResult = new AlignmentResult(mmScore=Some(mmScore), mmNextScore = Some(mmNextScore))
    def apply(mmScore: Int, mmNextScore: Int, fullNextBest: Int): AlignmentResult = new AlignmentResult(mmScore = Some(mmScore), mmNextScore = Some(mmNextScore))
  }

  /** Stores the results of running [[PrimerMatcher.matchWithMismatchAlignment()]] and [[PrimerMatcher.matchWithFullAlignment()]]*/
  protected case class AlignmentResult
  (
    mmScore: Option[Int]     = None,
    mmNextScore: Option[Int] = None,
    fullNextBest: Int        = minAlignmentScore
  ) {
    if (mmNextScore.isDefined) require(mmScore.isDefined)
  }

  /** The primers to test. */
  protected val primers = IndexedSeq(
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
  protected val results = IndexedSeq(
    AlignmentResult(0),                AlignmentResult(0),
    AlignmentResult(fullNextBest = 6), AlignmentResult(fullNextBest = 6),
    AlignmentResult(fullNextBest = 6), AlignmentResult(fullNextBest = 6),
    AlignmentResult(0, 3),             AlignmentResult(0, 3),
    AlignmentResult(0, 3),             AlignmentResult(0, 3),
    AlignmentResult(0),                AlignmentResult(0),
  )
  require(results.length == primers.length)

  protected val fragBuilder: SamBuilder = {
    val builder = new SamBuilder(readLength=10)

    // location match: use all the primers!
    this.primers.foreach { primer =>
      fromPrimer(builder, primer).foreach(r => r("lc") = 1)
    }

    // mismatch match: change reference and add mismatches to the bases
    Range.inclusive(0, maxMismatches+1).foreach { numMismatches =>
      this.primers.foreach { primer =>
        fromPrimer(builder, primer, Mismatch(numMismatches))
      }
    }

    // full alignment: change reference and delete a base
    Range.inclusive(1, 2).foreach { indelLength =>
      this.primers.foreach { primer =>
        fromPrimer(builder, primer, Deletion(indelLength))
      }
      this.primers.foreach { primer =>
        fromPrimer(builder, primer, Insertion(indelLength))
      }
    }

    // no matches
    builder.addFrag(bases = "G"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Plus)
    builder.addFrag(bases = "C"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Minus)

    // unmapped and no matches
    builder.addFrag(bases = "G"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Plus, unmapped = true)
    builder.addFrag(bases = "C"*builder.readLength, contig = 0, start = 100000, cigar = s"${builder.readLength}M", strand = Minus, unmapped = true)

    // non-canonical when converted to pairs
    {
      fromPrimer(builder, primers(0))
      fromPrimer(builder, primers(3))
    }

    builder
  }

  protected val pairs: Seq[SamRecord] = {
    // It's funny we have to write to a file and read back in to clone SamRecords!
    val path  = fragBuilder.toTempFile()
    val in    = SamSource(path)
    val pairs = in.iterator.grouped(2).flatMap { case Seq(r1, r2) =>
      r2.name         = r1.name
      r1.paired       = true
      r2.paired       = true
      r1.firstOfPair  = true
      r2.secondOfPair = true
      SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
      Seq(r1, r2)
    }.toList
    in.close()
    require(pairs.length == fragBuilder.toSeq.length, s"${pairs.length} == ${fragBuilder.toSeq.length}")
    pairs
  }


  // The type of edit to perform
  protected sealed trait EditType
  protected case class Mismatch(numMismatches: Int) extends EditType
  protected case class Insertion(length: Int) extends EditType
  protected case class Deletion(length: Int) extends EditType

  def addNsInTheMiddle(seq: String, numNs: Int, fwd: Boolean): String = {
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

  def deleteInTheMiddle(seq: String, indelLength: Int): String = {
    val remaining = seq.length - indelLength
    require(remaining > 0)
    val str = seq.take((remaining+1)/2) + seq.takeRight(remaining/2)
    require(str.length == remaining)
    str
  }

  def insertInTheMiddle(seq: String, indelLength: Int): String = {
    val str = seq.take((seq.length+1)/2) + ("N"*indelLength) + seq.take(seq.length/2)
    require(str.length == seq.length + indelLength)
    str
  }

  /** Creates a */
  protected def fromPrimer(builder: SamBuilder,
                           primer: Primer,
                           edit: Option[EditType] = None): Option[SamRecord] = {

    val dict       = builder.dict
    val newRefName = dict.getSequence(dict.getSequences.size() - 1).getSequenceName
    val newPrimer: Primer  = edit match {
      case None                 => primer
      case Some(mm: Mismatch)   => primer.copy(ref_name=newRefName, sequence = addNsInTheMiddle(primer.sequence, mm.numMismatches, primer.forward))
      case Some(ins: Insertion) => primer.copy(ref_name=newRefName, sequence = insertInTheMiddle(primer.sequence, ins.length), end=primer.end+ins.length)
      case Some(del: Deletion)  => primer.copy(ref_name=newRefName, sequence = deleteInTheMiddle(primer.sequence, del.length), start=primer.start+del.length)
      case _                    => unreachable(s"Unknown edit: $edit")
    }

    val refIndex = dict.getSequenceIndex(newPrimer.ref_name)
    val strand   = if (newPrimer.forward) Plus else Minus
    val quals    = (33 + builder.baseQuality).toChar.toString * newPrimer.length
    // NB: do not reverse complement the bases as the primer's bases is on the forward strand already!
    builder.addFrag(bases = newPrimer.sequence, quals = quals, contig = refIndex, start = newPrimer.start, cigar = s"${newPrimer.length}M", strand = strand).map { rec =>
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

  protected def fromPrimer(builder: SamBuilder,
                           primer: Primer,
                           edit: EditType): Option[SamRecord] = fromPrimer(builder, primer, Some(edit))

}

final class IdentifyPrimersTest extends UnitSpec with OptionValues with PrimerMatcherTestData {

  "IdentifyPrimers" should "run end to end" in {
    val input   = {
      val path = makeTempFile("input.", ".bam")
      val writer = SamWriter(path, fragBuilder.header)
      writer ++= this.pairs
      writer.close()
      path
    }
    val metrics = makeTempFile("metrics.", ".prefix")
    val output  = makeTempFile("output.", ".bam")
    val primers = makeTempFile("primers.", ".tab")

    Primer.write(primers, this.primers)

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

    {
      val actual = Metric.read[IdentifyPrimersMetric](PathUtil.pathTo(metrics + ".summary.txt")) match {
        case Seq(s)    => s
        case summaries => fail(s"Found ${summaries.length} summary metrics, should be only one.")
      }

      // relationships across metric groups
      actual.templates shouldBe actual.total_primer_pair_types
      actual.total_primer_pair_types   shouldBe (actual.canonical_primer_pair + actual.non_canonical_primer_pair + actual.single_primer_pair + actual.no_primer_pair)
      actual.match_attempts  shouldBe (actual.location + actual.mismatch + actual.full_alignment + actual.no_match)

      val expected = new IdentifyPrimersMetric(
        // template types
        templates      = 63,
        pairs          = 63,
        mapped_pairs   = 62,
        unmapped_pairs = 1,
        // primer pair match types
        total_primer_pair_types    = 63,
        canonical_primer_pair      = 41,
        non_canonical_primer_pair  = 1,
        single_primer_pair         = 4,
        no_primer_pair             = 17,
        // primer match types
        match_attempts = 92,
        location       = 14,
        mismatch       = 68,
        full_alignment = 6,
        no_match       = 4
      )

      actual.zip(expected).foreach { case (act, exp) => act shouldBe exp }  // NB: this helps show **which** metric is different
    }
  }

  // TODO: a few more tests, including but not limited to fragment reads.
}
*/