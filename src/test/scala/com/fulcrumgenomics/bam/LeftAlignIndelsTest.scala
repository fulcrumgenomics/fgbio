/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata, Topology}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.{BlockCompressedOutputStream, SequenceUtil}
import org.scalatest.OptionValues

import java.nio.file.{Files, Path, StandardCopyOption}
import java.util.zip.GZIPOutputStream
import scala.util.Using

class LeftAlignIndelsTest extends UnitSpec with OptionValues {

  /** Checks the CIGAR of the padded alignment, both with and without = and X operators. */
  private def checkCigar(alignment: IndelPaddedAlignment, cigarM: String, cigarX: String): Unit = {
    withClue (f"Hint: query: ${alignment.query} target: ${alignment.target} ::") {
      alignment.cigar(useEqualsAndX = false).toString shouldBe cigarM
      alignment.cigar(useEqualsAndX = true).toString shouldBe cigarX
    }
    ()
  }

  /** Builds a fragment with the given bases, CIGAR, and attributes, then left-aligns it given the reference bases,
    * checking the returned value and the resulting CIGAR.  Returns the record for further checks. */
  private def align(bases: String,
                    cigar: String,
                    ref: => String,
                    changed: Boolean,
                    expected: String,
                    useEqualsAndX: Boolean = false,
                    attrs: Map[String, Any] = Map.empty,
                    trim: Boolean = false,
                    expectedStart: Int = 100,
                    samtoolsMismatches: Boolean = false): SamRecord = {
    val builder = new SamBuilder(readLength=bases.length)
    val rec     = builder.addFrag(bases=bases, cigar=cigar, start=100, attrs=attrs).value
    withClue(f"Hint: bases: $bases cigar: $cigar ::") {
      aligner(useEqualsAndX=useEqualsAndX, trim=trim, samtoolsMismatches=samtoolsMismatches).align(rec, ref).changed shouldBe changed
      rec.cigar.toString shouldBe expected
      rec.start shouldBe expectedStart
    }
    rec
  }

  /** Pretends the record was left-aligned to the given CIGAR, start, and edit distance, storing its original alignment
    * as the tool does. */
  private def pretendAligned(rec: SamRecord, cigar: String, nm: Int, start: Int = -1): Unit = {
    IndelLeftAligner.setOriginalAlignment(rec, rec.start, rec.cigar, nm)
    rec.cigar = cigar
    if (start > 0) rec.start = start
  }

  /** Builds an aligner that does not trim leading indels unless `trim` is true. */
  private def aligner(useEqualsAndX: Boolean = false, trim: Boolean = false, samtoolsMismatches: Boolean = false): IndelLeftAligner =
    new IndelLeftAligner(useEqualsAndX=useEqualsAndX, trimLeadingIndels=trim, samtoolsMismatches=samtoolsMismatches)

  /** Removes the hard-clip, soft-clip, and gap characters from a padded string. */
  private def unpad(str: String): String = str.filter(b => SequenceUtil.isValidBase(b.toByte))

  // A complicated alignment, with hard-clips, soft-clips, insertions, and deletions
  private val ComplicatedQuery  = unpad("HHHHHAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT----------GGGGGGGGGG----------CCCCCCCCCCTTTHH")
  private val ComplicatedTarget = unpad("HHHHHXXXXTTTTTTTTTT----------TTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCXXXHH")
  private val ComplicatedCigar  = "5H4S10M10I10M10D10M10D10M3S2H"

  "IndelPaddedAlignment.cigar" should "return the cigar for a simple alignment" in {
    checkCigar(IndelPaddedAlignment("AAAAAAAAAA", "AAAAATAAAA"), "10M", "5=1X4=")
  }

  it should "return the cigar for internal deletion" in {
    checkCigar(IndelPaddedAlignment("AA--AAAAAA", "AATTTAAAAA"), "2M2D6M", "2=2D1X5=")
  }

  it should "return the cigar for an internal insertion" in {
    checkCigar(IndelPaddedAlignment("AATTAAAAAA", "AA--AAAAAG"), "2M2I6M", "2=2I5=1X")
  }

  it should "return the cigar for a combination of insertions, deletions, and matches" in {
    checkCigar(IndelPaddedAlignment("A-AAAA-A", "-AAAAAAA"), "1I1D4M1D1M", "1I1D4=1D1=")
  }

  it should "return the cigar for an alignment without matches" in {
    checkCigar(IndelPaddedAlignment("A-", "-A"), "1I1D", "1I1D")
  }

  it should "treat an X in the target as a base, not as clipping" in {
    checkCigar(IndelPaddedAlignment("AAA", "AXA"), "3M", "1=1X1=")
  }

  "IndelPaddedAlignment.apply(Cigar, String, String)" should "build a padded alignment" in {
    checkCigar(IndelPaddedAlignment(Cigar("5M"), "AAAAA", "AATAA"), "5M", "2=1X2=")
    checkCigar(IndelPaddedAlignment(Cigar("4M5I4M5D10M"), "A"*23, "A"*23), "4M5I4M5D10M", "4=5I4=5D10=")
  }

  it should "build a padded alignment for a complicated cigar" in {
    checkCigar(IndelPaddedAlignment(Cigar("10M10I10M10D10M10D10M"), ComplicatedQuery.drop(4).dropRight(3), ComplicatedTarget),
      "10M10I10M10D10M10D10M", "10=10I10=10D10=10D10=")
  }

  it should "reject reference skips, padding, and clipping" in {
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("2M2N2M"), "AAAA", "AAAAAA") }
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("2M1P2M"), "AAAA", "AAAA") }
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("1S3M"), "AAAA", "AAA") }
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("1H3M"), "AAA", "AAA") }
  }

  it should "reject the gap character in the query or target" in {
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("3M"), "AAA", "A-A") }
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment(Cigar("3M"), "A-A", "AAA") }
  }

  /** Checks the left-aligned CIGAR, and that the same CIGAR (with I and D swapped) results when swapping the query and
    * target. */
  private def testLeftAlign(query: String, target: String, cigar: String): Unit = {
    withClue (f"Hint: query: $query target: $target ::") {
      IndelPaddedAlignment(query, target).leftAligned.cigar().toString shouldBe cigar
      val swapped = cigar.map { case 'I' => 'D'; case 'D' => 'I'; case c => c }
      IndelPaddedAlignment(target, query).leftAligned.cigar().toString shouldBe swapped
    }
    ()
  }

  "IndelPaddedAlignment.leftAligned" should "not left align an alignment that has no indels" in {
    testLeftAlign("AAAAAAAAAA", "AAAAATAAAA", "10M")
  }

  it should "not left align when a deletion cannot be shifted" in {
    testLeftAlign("GGTT--TAAAAGGGTT", "GGTTAAAAAAAGGGTT", "4M2D10M")
  }

  it should "shift a deletion" in {
    // same as above, but now allows a single base to shift
    testLeftAlign("GGTT--TAAAAGGGTT", "GGTTATTAAAAGGGTT", "3M2D11M")
    // same as above, but now allows two bases to shift
    testLeftAlign("GGTT--TAAAAGGGTT", "GGTTTTTAAAAGGGTT", "2M2D12M")
  }

  it should "not modify an alignment that has a leading left-aligned indel" in {
    testLeftAlign("AAAAAAAAAA", "-AAAAAAAAA", "1I9M")
  }

  it should "left align a 1bp indel" in {
    testLeftAlign("AAAAAAAAAA", "AAAAA-AAAA", "1I9M")
  }

  it should "left align a trailing 1bp indel" in {
    testLeftAlign("AAAAAAAAAA", "AAAAAAAAA-", "1I9M")
  }

  it should "left align and join two 1bp indel" in {
    testLeftAlign("AAAAAAAAAA", "AA-AAAAAA-", "2I8M")
  }

  it should "not modify an already left aligned 4bp motif indel" in {
    testLeftAlign("TTTTCGCGCGCGTTTT", "TTTT----CGCGTTTT", "4M4I8M")
  }

  it should "left align a 4bp motif indel" in {
    testLeftAlign("TTTTTACGTACGTTTT", "TTTTTACG----TTTT", "4M4I8M")
  }

  it should "left align and join two separated 4bp motif indels" in {
    testLeftAlign("TTTTTACGTACGTACGTACGTTTT", "TTTTTACGT----ACGT----TTT", "4M8I12M")
  }

  it should "not shift an indel when doing so would change which bases are aligned" in {
    // shifting the insertion left would pair the reference A with the read A instead of the read C
    testLeftAlign("CAA", "A-A", "1M1I1M")
    testLeftAlign("GCAAT", "GA-AT", "2M1I2M")
  }

  it should "shift an indel across mismatches that are preserved by the shift" in {
    // the read base before the insertion equals the inserted base, so the C/G mismatch is kept
    testLeftAlign("TCAAAT", "TGA-AT", "2M1I3M")
    // the reference base before the deletion equals the deleted base, so the C/G mismatch is kept but shifted past
    testLeftAlign("TC-T", "TGGT", "1M1D2M")
  }

  it should "fully left-align indels that merge with a preceding indel" in {
    // a single pass would merge the two insertions (2M2I2M) but not shift the merged insertion
    testLeftAlign("TACAAT", "TA-A-T", "1M2I3M")
  }

  it should "place insertions before deletions in a run of adjacent indels" in {
    // shifting the second insertion left places it after the deletion, giving 1I1D1I without reordering
    IndelPaddedAlignment("GT-AAT", "G-CA-T").leftAligned.cigar().toString shouldBe "1M2I1D2M"
    // shifting the second deletion left places it after the insertion, giving 1D1I1D without reordering
    IndelPaddedAlignment("G-CA-T", "GT-AAT").leftAligned.cigar().toString shouldBe "1M1I2D2M"
  }

  it should "shift an insertion further after placing it before a deletion" in {
    // the insertion (C) can only be shifted past the preceding C once it is placed before the deletion (G)
    IndelPaddedAlignment("AC-CT", "ACG-T").leftAligned.cigar().toString shouldBe "1M1I1M1D1M"
  }

  it should "give the same result for adjacent insertions and deletions in either order" in {
    // the deletion (A) can be shifted past the preceding A even when it follows the insertion (C)
    IndelPaddedAlignment("GA-CT", "GAA-T").leftAligned.cigar().toString shouldBe "1M1D1M1I1M"
    IndelPaddedAlignment("GAC-T", "GA-AT").leftAligned.cigar().toString shouldBe "1M1D1M1I1M"
    // both the insertion and the deletion can be shifted: the insertion is shifted first
    IndelPaddedAlignment("TA-AT", "TAA-T").leftAligned shouldBe IndelPaddedAlignment("TAA-T", "TA-AT").leftAligned
    IndelPaddedAlignment("TA-AT", "TAA-T").leftAligned.cigar().toString shouldBe "1M1I1D2M"
  }

  it should "reject a column that is a gap in both the query and target" in {
    an[IllegalArgumentException] shouldBe thrownBy { IndelPaddedAlignment("A-A", "A-A") }
  }

  it should "not shift an indel across bases that differ, even ambiguous bases that htsjdk compares as equal" in {
    // htsjdk's SequenceUtil.basesEqual treats X and Z as the same base, which would pair the read A with Z instead of X
    testLeftAlign("A-A", "XZA", "1M1D1M")
    // lower-case repeats are still shifted across
    testLeftAlign("gaa-t", "gaaat", "1M1D3M")
  }

  it should "left-align long isolated pure indel runs without changing any converged columns" in {
    for (deleted <- Seq(false, true); size <- Seq(500, 25000)) {
      val bases = "A" * (size * 2)
      val gaps = "A" * size + "-" * size
      val input = if (deleted) IndelPaddedAlignment(gaps, bases) else IndelPaddedAlignment(bases, gaps)
      val expected = if (deleted) IndelPaddedAlignment("-" * size + "A" * size, bases)
                     else IndelPaddedAlignment(bases, "-" * size + "A" * size)
      input.leftAligned shouldBe expected
    }
  }

  it should "give the same result as shifting every run of indels in each pass" in {
    // The left-alignment before it visited only the runs changed in the previous pass
    def reference(alignment: IndelPaddedAlignment): IndelPaddedAlignment = {
      import IndelPaddedAlignment.{Indel, sameBase}
      val queryArr  = alignment.query.toCharArray
      val targetArr = alignment.target.toCharArray
      val length    = alignment.length
      def isIndel(index: Int): Boolean = queryArr(index) == Indel || targetArr(index) == Indel
      def runEnd(start: Int): Int = { var end = start; while (end < length && isIndel(end)) end += 1; end }
      def reorderAll(): Unit = {
        var index = 0
        while (index < length) {
          if (!isIndel(index)) index += 1 else {
            val end      = runEnd(index)
            val inserted = (index until end).filter(i => targetArr(i) == Indel).map(queryArr)
            val deleted  = (index until end).filter(i => queryArr(i) == Indel).map(targetArr)
            inserted.indices.foreach { i => queryArr(index + i) = inserted(i); targetArr(index + i) = Indel }
            deleted.indices.foreach  { i => queryArr(index + inserted.length + i) = Indel; targetArr(index + inserted.length + i) = deleted(i) }
            index = end
          }
        }
      }
      def shiftAll(): Boolean = {
        var shifted = false
        var start   = 0
        while (start < length) {
          if (!isIndel(start)) start += 1 else {
            val end         = runEnd(start)
            val numInserted = (start until end).count(i => targetArr(i) == Indel)
            val numDeleted  = end - start - numInserted
            val prev        = start - 1
            val inserted    = queryArr.slice(start, start + numInserted)
            val deleted     = targetArr.slice(start + numInserted, end)
            if (prev >= 0 && numInserted > 0 && sameBase(queryArr(prev), inserted.last)) {
              val (prevQuery, prevTarget) = (queryArr(prev), targetArr(prev))
              (0 until numInserted).foreach { i => queryArr(prev + i) = if (i == 0) prevQuery else inserted(i - 1); targetArr(prev + i) = Indel }
              queryArr(prev + numInserted)  = inserted.last
              targetArr(prev + numInserted) = prevTarget
              shifted = true
            }
            else if (prev >= 0 && numDeleted > 0 && sameBase(targetArr(prev), deleted.last)) {
              val (prevQuery, prevTarget) = (queryArr(prev), targetArr(prev))
              (0 until numDeleted).foreach { i => queryArr(prev + i) = Indel; targetArr(prev + i) = if (i == 0) prevTarget else deleted(i - 1) }
              queryArr(prev + numDeleted)  = prevQuery
              targetArr(prev + numDeleted) = deleted.last
              (0 until numInserted).foreach { i => queryArr(prev + numDeleted + 1 + i) = inserted(i); targetArr(prev + numDeleted + 1 + i) = Indel }
              shifted = true
            }
            start = end
          }
        }
        shifted
      }
      var shifted = true
      while (shifted) { reorderAll(); shifted = shiftAll() }
      IndelPaddedAlignment(new String(queryArr), new String(targetArr))
    }

    // Repeats, mixed runs that merge, and case-insensitive shifts must preserve the exact converged columns.
    val random = new scala.util.Random(42)
    Seq("AC", "ACac").foreach { alphabet =>
      forloop(from=0, until=20000) { _ =>
        val columns = IndexedSeq.fill(1 + random.nextInt(30)) {
          def base: Char = alphabet(random.nextInt(alphabet.length))
          random.nextInt(4) match {
            case 0 => (base, '-')
            case 1 => ('-', base)
            case _ => (base, base)
          }
        }
        val alignment = IndelPaddedAlignment(columns.map(_._1).mkString, columns.map(_._2).mkString)
        withClue(s"query: ${alignment.query} target: ${alignment.target} ::") {
          alignment.leftAligned shouldBe reference(alignment)
        }
        ()
      }
    }
  }

  "IndelLeftAligner.align" should "convert M to = and X when useEqualsAndX is true" in {
    align("A"*10, "10M", "A"*10, changed=true, expected="10=", useEqualsAndX=true)
  }

  Seq("1M1I1D1I2M", "1M1I1I1D2M").foreach { cigar =>
    it should s"canonicalize adjacent indels in Keep mode from $cigar" in {
      val rec = align("TCCAA", cigar, "TGAA", changed=true, expected="1M2I1D2M",
        attrs=Map("MD" -> "1^G2", "NM" -> 3, "de" -> 0.5f))
      rec[String]("MD") shouldBe "1^G2"
      rec[Int]("NM") shouldBe 3
      rec.contains("de") shouldBe false
      aligner().align(rec, "TGAA").changed shouldBe false
    }
  }

  it should "return false when the alignment is unchanged" in {
    align("A"*10, "10M", "A"*10, changed=false, expected="10M")
  }

  it should "return true when an indel was shifted" in {
    align("CAAG", "3M1D1M", "CAAAG", changed=true, expected="1M1D3M")
  }

  it should "keep a leading insertion in place when not trimming" in {
    align("A"*10, "4I6M", "A"*6, changed=false, expected="4I6M")
    align("A"*10, "4I6M", "A"*6, changed=true, expected="4I6=", useEqualsAndX=true)
  }

  it should "keep a leading deletion in place when not trimming" in {
    align("A"*10, "4D10M", "A"*14, changed=false, expected="4D10M")
    align("A"*10, "4D10M", "A"*14, changed=true, expected="4D10=", useEqualsAndX=true)
  }

  it should "trim a leading insertion by soft-clipping it" in {
    align("A"*10, "6M4I", "A"*6, changed=true, expected="4S6M", trim=true)
    align("A"*10, "4I6M", "A"*6, changed=true, expected="4S6M", trim=true)
    align("A"*12, "2H1S1M1I6M1I2S", "A"*7, changed=true, expected="2H3S7M2S", trim=true)
    align("A"*10, "2H6M4I", "A"*6, changed=true, expected="2H4S6M", trim=true)
  }

  it should "trim a leading deletion by moving the start of the alignment" in {
    align("AAAA", "2M1D2M", "AAAAA", changed=true, expected="4M", trim=true, expectedStart=101)
    align("A"*10, "4D10M", "A"*14, changed=true, expected="10M", trim=true, expectedStart=104)
    align("A"*10, "2S4D8M", "A"*12, changed=true, expected="2S8M", trim=true, expectedStart=104)
  }

  it should "trim leading insertions and deletions together, and re-compute the MD, NM, and UQ tags" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="TACG", quals="IIII", cigar="1D1I3M", start=100,
      attrs=Map("MD" -> "0^C3", "NM" -> 2, "UQ" -> 0)).value
    aligner(trim=true).align(rec, "CACG") shouldBe LeftAlignResult(changed=true, trimmed=true, nm=0)
    rec.cigar.toString shouldBe "1S3M"
    rec.start shouldBe 101
    rec[String]("MD") shouldBe "3"
    rec[Int]("NM") shouldBe 0
    rec[Int]("UQ") shouldBe 0
  }

  it should "not trim a leading indel when no aligned bases follow it in the first block" in {
    align("AAAA", "2I10N2M", "AA", changed=false, expected="2I10N2M", trim=true)
  }

  it should "count a record as trimmed only when an indel was trimmed" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CAAG", cigar="3M1D1M", start=100).value
    aligner(trim=true).align(rec, "CAAAG") shouldBe LeftAlignResult(changed=true, trimmed=false, nm=1)
  }

  it should "left-align an insertion to the start" in {
    align("A"*10, "6M4I", "A"*6, changed=true, expected="4I6=", useEqualsAndX=true)
  }

  it should "left-align a deletion to the start" in {
    align("A"*10, "10M4D", "A"*14, changed=true, expected="4D10=", useEqualsAndX=true)
  }

  it should "not shift indels into leading or trailing soft-clips" in {
    align("A"*12, "1S1M1I6M1I2S", "A"*7, changed=true, expected="1S2I7M2S")
  }

  it should "treat an X in the reference as a base, not as clipping" in {
    align("AAA", "3M", "AXA", changed=true, expected="1=1X1=", useEqualsAndX=true)
  }

  it should "preserve hard-clips" in {
    align("A"*10, "5H10M4H", "A"*10, changed=true, expected="5H10=4H", useEqualsAndX=true)
  }

  it should "left-align a complicated case" in {
    align(ComplicatedQuery, ComplicatedCigar, ComplicatedTarget, changed=true, expected="5H4S10I20=20D20=3S2H", useEqualsAndX=true)
  }

  it should "left-align indels within each block separated by reference skips" in {
    val rec = align("TTTTT" + "CAAG", "5M100N3M1D1M", "TTTTTCAAAG", changed=true, expected="5M100N1M1D3M",
      attrs=Map("MD" -> "8^A1", "NM" -> 1))
    rec[String]("MD") shouldBe "6^A3"
    rec[Int]("NM") shouldBe 1
  }

  it should "not shift an indel across a reference skip" in {
    align("AAAA" + "AAAA", "4M10N1M1I2M", "AAAAAAA", changed=true, expected="4M10N1I3M")
  }

  it should "keep using = and X operators when the input CIGAR uses them" in {
    align("CAAG", "3=1D1=", "CAAAG", changed=true, expected="1=1D3=")
  }

  it should "update the MD and NM tags when the alignment changes" in {
    val rec = align("AAAA", "2M1D2M", "AAAAA", changed=true, expected="1D4M", attrs=Map("MD" -> "2^A2", "NM" -> 1))
    rec[String]("MD") shouldBe "0^A4"
    rec[Int]("NM") shouldBe 1
  }

  it should "update the MD and NM tags with mismatches, insertions, and deletions" in {
    // read: TG C C T T AAG - A (two soft-clipped bases, a C/A mismatch, two 1bp insertions, and a 1bp deletion)
    // ref:  XX A - T - AAG G A
    val rec = align("TGCCTTAAGA", "2S1M1I1M1I3M1D1M", "ATAAGGA", changed=true, expected="2S1I1M1I3M1D2M",
      attrs=Map("MD" -> "0A4^G1", "NM" -> 4))
    val expected = rec.asSam.deepCopy()
    SequenceUtil.calculateMdAndNmTags(expected, ("N"*99 + "ATAAGGA").getBytes, true, true)
    rec[String]("MD") shouldBe expected.getStringAttribute("MD")
    rec[Int]("NM") shouldBe expected.getIntegerAttribute("NM").intValue()
  }

  it should "use the same base comparison for = and X as for MD and NM" in {
    // the read N is a mismatch against the reference A in the CIGAR, MD, and NM
    val rec = align("CANG", "3M1D1M", "CAAAG", changed=true, expected="1=1D1=1X1=", useEqualsAndX=true,
      attrs=Map("MD" -> "2A0^A1", "NM" -> 2))
    rec[String]("MD") shouldBe "1^A1A1"
    rec[Int]("NM") shouldBe 2
  }

  it should "not add the MD or NM tags when they are absent" in {
    val rec = align("CAAG", "3M1D1M", "CAAAG", changed=true, expected="1M1D3M")
    rec.get[String]("MD") shouldBe None
    rec.get[Int]("NM") shouldBe None
  }

  it should "count an N in both the read and reference as a match in the CIGAR, MD, and NM, as HTSJDK does" in {
    val rec = align("NAC", "2M1D1M", "NAAC", changed=true, expected="1=1D2=", useEqualsAndX=true,
      attrs=Map("MD" -> "1A^A1", "NM" -> 1))
    rec[String]("MD") shouldBe "1^A2"
    rec[Int]("NM") shouldBe 1
    // the same as HTSJDK, as used by other fgbio tools
    val expected = rec.asSam.deepCopy()
    SequenceUtil.calculateMdAndNmTags(expected, ("N" * 99 + "NAAC").getBytes, true, true)
    rec[String]("MD") shouldBe expected.getStringAttribute("MD")
    rec[Int]("NM") shouldBe expected.getIntegerAttribute("NM").intValue
  }

  it should "count an N in both the read and reference as a mismatch in the CIGAR, MD, and NM with samtoolsMismatches" in {
    val rec = align("NAC", "2M1D1M", "NAAC", changed=true, expected="1X1D2=", useEqualsAndX=true,
      attrs=Map("MD" -> "0N1^A1", "NM" -> 2), samtoolsMismatches=true)
    rec[String]("MD") shouldBe "0N0^A2"
    rec[Int]("NM") shouldBe 2
  }

  Seq(false, true).foreach { samtoolsMismatches =>
    it should s"count the same IUPAC code in the read and reference as a match with samtoolsMismatches=$samtoolsMismatches" in {
      val rec = align("RAC", "2M1D1M", "RAAC", changed=true, expected="1=1D2=", useEqualsAndX=true,
        attrs=Map("MD" -> "1A^A1", "NM" -> 1), samtoolsMismatches=samtoolsMismatches)
      rec[String]("MD") shouldBe "1^A2"
      rec[Int]("NM") shouldBe 1
    }
  }

  it should "re-compute the NM tag when the MD tag is absent" in {
    val rec = align("AAAA", "2M1D2M", "AAAAA", changed=true, expected="4M", trim=true, expectedStart=101, attrs=Map("NM" -> 1))
    rec[Int]("NM") shouldBe 0
    rec.get[String]("MD") shouldBe None
  }

  it should "keep a CIGAR that uses = and X operators when no indel is moved" in {
    // the labels are not re-derived, even though the first base is a match
    align("CTAG", "1X1I2=", "CAG", changed=false, expected="1X1I2=")
    align("CTAG", "1X1I2=", "CAG", changed=false, expected="1X1I2=", useEqualsAndX=true)
  }

  it should "remove the UQ tag when the read has no base qualities" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CCAT", cigar="1M1I2M", start=100, attrs=Map("UQ" -> 10)).value
    rec.quals   = Array.emptyByteArray
    aligner().align(rec, "GAT").changed shouldBe true
    rec.get[Int]("UQ") shouldBe None
  }

  it should "re-compute the UQ tag when a shift changes which read base is mismatched" in {
    // the C/G mismatch moves from the first C (Q10) to the second C (Q40)
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CCAT", quals="+I??", cigar="1M1I2M", start=100, attrs=Map("UQ" -> 10)).value
    aligner().align(rec, "GAT").changed shouldBe true
    rec.cigar.toString shouldBe "1I3M"
    rec[Int]("UQ") shouldBe 40
  }

  Seq("BQ", "cs", "ds", "de").foreach { tag =>
    it should s"invalidate $tag on a changed path but retain it on notation-only changes" in {
      val value: Any = if (tag == "de") 0.3333f else "@@@@"
      val attrs = Map(tag -> value)
      val shifted = align("CCAT", "1M1I2M", "GAT", changed=true, expected="1I3M", attrs=attrs)
      shifted.contains(tag) shouldBe false
      val trimmed = align("AAAA", "2M1D2M", "AAAAA", changed=true, expected="4M", trim=true,
        expectedStart=101, attrs=attrs)
      trimmed.contains(tag) shouldBe false
      val notation = align("ACAT", "4M", "AGAT", changed=true, expected="1=1X2=", useEqualsAndX=true, attrs=attrs)
      notation[Any](tag) shouldBe value
      val unchanged = align("ACAT", "4M", "AGAT", changed=false, expected="4M", attrs=attrs)
      unchanged[Any](tag) shouldBe value
    }
  }

  it should "invalidate de when insertion events merge without changing NM" in {
    val rec = align("AAAAAA", "1M1I1M1I2M", "AAAA", changed=true, expected="2I4M",
      attrs=Map("de" -> 0.3333f, "NM" -> 2))
    rec.contains("de") shouldBe false
    rec[Int]("NM") shouldBe 2
  }

  it should "reject an unrepresentable resolved read base before changing the record" in {
    val rec = new SamBuilder(readLength=3).addFrag(name="uracil", bases="=AA", cigar="3M", start=100,
      attrs=Map("MD" -> "3", "NM" -> 0)).value
    val ex = intercept[InvalidRecordException] { aligner(useEqualsAndX=true).align(rec, "UAA") }
    ex.getMessage should include ("uracil")
    rec.basesString shouldBe "=AA"
    rec.cigar.toString shouldBe "3M"
    rec[String]("MD") shouldBe "3"
  }

  for (tag <- Seq("MM", "Mm"); negative <- Seq(false, true)) {
    it should s"reject sequence materialization before mutation with $tag and negative=$negative" in {
      val rec = new SamBuilder(readLength=3).addFrag(name="modified", bases="A=A", cigar="3M", start=100,
        attrs=Map(tag -> "A+a,1;", "ML" -> Array(200.toByte), "MN" -> 3, "MD" -> "3", "NM" -> 0)).value
      rec.negativeStrand = negative
      val quals = rec.quals.clone()
      val ex = intercept[InvalidRecordException] { aligner(useEqualsAndX=true).align(rec, "AAA") }
      ex.getMessage should include ("modified")
      ex.getMessage should include ("base-modification")
      rec.start shouldBe 100
      rec.cigar.toString shouldBe "3M"
      rec.basesString shouldBe "A=A"
      rec.quals shouldBe quals
      rec[String](tag) shouldBe "A+a,1;"
      rec[Array[Byte]]("ML") shouldBe Array(200.toByte)
      rec[Int]("MN") shouldBe 3
      rec[String]("MD") shouldBe "3"
      rec[Int]("NM") shouldBe 0
    }
  }

  Seq("MM", "Mm").foreach { tag =>
    it should s"preserve $tag and probabilities when the sequence is unchanged" in {
      val attrs: Map[String, Any] = Map(tag -> "A+a,1;", "ML" -> Array(200.toByte), "MN" -> 4)
      val rec = align("AAAA", "2M1D2M", "AAAAA", changed=true, expected="4M", trim=true,
        expectedStart=101, attrs=attrs)
      rec.basesString shouldBe "AAAA"
      rec[String](tag) shouldBe "A+a,1;"
      rec[Array[Byte]]("ML") shouldBe Array(200.toByte)
      rec[Int]("MN") shouldBe 4
      val marker = align("CA=G", "1M1I2M", "CAG", changed=false, expected="1M1I2M", attrs=attrs.updated(tag, "A+a,0;"))
      marker.basesString shouldBe "CA=G"
      marker[String](tag) shouldBe "A+a,0;"
    }
  }

  it should "resolve '=' read bases against the reference" in {
    val rec = align("A=A", "3M", "ACA", changed=true, expected="3=", useEqualsAndX=true, attrs=Map("MD" -> "3", "NM" -> 0))
    rec.basesString shouldBe "ACA"
    rec[String]("MD") shouldBe "3"
    rec[Int]("NM") shouldBe 0
  }

  it should "write the resolved base for an '=' read base that is shifted into an insertion" in {
    val rec = align("CA=AG", "3M1I1M", "CAAG", changed=true, expected="1M1I3M", attrs=Map("MD" -> "4"))
    rec.basesString shouldBe "CAAAG"
    rec[String]("MD") shouldBe "4"
  }

  it should "not change the bases of a record with '=' read bases when the CIGAR is unchanged" in {
    val rec = align("CA=G", "1M1I2M", "CAG", changed=false, expected="1M1I2M")
    rec.basesString shouldBe "CA=G"
  }

  it should "fully left-align indels that merge with a preceding indel, and be idempotent" in {
    val rec = align("TACAAT", "2M1I1M1I1M", "TAAT", changed=true, expected="1M2I3M")
    aligner().align(rec, "TAAT").changed shouldBe false
  }

  it should "not change an unmapped record" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CAAG", unmapped=true).value
    aligner(useEqualsAndX=true).align(rec, fail("reference bases should not be retrieved")).changed shouldBe false
  }

  it should "not change a record without bases (SEQ of '*')" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CAAG", cigar="3M1D1M", start=100).value
    rec.bases   = Array.emptyByteArray
    rec.quals   = Array.emptyByteArray
    aligner(useEqualsAndX=true).align(rec, fail("reference bases should not be retrieved")).changed shouldBe false
    rec.cigar.toString shouldBe "3M1D1M"
  }

  it should "not change a record with padding (P) operators" in {
    align("CAAG", "3M1P1D1M", fail("reference bases should not be retrieved"), changed=false, expected="3M1P1D1M")
  }

  it should "not retrieve the reference bases for a record that cannot change" in {
    align("A"*10, "2S8M", fail("reference bases should not be retrieved"), changed=false, expected="2S8M")
    align("CAAG", "3=1X", fail("reference bases should not be retrieved"), changed=false, expected="3=1X")
  }

  it should "fail when the reference bases do not span the alignment" in {
    intercept[InvalidRecordException] { align("CAAG", "3M1D1M", "CAAA", changed=true, expected="") }
      .getMessage should include ("Reference length mismatch")
  }

  it should "fail when the reference bases contain the gap character" in {
    intercept[InvalidRecordException] { align("CAAG", "3M1D1M", "CA-AG", changed=true, expected="") }
      .getMessage should include ("contain the gap character")
  }

  "IndelLeftAligner.fixMateAndSupplementaryInfo" should "update the mate cigar on secondary records and the cigars in SA tags" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(name="q1", start1=100, start2=300, cigar1="10M", cigar2="4M2I4M")
    val supp = builder.addFrag(name="q1", start=1000, cigar="5S5M").value
    supp.paired = true; supp.firstOfPair = true; supp.supplementary = true
    val sec  = builder.addFrag(name="q1", start=2000, cigar="10M").value
    sec.paired = true; sec.firstOfPair = true; sec.secondary = true
    sec.mateRefIndex = r2.refIndex; sec.mateStart = r2.start; sec.mateNegativeStrand = r2.negativeStrand
    sec("MC") = "4M2I4M"
    r1("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    supp("SA") = s"${r1.refName},100,+,10M,60,0;"
    pretendAligned(r2, "2I8M", nm=2)
    pretendAligned(supp, "5S5=", nm=0)

    val template = Template(Iterator(r1, r2, supp, sec))
    IndelLeftAligner.fixMateAndSupplementaryInfo(template)
    r1[String]("MC") shouldBe "2I8M"
    sec[String]("MC") shouldBe "2I8M"
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S5=,60,0;"
    supp[String]("SA") shouldBe s"${r1.refName},100,+,10M,60,0;"
  }

  it should "update SA entries for alignments that share a position and strand" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="5M5S").value
    val supp    = builder.addFrag(name="q1", start=100, cigar="5S2M1I2M").value
    supp.supplementary = true
    r1("SA")   = s"${supp.refName},100,+,5S2M1I2M,60,1;"
    supp("SA") = s"${r1.refName},100,+,5M5S,60,0;"
    pretendAligned(supp, "5S1I4M", nm=1)

    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp)))
    r1[String]("SA") shouldBe s"${supp.refName},100,+,5S1I4M,60,1;"
    supp[String]("SA") shouldBe s"${r1.refName},100,+,5M5S,60,0;"
  }

  it should "keep the clipping of an SA entry for a hard-clipped alignment" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="5M5S").value
    val supp    = builder.addFrag(name="q1", start=1000, bases="AAAAA", quals="IIIII", cigar="5H2M1I2M").value
    supp.supplementary = true
    r1("SA") = s"${supp.refName},1000,+,5S2M1I2M,60,1;"
    pretendAligned(supp, "5H1I4M", nm=1)

    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp)))
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S1I4M,60,1;"
  }

  it should "not update an SA entry that does not match exactly one alignment" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="10M").value
    val supp1   = builder.addFrag(name="q1", start=1000, cigar="5S1I4M").value
    val supp2   = builder.addFrag(name="q1", start=1000, cigar="5S1I4M").value
    Seq(supp1, supp2).foreach { s => s.supplementary = true; IndelLeftAligner.setOriginalAlignment(s, 1000, Cigar("5S2M1I2M"), 1) }
    // two alignments match the first entry, and no alignment matches the CIGAR of the second entry
    r1("SA") = s"${supp1.refName},1000,+,5S2M1I2M,60,1;${supp1.refName},1000,+,10M,60,0;"

    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp1, supp2))) shouldBe 2
    r1[String]("SA") shouldBe s"${supp1.refName},1000,+,5S2M1I2M,60,1;${supp1.refName},1000,+,10M,60,0;"
  }

  it should "not update an SA entry from an alignment of the other read of the pair" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(name="q1", start1=100, start2=300)
    val supp = builder.addFrag(name="q1", start=1000, cigar="5S5M").value
    supp.paired = true; supp.secondOfPair = true; supp.supplementary = true
    r1("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    pretendAligned(supp, "5S5=", nm=0)
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, r2, supp))) shouldBe 1
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S5M,60,0;"
  }

  it should "not update an SA entry from an alignment on the other strand" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100).value
    val supp    = builder.addFrag(name="q1", start=1000, cigar="5S5M", strand=Plus).value
    supp.supplementary = true
    r1("SA") = s"${supp.refName},1000,-,5S5M,60,0;"
    pretendAligned(supp, "5S5=", nm=0)
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp))) shouldBe 1
    r1[String]("SA") shouldBe s"${supp.refName},1000,-,5S5M,60,0;"
  }

  it should "not update an SA entry from the record itself" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100).value
    val supp    = builder.addFrag(name="q1", start=1000, cigar="5S5M").value
    supp.supplementary = true
    supp("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    pretendAligned(supp, "5S5=", nm=0)
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp))) shouldBe 1
    supp[String]("SA") shouldBe s"${supp.refName},1000,+,5S5M,60,0;"
  }

  it should "update the position of SA entries for alignments whose leading indel was trimmed" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="5M5S").value
    val supp1   = builder.addFrag(name="q1", start=1000, cigar="5S2D5M").value
    val supp2   = builder.addFrag(name="q1", start=2000, cigar="5S1I4M").value
    Seq(supp1, supp2).foreach(_.supplementary = true)
    r1("SA") = s"${supp1.refName},1000,+,5S2D5M,60,2;${supp2.refName},2000,+,5S1I4M,60,1;"
    // pretend the leading indels of the supplementary records were trimmed
    pretendAligned(supp1, "5S5M", nm=0, start=1002)
    pretendAligned(supp2, "6S4M", nm=0)
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp1, supp2))) shouldBe 0
    // the edit distance is that of the trimmed alignments
    r1[String]("SA") shouldBe s"${supp1.refName},1002,+,5S5M,60,0;${supp2.refName},2000,+,6S4M,60,0;"
  }

  it should "not update an SA entry for an absent alignment that left-aligning could not have changed into a present one" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="5M5S").value
    val supp    = builder.addFrag(name="q1", start=1002, cigar="5S5M").value
    supp.supplementary = true
    // the entry has the same end and aligned, clipped, and skipped bases as the (unchanged) supplementary alignment
    r1("SA") = s"${supp.refName},1000,+,5S2D5M,60,2;"
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp))) shouldBe 1
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S2D5M,60,2;"
  }

  Seq(
    ("5H2=1I2X1D4M3H", "5S9M3S"),
    ("5H2M3I5M3H", "5S7M3I3S"),
    ("5H2M3D5M3H", "5S7M3D3S"),
    ("5H2M10N5M3H", "5S7M10D3S")
  ).foreach { case (original, summary) =>
    it should s"uniquely match the minimap2 SA summary $summary against original CIGAR $original" in {
      val builder = new SamBuilder(readLength=20)
      val primary = builder.addFrag(name="summary", start=100).value
      val supp = builder.addFrag(name="summary", start=1000, cigar=original,
        bases="A" * Cigar(original).lengthOnQuery, quals="I" * Cigar(original).lengthOnQuery).value
      supp.supplementary = true
      primary("SA") = s"${supp.refName},1000,+,$summary,60,99;"
      pretendAligned(supp, "5H9M3H", nm=2, start=1001)
      IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(primary, supp))) shouldBe 0
      primary[String]("SA") shouldBe s"${supp.refName},1001,+,5S9M3S,60,2;"
    }
  }

  it should "leave ambiguous minimap2 summaries and non-summary length-only matches unchanged and counted" in {
    val builder = new SamBuilder(readLength=20)
    val primary = builder.addFrag(name="summary", start=100).value
    val supp1 = builder.addFrag(name="summary", start=1000, cigar="5S2M1I2M1D4M3S", bases="A" * Cigar("5S2M1I2M1D4M3S").lengthOnQuery, quals="I" * Cigar("5S2M1I2M1D4M3S").lengthOnQuery).value
    val supp2 = builder.addFrag(name="summary", start=1000, cigar="5S3M1D2M1I3M3S", bases="A" * Cigar("5S3M1D2M1I3M3S").lengthOnQuery, quals="I" * Cigar("5S3M1D2M1I3M3S").lengthOnQuery).value
    Seq(supp1, supp2).foreach { supp =>
      supp.supplementary = true
      pretendAligned(supp, "5S9M3S", nm=0)
    }
    val sa = s"${supp1.refName},1000,+,5S9M3S,60,9;${supp1.refName},1000,+,5S4M1I4M1D3S,60,9;"
    primary("SA") = sa
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(primary, supp1, supp2))) shouldBe 2
    primary[String]("SA") shouldBe sa
  }

  it should "not match an absent alignment using the current trimmed alignment's summary" in {
    val builder = new SamBuilder(readLength=10)
    val primary = builder.addFrag(name="summary", start=100).value
    val supp = builder.addFrag(name="summary", start=1000, cigar="5H2M2D3M", bases="A" * Cigar("5H2M2D3M").lengthOnQuery, quals="I" * Cigar("5H2M2D3M").lengthOnQuery).value
    supp.supplementary = true
    pretendAligned(supp, "5H5M", nm=0, start=1002)
    val sa = s"${supp.refName},1002,+,5S5M,60,2;"
    primary("SA") = sa
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(primary, supp))) shouldBe 1
    primary[String]("SA") shouldBe sa
  }

  it should "repair a minimap2 SA summary after actually trimming a supplementary deletion" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname))
    val primary = builder.addFrag(name="summary", start=100, cigar="5M5S").value
    val supp = builder.addFrag(name="summary", start=1000, cigar="5H2M1D3M", bases="AAAAA", quals="IIIII",
      attrs=Map("MD" -> "2^A3")).value
    supp.supplementary = true
    primary("SA") = s"${supp.refName},1000,+,5S5M1D,60,1;"
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output))
    val records = readBamRecs(output)
    val finalSupp = records.find(_.supplementary).value
    finalSupp.start shouldBe 1001
    finalSupp.cigar.toString shouldBe "5H5M"
    finalSupp.contains("NM") shouldBe false
    records.find(!_.supplementary).value[String]("SA") shouldBe s"${supp.refName},1001,+,5S5M,60,0;"
  }

  it should "count exact and summarized candidates together rather than prefer an exact match" in {
    val builder = new SamBuilder(readLength=20)
    val primary = builder.addFrag(name="summary", start=100).value
    val exact = builder.addFrag(name="summary", start=1000, cigar="5S9M3S", bases="A" * Cigar("5S9M3S").lengthOnQuery, quals="I" * Cigar("5S9M3S").lengthOnQuery).value
    val summarized = builder.addFrag(name="summary", start=1000, cigar="5S2M1I2M1D4M3S", bases="A" * Cigar("5S2M1I2M1D4M3S").lengthOnQuery, quals="I" * Cigar("5S2M1I2M1D4M3S").lengthOnQuery).value
    Seq(exact, summarized).foreach(_.supplementary = true)
    val sa = s"${exact.refName},1000,+,5S9M3S,60,9;"
    primary("SA") = sa
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(primary, exact, summarized))) shouldBe 1
    primary[String]("SA") shouldBe sa
  }

  Seq(true, false).foreach { hasScore =>
    it should s"repair mate fields without changing or adding ms tags with existing ms=$hasScore" in {
      val builder = new SamBuilder(readLength=10)
      val pair = builder.addPair(name="scores", start1=100, start2=300)
      val supp = pair.head.clone()
      supp.supplementary = true
      supp.start = 1000
      val secondary = pair.last.clone()
      secondary.secondary = true
      val records = pair ++ Seq(supp, secondary)
      records.zipWithIndex.foreach { case (r, i) =>
        r("AS") = 20 + i
        if (hasScore) r("ms") = 400 + i else r.remove("ms")
      }
      val expected = records.map(_.get[Int]("ms"))
      pair.head("MC") = "1M"
      IndelLeftAligner.fixMateAndSupplementaryInfo(Template(records.iterator))
      records.map(_.get[Int]("ms")) shouldBe expected
      pair.head[String]("MC") shouldBe "10M"
    }
  }

  Seq("2S2D8M" -> true, "2H2D8M" -> false, "2S2D8=" -> false, (null, true)).foreach { case (mc, repair) =>
    it should s"require exact secondary MC agreement when MC=$mc" in {
      val builder = new SamBuilder(readLength=10)
      val Seq(r1, r2) = builder.addPair(name="clipping", start1=100, start2=300, cigar2="2S2D8M")
      val secondary = r1.clone()
      secondary.secondary = true
      secondary.start = 2000
      secondary("MC") = mc
      secondary.insertSize = 777
      pretendAligned(r2, "2S8M", nm=0, start=302)
      IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, r2, secondary)))
      secondary.mateStart shouldBe (if (repair) 302 else 300)
      secondary.get[String]("MC") shouldBe (if (repair) Some("2S8M") else Option(mc))
      if (repair) secondary.insertSize should not be 777 else secondary.insertSize shouldBe 777
    }
  }

  Seq(false, true).foreach { negativeMate =>
    it should s"repair eligible secondary mate fields after trimming with negativeMate=$negativeMate" in {
      val builder = new SamBuilder(readLength=10)
      val Seq(r1, r2) = builder.addPair(name="q1", start1=100, start2=300, cigar1="10M", cigar2="2D10M")
      r2.negativeStrand = negativeMate
      def secondary(mateStart: Int, mateCigar: String): SamRecord = {
        val sec = builder.addFrag(name="q1", start=2000, cigar="10M").value
        sec.paired = true; sec.firstOfPair = true; sec.secondary = true
        sec.mateRefIndex = r2.refIndex; sec.mateStart = mateStart; sec.mateNegativeStrand = r2.negativeStrand
        sec("MC") = mateCigar
        sec
      }
      val sec    = secondary(mateStart=300, mateCigar="2D10M")
      val noMc   = secondary(mateStart=300, mateCigar=null) // no mate cigar, so only its mate start must match R2's
      val decoy1 = secondary(mateStart=310, mateCigar="10M") // its mate start does not match R2's
      val decoy2 = secondary(mateStart=300, mateCigar="10M") // its mate cigar does not match R2's
      val opposite = secondary(mateStart=300, mateCigar="2D10M")
      opposite.mateNegativeStrand = !r2.negativeStrand
      val oppositeNoMc = secondary(mateStart=300, mateCigar=null)
      oppositeNoMc.mateNegativeStrand = !r2.negativeStrand
      // pretend R2's leading deletion was trimmed
      pretendAligned(r2, "10M", nm=0, start=302)
      IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, r2, sec, noMc, decoy1, decoy2, opposite, oppositeNoMc)))
      Seq(sec, noMc).foreach { r =>
        r.mateStart shouldBe 302
        r.insertSize shouldBe htsjdk.samtools.SamPairUtil.computeInsertSize(r.asSam, r2.asSam)
        r.insertSize should not be 0
        r[String]("MC") shouldBe "10M"
      }
      decoy1.mateStart shouldBe 310
      decoy2.mateStart shouldBe 300
      decoy2[String]("MC") shouldBe "10M"
      Seq(opposite, oppositeNoMc).foreach { r =>
        r.mateStart shouldBe 300
        r.mateNegativeStrand shouldBe !r2.negativeStrand
        r.insertSize shouldBe 0
      }
      opposite[String]("MC") shouldBe "2D10M"
      oppositeNoMc.contains("MC") shouldBe false
    }
  }

  behavior of "LeftAlignIndels"

  // A non-repetitive reference with a dinucleotide repeat (ACACACAC) at 0-based offsets 200-207, preceded by a G and
  // followed by a T, and with a soft-masked (lower-case) stretch after it that has a homopolymer (tttt) at 0-based
  // offsets 260-263, preceded by a g and followed by a c.
  private val RefBases: String = {
    val random = new scala.util.Random(42)
    def randomBases(n: Int): String = {
      val bases = Iterator.continually("ACGT"(random.nextInt(4))).sliding(2).collect { case Seq(a, b) if a != b => b }
      bases.take(n).mkString
    }
    val left  = randomBases(199) + "G"
    val right = "T" + randomBases(50).toLowerCase + "gttttc" + randomBases(43).toLowerCase + randomBases(500)
    left + "ACACACAC" + right
  }

  /** The reference bases in upper-case, for computing the expected MD and NM tags. */
  private val RefBytesUpper: Array[Byte] = RefBases.toUpperCase.getBytes

  private case class Inputs(sam: SamBuilder, ref: Path)

  /** Builds a reference FASTA and reads: a fragment with no indels, a fragment with a deletion in the repeat, a
    * fragment with a deletion in the soft-masked homopolymer, a pair with an insertion in the repeat, and an unmapped
    * fragment. */
  private def inputs(samOrder: SamOrder, md: Boolean = true): Inputs = {
    val ref = fasta(Seq("chr1" -> RefBases))

    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=100, sort=Some(samOrder), sd=Some(dict))
    def bases(from: Int, until: Int): String = RefBases.substring(from, until).toUpperCase
    // no indels
    builder.addFrag(name="frag1", bases=bases(10, 110), start=11, cigar="100M")
    // deletion of the last "AC" in the repeat (offsets 206-207), read starts at offset 150
    builder.addFrag(name="frag2", bases=bases(150, 206) + bases(208, 252), start=151, cigar="56M2D44M")
    // insertion of "AC" after the repeat (after offset 207), read starts at offset 150
    builder.addPair(name="pair1", bases1=bases(150, 208) + "AC" + bases(208, 248), start1=151, cigar1="58M2I40M",
      bases2=bases(300, 400), start2=301, cigar2="100M", strand1=Plus, strand2=Minus)
    // unmapped
    builder.addFrag(name="frag3", bases=bases(10, 110), unmapped=true)
    // deletion of the third "t" in the soft-masked homopolymer (offset 262), read starts at offset 230
    builder.addFrag(name="frag4", bases=bases(230, 262) + bases(263, 331), start=231, cigar="32M1D68M")
    if (md) {
      builder.iterator.filter(_.mapped).foreach(rec => SequenceUtil.calculateMdAndNmTags(rec.asSam, RefBytesUpper, true, true))
    }
    Inputs(builder, ref)
  }

  Seq(SamOrder.Coordinate, SamOrder.Queryname, SamOrder.RandomQuery, SamOrder.Random).foreach { samOrder =>
    Seq(true, false).foreach { useRef =>
      val testName = f"on a $samOrder sorted BAM " + (if (useRef) "using the ref" else "using the MD tag")
      it should f"run end to end $testName" in {
        val in     = inputs(samOrder)
        val output = makeTempFile("test.", ".bam")
        val tool   = new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=if (useRef) Some(in.ref) else None, useEqualsAndX=true)
        executeFgbioTool(tool)

        // the output is in the input's order
        val source  = SamSource(output)
        val ordered = source.toIndexedSeq
        source.safelyClose()
        SamOrder(source.header) shouldBe Some(samOrder)
        samOrder match {
          case SamOrder.Coordinate => ordered.map(_.name) shouldBe ordered.sortBy(SamOrder.Coordinate.sortkey).map(_.name)
          case SamOrder.Queryname  => ordered.map(_.name) shouldBe ordered.sortBy(SamOrder.Queryname.sortkey).map(_.name)
          case _                   => ()
        }

        val records = ordered.map(r => (r.name, r.paired && r.firstOfPair, r)).sortBy(t => (t._1, !t._2)).map(_._3)
        records.map(_.name) shouldBe Seq("frag1", "frag2", "frag3", "frag4", "pair1", "pair1")
        val Seq(frag1, frag2, frag3, frag4, r1, r2) = records
        frag1.start shouldBe 11
        frag1.cigar.toString shouldBe "100="
        frag2.start shouldBe 151
        frag2.cigar.toString shouldBe "50=2D50="
        frag3.unmapped shouldBe true
        frag4.start shouldBe 231
        frag4.cigar.toString shouldBe "30=1D70="
        r1.cigar.toString shouldBe "50=2I48="
        r2.cigar.toString shouldBe "100="
        r1.mateCigar.value.toString shouldBe "100="
        r2.mateCigar.value.toString shouldBe "50=2I48="

        // the MD and NM tags should be consistent with the left-aligned CIGARs
        records.filter(_.mapped).foreach { rec =>
          val expected = rec.asSam.deepCopy()
          SequenceUtil.calculateMdAndNmTags(expected, RefBytesUpper, true, true)
          withClue(rec.name) {
            rec[String]("MD") shouldBe expected.getStringAttribute("MD")
            rec[Int]("NM") shouldBe expected.getIntegerAttribute("NM").intValue()
          }
        }
      }
    }
  }

  it should "use M operators when useEqualsAndX is false" in {
    val in     = inputs(SamOrder.Coordinate)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    val cigars = readBamRecs(output).filter(_.mapped).map(r => (r.name, r.paired && r.firstOfPair) -> r.cigar.toString).toMap
    cigars shouldBe Map(
      ("frag1", false) -> "100M",
      ("frag2", false) -> "50M2D50M",
      ("frag4", false) -> "30M1D70M",
      ("pair1", true)  -> "50M2I48M",
      ("pair1", false) -> "100M"
    )
  }

  it should "accept a plain FASTA without a standard filename extension" in {
    val in = inputs(SamOrder.Queryname, md=false)
    val ref = makeTempFile("test.", ".txt")
    Files.copy(in.ref, ref, StandardCopyOption.REPLACE_EXISTING)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
    readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
  }

  Seq(false, true).foreach { useRef =>
    it should s"omit skipped bases from CLI reference lookup and preserve spliced alignments with FASTA=$useRef" in {
      val skip = 1000000
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=skip + 106))
      val builder = new SamBuilder(readLength=7, sort=Some(SamOrder.Queryname), sd=Some(dict))
      val rec = builder.addFrag(name="spliced", bases="CAAAAAG", start=1, cigar=s"2M${skip}N2M1I1M100N1M",
        attrs=Map("MD" -> "6", "NM" -> 1, "UQ" -> 0)).value
      val ref = if (useRef) Some(fasta(Seq("chr1" -> ("CA" + "N" * skip + "AAA" + "N" * 100 + "G")))) else None
      val input = builder.toTempFile()
      val output = makeTempFile("test.", ".bam")
      val tool = new LeftAlignIndels(input=input, output=output, ref=ref, leadingIndelMode=LeadingIndelMode.Keep)
      tool.refBasesFor(builder.header)(rec) shouldBe "CAAAAG"
      executeFgbioTool(tool)
      val records = readBamRecs(output)
      records.size shouldBe 1
      val result = records.head
      result.cigar.toString shouldBe s"2M${skip}N1I3M100N1M"
      result.basesString shouldBe "CAAAAAG"
      result[String]("MD") shouldBe "6"
      result[Int]("NM") shouldBe 1
      result[Int]("UQ") shouldBe 0
    }
  }

  Seq(false, true).foreach { useRef =>
    it should s"trim and regenerate tags across multiple skips and clipping with FASTA=$useRef" in {
      val target = "AAG" + "N" * 5 + "CGG" + "N" * 3 + "AT"
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=target.length))
      val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
      val bases = if (useRef) "=AGCCGG=T=" else "=AGCCGGAT="
      builder.addFrag(name="spliced", bases=bases, start=1, cigar="1S1D2M5N1M1I2M3N2M1S",
        attrs=Map("MD" -> "0^A7", "NM" -> 2, "UQ" -> 0))
      val ref = if (useRef) Some(fasta(Seq("chr1" -> target))) else None
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=ref, useEqualsAndX=true))
      val result = readBamRecs(output).head
      result.start shouldBe 2
      result.cigar.toString shouldBe "1S2=5N1I3=3N2=1S"
      result.basesString shouldBe "=AGCCGGAT="
      result[String]("MD") shouldBe "7"
      result[Int]("NM") shouldBe 1
      result[Int]("UQ") shouldBe 0
    }
  }

  it should "report a read-specific error for an unrepresentable resolved FASTA base" in {
    val builder = new SamBuilder(readLength=3, sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=3))))
    builder.addFrag(name="uracil", bases="=AA", cigar="3M", start=1)
    val output = makeTempFile("test.", ".bam")
    val ref = fasta(Seq("chr1" -> "UAA"))
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref), useEqualsAndX=true))
    }
    ex.message.value should include ("uracil")
    ex.message.value should include ("Cannot store the resolved '=' read bases unchanged")
  }

  Seq(false, true).foreach { useRef =>
    it should s"left-align across a circular origin with FASTA=$useRef" in {
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=10, topology=Some(Topology.Circular)))
      val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Queryname), sd=Some(dict))
      val rec = builder.addFrag(name="circular", bases="TACG", start=9, cigar="2M1D2M",
        attrs=Map("MD" -> "2^A2", "NM" -> 1)).value
      val ref = if (useRef) Some(fasta(Seq("chr1" -> "ACGTCAGTTA"))) else None
      val input = builder.toTempFile()
      readBamRecs(input).size shouldBe 1
      val output = makeTempFile("test.", ".bam")
      val tool = new LeftAlignIndels(input=input, output=output, ref=ref,
        useEqualsAndX=true, leadingIndelMode=LeadingIndelMode.Keep)
      tool.refBasesFor(builder.header)(rec) shouldBe "TAACG"
      executeFgbioTool(tool)
      val result = readBamRecs(output).head
      result.start shouldBe 9
      result.cigar.toString shouldBe "1=1D3="
      result[String]("MD") shouldBe "1^A3"
      result[Int]("NM") shouldBe 1
    }

    it should s"wrap a trimmed circular alignment start with FASTA=$useRef" in {
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=10, topology=Some(Topology.Circular)))
      val builder = new SamBuilder(readLength=2, sort=Some(SamOrder.Queryname), sd=Some(dict))
      builder.addFrag(name="circular-trim", bases="AC", start=10, cigar="1M1D1M",
        attrs=Map("MD" -> "1^A1", "NM" -> 1))
      val ref = if (useRef) Some(fasta(Seq("chr1" -> "ACGTACGTAA"))) else None
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=ref, useEqualsAndX=true))
      val result = readBamRecs(output).head
      result.start shouldBe 1
      result.cigar.toString shouldBe "2="
      result[String]("MD") shouldBe "2"
      result[Int]("NM") shouldBe 0
    }
  }

  it should "skip multiple circular turns without allocating the skipped reference bases" in {
    val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=10, topology=Some(Topology.Circular)))
    val builder = new SamBuilder(readLength=6, sort=Some(SamOrder.Queryname), sd=Some(dict))
    val rec = builder.addFrag(name="circular-skip", bases="TAAAGT", start=9, cigar="2M25N1M1I2M",
      attrs=Map("MD" -> "5", "NM" -> 1)).value
    val ref = fasta(Seq("chr1" -> "ACGTCAGTTA"))
    val output = makeTempFile("test.", ".bam")
    val tool = new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref), leadingIndelMode=LeadingIndelMode.Keep)
    tool.refBasesFor(builder.header)(rec) shouldBe "TAAGT"
    executeFgbioTool(tool)
    val result = readBamRecs(output).head
    result.cigar.toString shouldBe "2M25N1I3M"
    result[String]("MD") shouldBe "5"
    result[Int]("NM") shouldBe 1
  }

  it should "wrap several times within a single reference-consuming CIGAR element" in {
    val contig = "ACGTCAGTTA"
    val target = (contig * 4).substring(8, 33)
    val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=10, topology=Some(Topology.Circular)))
    val builder = new SamBuilder(readLength=25, sort=Some(SamOrder.Queryname), sd=Some(dict))
    val rec = builder.addFrag(name="circular-turns", bases=target, start=9, cigar="25M",
      attrs=Map("MD" -> "25", "NM" -> 0)).value
    val ref = fasta(Seq("chr1" -> contig))
    val output = makeTempFile("test.", ".bam")
    val tool = new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref), useEqualsAndX=true)
    tool.refBasesFor(builder.header)(rec) shouldBe target
    executeFgbioTool(tool)
    val result = readBamRecs(output).head
    result.start shouldBe 9
    result.cigar.toString shouldBe "25="
    result[String]("MD") shouldBe "25"
    result[Int]("NM") shouldBe 0
  }

  it should "still reject an alignment past the end of a linear contig" in {
    val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Queryname),
      sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=10))))
    val rec = builder.addFrag(name="linear", bases="TACG", start=9, cigar="2M1D2M").value
    val ref = fasta(Seq("chr1" -> "ACGTCAGTTA"))
    val tool = new LeftAlignIndels(input=builder.toTempFile(), output=makeTempFile("test.", ".bam"), ref=Some(ref))
    intercept[FailureException] { tool.refBasesFor(builder.header)(rec) }.message.value should include ("past the end")
  }

  it should "fail contextually when the reference is not valid FASTA" in {
    val in = inputs(SamOrder.Queryname)
    val ref = makeTempFile("test.", ".txt")
    Files.write(ref, "This is not a FASTA file.\n".getBytes)
    val output = makeTempFile("test.", ".bam")
    Files.write(output, "previous output".getBytes)
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
    }
    ex.message.value should include ("Could not read reference FASTA")
    ex.message.value should include (ref.toString)
    new String(Files.readAllBytes(output)) shouldBe "previous output"
  }

  it should "fail with a helpful message when the MD tag is needed but missing" in {
    val in     = inputs(SamOrder.Queryname, md=false)
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output))
    }
    ex.message.value should include ("has no MD tag")
  }

  Seq(SamOrder.Coordinate, SamOrder.Queryname, SamOrder.RandomQuery, SamOrder.Random, SamOrder.Unsorted, SamOrder.TemplateCoordinate).foreach { order =>
    it should s"use a sidecarless reference with $order input without MD tags" in {
      val in     = inputs(order, md=false)
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
      readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
    }
  }

  /** A well-formed sequence dictionary whose chr1 length (1) differs from the reference, which would truncate chr1 if
    * it were used when reading the FASTA. */
  private val StaleDictionary: String = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1\n"

  /** Writes a reference FASTA with the given contigs and no sidecars. */
  private def fasta(contigs: Seq[(String, String)]): Path = {
    val ref = makeTempFile("test.", ".fa")
    Files.write(ref, contigs.map { case (name, bases) => s">$name\n$bases\n" }.mkString.getBytes)
    ref
  }

  /** Compresses a generated FASTA without creating index or dictionary sidecars. */
  private def compressedFasta(ref: Path, bgzf: Boolean): Path = {
    val compressed = makeTempFile("test.", ".fa.gz")
    val stream = if (bgzf) new BlockCompressedOutputStream(compressed)
                 else new GZIPOutputStream(Files.newOutputStream(compressed))
    Using.resource(stream)(_.write(Files.readAllBytes(ref)))
    compressed
  }

  it should "fail when the reference FASTA is missing a contig in the BAM's sequence dictionary" in {
    val in     = inputs(SamOrder.Queryname)
    val output = makeTempFile("test.", ".bam")
    val ref    = fasta(Seq("chr2" -> RefBases))
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
    }
    ex.message.value should include ("Contig chr1 is not in the reference FASTA")
  }

  it should "fail when the reference FASTA has a contig length that differs from the BAM's sequence dictionary" in {
    val in     = inputs(SamOrder.Queryname)
    val output = makeTempFile("test.", ".bam")
    val ref    = fasta(Seq("chr1" -> (RefBases + "ACGT")))
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
    }
    ex.message.value should include ("Contig chr1 has length")
  }

  it should "fail with a helpful message when a read extends past the end of its contig" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Queryname), sd=Some(dict))
    builder.addFrag(name="frag1", bases="A"*100, start=RefBases.length - 50, cigar="50M1D50M")
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(fasta(Seq("chr1" -> RefBases)))))
    }
    ex.message.value should include ("past the end of contig chr1")
  }
  it should "use a reference FASTA whose contigs are in a different order than the BAM's sequence dictionary" in {
    val in     = inputs(SamOrder.Coordinate, md=false)
    val output = makeTempFile("test.", ".bam")
    val ref    = fasta(Seq("chr0" -> "ACGT", "chr1" -> RefBases))
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
    readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
  }

  it should "write a query name sorted output when the input order is not known" in {
    val in     = inputs(SamOrder.Unsorted)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output))
    val source = SamSource(output)
    SamOrder(source.header) shouldBe Some(SamOrder.Queryname)
    source.safelyClose()
  }

  it should "not re-sort an input that is query name sorted with a sub-sort fgbio does not know" in {
    val in     = inputs(SamOrder.Queryname)
    // samtools' natural query name sub-sort, which SamOrder does not model
    val header = SamOrder.Queryname.applyTo(in.sam.header.clone())
    header.setAttribute("SS", "queryname:natural")
    SamOrder(header) shouldBe None
    // read2 precedes read10 in natural order, but not in lexicographic order
    val read2  = in.sam.addFrag(name="read2", unmapped=true).value
    val read10 = in.sam.addFrag(name="read10", unmapped=true).value
    val records = in.sam.toIndexedSeq.filterNot(r => r.name.startsWith("read")).sortBy(SamOrder.Queryname.sortkey) ++ Seq(read2, read10)
    val input = writeRecords(header, records)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output))
    val source = SamSource(output)
    source.header.getAttribute("SS") shouldBe "queryname:natural"
    source.safelyClose()
    readBamRecs(output).map(_.name) shouldBe Seq("frag1", "frag2", "frag3", "frag4", "pair1", "pair1", "read2", "read10")
  }

  it should "keep the input order of a template-coordinate sorted input" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sd=Some(dict))
    // R2 is at the lower coordinate, so precedes R1 in template-coordinate order
    val Seq(r1, r2) = builder.addPair(name="pair1", start1=101, start2=11, strand1=Minus, strand2=Plus)
    val header = SamOrder.TemplateCoordinate.applyTo(builder.header.clone())
    val input = writeRecords(header, Seq(r2, r1))

    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output))
    val source = SamSource(output)
    SamOrder(source.header) shouldBe Some(SamOrder.TemplateCoordinate)
    source.safelyClose()
    readBamRecs(output).map(r => (r.name, r.firstOfPair)) shouldBe Seq(("pair1", false), ("pair1", true))
  }

  it should "not group records without a query name ('*')" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
    builder.addFrag(name="*", unmapped=true)
    builder.addFrag(name="*", unmapped=true)
    builder.addPair(name="pair1", start1=11, start2=101)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output))
    readBamRecs(output).map(_.name) shouldBe Seq("*", "*", "pair1", "pair1")
  }

  it should "fail when the output is the same file as the input" in {
    val in    = inputs(SamOrder.Coordinate)
    val input = in.sam.toTempFile()
    an[ValidationException] shouldBe thrownBy { new LeftAlignIndels(input=input, output=input) }
  }

  it should "fail with a helpful message when using the MD tag with '=' read bases" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Coordinate), sd=Some(dict))
    builder.addFrag(name="frag1", bases="CA=G", cigar="3M1D1M", start=11, attrs=Map("MD" -> "3^A1"))
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output))
    }
    ex.message.value should include ("for an '=' read base")
    ex.message.value should include ("Provide --ref")
  }

  Seq(false, true).foreach { secondary =>
    it should s"fail on a paired read without its mate unless allowed, with secondary=$secondary" in {
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
      val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
      val Seq(r1, _) = builder.addPair(name="pair1", start1=11, start2=101)
      r1.secondary = secondary
      val input = writeRecords(builder.header, Seq(r1))
      val output = makeTempFile("test.", ".bam")
      val ex = intercept[FailureException] { executeFgbioTool(new LeftAlignIndels(input=input, output=output)) }
      ex.message.value should include ("is paired but its mate was not found")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true))
      readBamRecs(output).map(_.name) shouldBe Seq("pair1")
    }
  }

  it should "not use an unpaired record with the same query name as the mate of a paired read" in {
    val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
    val Seq(_, r2) = builder.addPair(name="q", start1=11, start2=101)
    val unpaired = builder.addFrag(name="q", start=201, cigar="5M5S").value
    val (mateStart, mateCigar) = (r2.mateStart, r2.get[String]("MC"))
    val input = writeRecords(builder.header, Seq(unpaired, r2))
    val output = makeTempFile("test.", ".bam")
    val ex = intercept[FailureException] { executeFgbioTool(new LeftAlignIndels(input=input, output=output)) }
    ex.message.value should include ("is paired but its mate was not found")

    executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true))
    val records = readBamRecs(output)
    records should have size 2
    val outR2 = records.find(_.paired).value
    outR2.mateStart shouldBe mateStart
    outR2.get[String]("MC") shouldBe mateCigar
    records.find(!_.paired).value.start shouldBe 201
  }

  /** Builds a reference FASTA and coordinate sorted reads: a pair whose R1 has a deletion in the repeat that is shifted
    * to the start of its alignment, and a fragment with a supplementary alignment whose deletion is shifted. */
  private def leadingIndelInputs(samOrder: SamOrder = SamOrder.Coordinate): Inputs = {
    val ref     = fasta(Seq("chr1" -> RefBases))
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=100, sort=Some(samOrder), sd=Some(dict))
    def bases(from: Int, until: Int): String = RefBases.substring(from, until).toUpperCase
    // deletion of "AC" at offsets 204-205 of the repeat, read starts at offset 200 (the start of the repeat)
    builder.addPair(name="pair2", bases1=bases(200, 204) + bases(206, 302), start1=201, cigar1="4M2D96M",
      bases2=bases(400, 500), start2=401, cigar2="100M", strand1=Plus, strand2=Minus)
    // a primary alignment of the first 50 bases, and a supplementary alignment of the last 50 bases with a deletion of
    // "AC" at offsets 204-205 of the repeat
    val suppBases = bases(180, 204) + bases(206, 232)
    val primary   = builder.addFrag(name="frag5", bases=bases(500, 550) + suppBases, start=501, cigar="50M50S").value
    val supp      = builder.addFrag(name="frag5", bases=suppBases, quals="I"*50, start=181, cigar="50H24M2D26M").value
    supp.supplementary = true
    primary("SA") = "chr1,181,+,50S24M2D26M,60,2;"
    supp("SA")    = "chr1,501,+,50M50S,60,0;"
    builder.iterator.filter(_.mapped).foreach(rec => SequenceUtil.calculateMdAndNmTags(rec.asSam, RefBytesUpper, true, true))
    Inputs(builder, ref)
  }

  it should "trim a leading deletion by default, updating the mate information and tags" in {
    val in     = leadingIndelInputs()
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    val records = readBamRecs(output)
    val r1      = records.find(r => r.name == "pair2" && r.firstOfPair).value
    val r2      = records.find(r => r.name == "pair2" && r.secondOfPair).value
    r1.start shouldBe 203
    r1.cigar.toString shouldBe "100M"
    r2.mateStart shouldBe 203
    r2[String]("MC") shouldBe "100M"
    records.filter(_.mapped).foreach { rec =>
      val expected = rec.asSam.deepCopy()
      SequenceUtil.calculateMdAndNmTags(expected, RefBytesUpper, true, true)
      withClue(rec.name) {
        rec[String]("MD") shouldBe expected.getStringAttribute("MD")
        rec[Int]("NM") shouldBe expected.getIntegerAttribute("NM").intValue()
      }
    }
    // the output is still coordinate sorted
    records.map(_.start) shouldBe records.map(_.start).sorted
  }

  it should "keep a leading deletion when the leading indel mode is Keep" in {
    val in     = leadingIndelInputs()
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref), leadingIndelMode=LeadingIndelMode.Keep))
    val r1 = readBamRecs(output).find(r => r.name == "pair2" && r.firstOfPair).value
    r1.start shouldBe 201
    r1.cigar.toString shouldBe "2D100M"
  }

  it should "remove every record of a template with a leading indel when the leading indel mode is Filter" in {
    val in     = leadingIndelInputs()
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref), leadingIndelMode=LeadingIndelMode.Filter))
    readBamRecs(output).map(_.name) shouldBe Seq("frag5", "frag5")
  }

  it should "update the CIGARs in SA tags end to end" in {
    val in     = leadingIndelInputs()
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    val Seq(primary, supp) = readBamRecs(output).filter(_.name == "frag5").sortBy(_.supplementary)
    supp.cigar.toString shouldBe "50H20M2D30M"
    primary[String]("SA") shouldBe "chr1,181,+,50S20M2D30M,60,2;"
    supp[String]("SA") shouldBe "chr1,501,+,50M50S,60,0;"
  }

  it should "re-sort a template-coordinate input when trimming leading indels" in {
    val in     = leadingIndelInputs()
    // a pair whose lower coordinate is between that of pair2 before (201) and after (203) its leading deletion is trimmed
    val pair3  = in.sam.addPair(name="pair3", bases1=RefBases.substring(201, 301).toUpperCase, start1=202, cigar1="100M",
      bases2=RefBases.substring(400, 500).toUpperCase, start2=401, cigar2="100M", strand1=Plus, strand2=Minus)
    pair3.foreach(r => SequenceUtil.calculateMdAndNmTags(r.asSam, RefBytesUpper, true, true))
    val header = SamOrder.TemplateCoordinate.applyTo(in.sam.header.clone())
    val input = writeRecords(header, in.sam.toIndexedSeq.sortBy(SamOrder.TemplateCoordinate.sortkey))
    readBamRecs(input).map(_.id).filter(_.startsWith("pair")) shouldBe Seq("pair2/1", "pair2/2", "pair3/1", "pair3/2")

    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output, ref=Some(in.ref)))
    val records = readBamRecs(output)
    records.find(r => r.name == "pair2" && r.firstOfPair).value.start shouldBe 203
    records.map(_.id).filter(_.startsWith("pair")) shouldBe Seq("pair3/1", "pair3/2", "pair2/1", "pair2/2")
    records.map(_.id) shouldBe records.sortBy(SamOrder.TemplateCoordinate.sortkey).map(_.id)
  }

  it should "write unpaired records that share a query name without updating their SA tags" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate), sd=Some(dict))
    val sa = "chr1,101,+,10M,60,0;"
    builder.addFrag(name="dup", start=11, attrs=Map("SA" -> sa))
    builder.addFrag(name="dup", start=101)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output))
    val records = readBamRecs(output)
    records.map(_.name) shouldBe Seq("dup", "dup")
    records.find(_.start == 11).value[String]("SA") shouldBe sa
  }

  behavior of "LeftAlignIndels"

  it should "fail with a helpful message when paired records share a query name with another pair" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate), sd=Some(dict))
    builder.addPair(name="dup", start1=11, start2=101)
    builder.addPair(name="dup", start1=201, start2=301)
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] { executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output)) }
    ex.message.value should include ("Could not group the records with query name dup")
  }

  for (protectedInput <- Seq(false, true); alias <- Seq("direct", "hard link", "symlink")) {
    it should s"reject a selected BAM index collision with input=$protectedInput through $alias" in {
      val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Coordinate),
        sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=4))))
      val originalInput = builder.toTempFile()
      val protectedPath = if (protectedInput) originalInput else fasta(Seq("chr1" -> "ACGT"))
      val output = makeTempFile("test.", ".bam")
      val index = PathUtil.replaceExtension(output, ".bai")
      val source = alias match {
        case "direct" => Files.copy(protectedPath, index)
        case "hard link" => Files.createLink(index, protectedPath); protectedPath
        case "symlink" => Files.createSymbolicLink(index, protectedPath); protectedPath
      }
      val originalBytes = Files.readAllBytes(source)
      val oldOutput = "previous output".getBytes
      Files.write(output, oldOutput)
      val input = if (protectedInput) source else originalInput
      val ref = if (protectedInput) None else Some(source)
      val ex = intercept[FailureException] { executeFgbioTool(new LeftAlignIndels(input=input, output=output, ref=ref)) }
      ex.message.value should include ("output index")
      Files.readAllBytes(source) shouldBe originalBytes
      Files.readAllBytes(output) shouldBe oldOutput
    }
  }

  Seq(".data" -> ".bai", ".sam.gz" -> ".csi", ".cram" -> ".crai").foreach { case (suffix, indexSuffix) =>
    it should s"reject the selected $suffix writer index when it is the reference" in {
      val builder = new SamBuilder(sort=Some(SamOrder.Coordinate),
        sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=4))))
      val output = makeTempFile("test.", suffix)
      val index = output.resolveSibling(output.getFileName.toString + indexSuffix)
      Files.copy(fasta(Seq("chr1" -> "ACGT")), index)
      val original = Files.readAllBytes(index)
      val ex = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(index)))
      }
      ex.message.value should include ("output index")
      Files.readAllBytes(index) shouldBe original
    }
  }

  it should "reject the selected index of an output named only .bam when it is the reference" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate),
      sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=4))))
    val dir = Io.makeTempDir("left-align-indels")
    val output = dir.resolve(".bam")
    val index = dir.resolve(".bai") // HTSJDK strips the .bam suffix, leaving only the index extension
    Files.copy(fasta(Seq("chr1" -> "ACGT")), index)
    val original = Files.readAllBytes(index)
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(index)))
    }
    ex.message.value should include ("output index")
    Files.readAllBytes(index) shouldBe original
    Files.exists(output) shouldBe false
  }

  it should "leave an unselected alternate index name alone during ordinary indexed output" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate),
      sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=4))))
    val output = makeTempFile("test.", ".bam")
    val reference = output.resolveSibling(output.getFileName.toString + ".bai")
    Files.copy(fasta(Seq("chr1" -> "ACGT")), reference)
    val original = Files.readAllBytes(reference)
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(reference)))
    Files.readAllBytes(reference) shouldBe original
    Files.exists(PathUtil.replaceExtension(output, ".bai")) shouldBe true
  }

  Seq(false, true).foreach { longContig =>
    it should s"allow an input at the would-be index path when indexing is disabled, longContig=$longContig" in {
      val length = if (longContig) 600000000 else 4
      val order = if (longContig) SamOrder.Coordinate else SamOrder.Queryname
      val builder = new SamBuilder(sort=Some(order), sd=Some(SequenceDictionary(SequenceMetadata(name="chr1", length=length))))
      val output = makeTempFile("test.", ".bam")
      val input = PathUtil.replaceExtension(output, ".bai")
      Files.copy(builder.toTempFile(), input)
      val original = Files.readAllBytes(input)
      executeFgbioTool(new LeftAlignIndels(input=input, output=output))
      Files.readAllBytes(input) shouldBe original
    }
  }

  it should "fail when the output is the same file as the reference FASTA" in {
    val in = inputs(SamOrder.Coordinate)
    an[ValidationException] shouldBe thrownBy { new LeftAlignIndels(input=in.sam.toTempFile(), output=in.ref, ref=Some(in.ref)) }
  }

  it should "filter only the unpaired records with a leading indel that share a query name with another unpaired record" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate), sd=Some(dict))
    builder.addFrag(name="dup", start=11, cigar="1D10M", attrs=Map("MD" -> "0^A10"))
    builder.addFrag(name="dup", start=101)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, leadingIndelMode=LeadingIndelMode.Filter))
    readBamRecs(output).map(r => (r.name, r.start)) shouldBe Seq(("dup", 101))
  }

  Seq(SamOrder.Coordinate, SamOrder.Queryname).foreach { order =>
    it should s"preserve genuine lA tags on changed and unchanged records with $order input" in {
      val in = leadingIndelInputs(order)
      in.sam.iterator.foreach(r => r("lA") = s"user:${r.name}:${r.start}:${r.cigar}")
      val expected = in.sam.iterator.map(r => (r.name, r.flags) -> r[String]("lA")).toMap
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
      val records = readBamRecs(output)
      records.map(r => (r.name, r.flags) -> r[String]("lA")).toMap shouldBe expected
      records.find(r => r.name == "pair2" && r.firstOfPair).value.start shouldBe 203
      records.find(r => r.name == "frag5" && !r.supplementary).value[String]("SA") shouldBe "chr1,181,+,50S20M2D30M,60,2;"
    }
  }

  it should "repair SA summaries when keeping leading indels" in {
    val in     = leadingIndelInputs()
    in.sam.iterator.find(r => r.name == "frag5" && !r.supplementary).value("SA") = "chr1,181,+,50S50M2D,60,2;"
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref), leadingIndelMode=LeadingIndelMode.Keep))
    val records = readBamRecs(output)
    // the SA entry of the supplementary alignment whose deletion was shifted is updated
    records.find(r => r.name == "frag5" && !r.supplementary).value[String]("SA") shouldBe "chr1,181,+,50S20M2D30M,60,2;"
  }

  for (hasNm <- Seq(true, false); samtoolsMismatches <- Seq(false, true)) {
    it should s"write the re-computed edit distance to SA entries for a changed alignment with NM present=$hasNm samtoolsMismatches=$samtoolsMismatches" in {
      // An N/N pair is a match by default (as in HTSJDK), but a mismatch with samtoolsMismatches (as in samtools)
      val nm      = if (samtoolsMismatches) 2 else 1
      val ref     = "GGGGGGGGGGTNAAC" + "G" * 85
      val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=ref.length))
      val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Coordinate), sd=Some(dict))
      val primary = builder.addFrag(name="frag", bases="TNAC", start=11, cigar="1M3S").value
      val supp    = builder.addFrag(name="frag", bases="NAC", quals="III", start=12, cigar="1H2M1D1M").value
      supp.supplementary = true
      primary("SA") = "chr1,12,+,1S2M1D1M,60,1;"
      supp("SA")    = "chr1,11,+,1M3S,60,0;"
      if (hasNm) supp("NM") = 1
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(fasta(Seq("chr1" -> ref))),
        samtoolsMismatches=samtoolsMismatches))
      val records = readBamRecs(output)
      val outSupp = records.find(_.supplementary).value
      outSupp.cigar.toString shouldBe "1H1M1D2M"
      if (hasNm) outSupp[Int]("NM") shouldBe nm else outSupp.contains("NM") shouldBe false
      records.find(!_.supplementary).value[String]("SA") shouldBe s"chr1,12,+,1S1M1D2M,60,$nm;"
    }
  }

  /** Writes records in their supplied order without a sorter. */
  private def writeRecords(header: SAMFileHeader, records: Iterable[SamRecord]): Path = {
    val input = makeTempFile("test.", ".bam")
    val writer = SamWriter(input, header)
    writer ++= records
    writer.close()
    input
  }

  /** Includes alignment identity when checking order, not just the query name. */
  private def identities(records: Seq[SamRecord]): Seq[(String, Int, Int, String)] =
    records.map(r => (r.name, r.flags, r.start, r.cigar.toString))

  private def assertSorted(records: Seq[SamRecord], order: SamOrder): Unit = {
    identities(records) shouldBe identities(records.sortBy(order.sortkey))
    ()
  }

  // The reference is loaded before, and independently of, the input order, which other tests cover with a reference
  Seq("plain", "gzip", "BGZF").foreach { encoding =>
    it should s"load all contigs by name from a sidecarless $encoding reference" in {
      val dict = SequenceDictionary(SequenceMetadata("chr1", 5), SequenceMetadata("chr2", 5))
      val builder = new SamBuilder(readLength=4, sort=Some(SamOrder.Coordinate), sd=Some(dict))
      builder.addFrag(name="a", bases="CAAG", start=1, cigar="3M1D1M", contig=0)
      builder.addFrag(name="b", bases="GTTC", start=1, cigar="3M1D1M", contig=1)
      builder.addFrag(name="c", bases="CAAG", start=1, cigar="3M1D1M", contig=0)
      val plain = fasta(Seq("extra" -> "ACGT", "chr2" -> "gtttc", "chr1" -> "caaag"))
      // Names must be truncated at whitespace, independently of dictionary/FASTA order.
      Files.write(plain, new String(Files.readAllBytes(plain)).replace(">chr1\n", ">chr1 description\n").getBytes)
      val ref = encoding match {
        case "plain" => plain
        case other   => compressedFasta(plain, bgzf=other == "BGZF")
      }
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref), useEqualsAndX=true))
      val records = readBamRecs(output)
      records.map(r => (r.name, r.refName, r.basesString, r.cigar.toString)).sortBy(_._1) shouldBe Seq(
        ("a", "chr1", "CAAG", "1=1D3="), ("b", "chr2", "GTTC", "1=1D3="), ("c", "chr1", "CAAG", "1=1D3=")
      )
    }
  }

  it should "reject a duplicated FASTA contig that is not in the input's sequence dictionary" in {
    val builder = new SamBuilder(readLength=4, sd=Some(SequenceDictionary(SequenceMetadata("chr1", 4))))
    builder.addFrag(name="unchangeable", start=1, bases="ACGT")
    val ref = makeTempFile("test.", ".fa")
    Files.write(ref, ">chr1\nACGT\n>extra\nAAAA\n>extra\nCCCC\n".getBytes)
    val output = makeTempFile("test.", ".bam")
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref)))
    }
    ex.message.value should include ("Contig extra occurs more than once in the reference FASTA")
  }

  // A malformed dictionary is ignored by htsjdk, so only a well-formed stale dictionary shows that none is used
  Seq("matching", "stale", "malformed", "stale dictionary only").foreach { sidecars =>
    it should s"read actual FASTA bases independently of $sidecars sidecars" in {
      val in = inputs(SamOrder.Coordinate, md=false)
      val referenceBuilder = new ReferenceSetBuilder()
      referenceBuilder.add("chr1").add(RefBases)
      val ref = referenceBuilder.toTempFile()
      val fai = ReferenceSequenceFileFactory.getFastaIndexFileName(ref)
      val dict = PathUtil.replaceExtension(ref, ".dict")
      if (sidecars != "matching") {
        val malformed = sidecars.startsWith("malformed")
        Files.write(dict, (if (malformed) "not a dictionary\n" else StaleDictionary).getBytes)
        if (sidecars.endsWith("only")) Files.delete(fai)
        else Files.write(fai, (if (malformed) "not an index\n" else "chr1\t1\t999999\t80\t81\n").getBytes)
      }
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
      readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
      readBamRecs(output).find(_.name == "frag4").value.cigar.toString shouldBe "30M1D70M"
    }
  }

  Seq(false, true).foreach { bgzf =>
    it should s"ignore malformed index and stale dictionary sidecars for compressed references with BGZF=$bgzf" in {
      val in = inputs(SamOrder.Queryname, md=false)
      val ref = compressedFasta(in.ref, bgzf=bgzf)
      Seq(".fai", ".gzi").foreach { suffix =>
        val sidecar = ref.resolveSibling(ref.getFileName.toString + suffix)
        Files.write(sidecar, "malformed sidecar".getBytes)
        sidecar.toFile.deleteOnExit()
      }
      // Both places that htsjdk looks for a dictionary: e.g. `ref.dict` and `ref.fa.gz.dict`
      Seq(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(ref),
          ref.resolveSibling(ref.getFileName.toString + ".dict")).foreach { dict =>
        Files.write(dict, StaleDictionary.getBytes)
        dict.toFile.deleteOnExit()
      }
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
      val records = readBamRecs(output)
      records.find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
      records.find(_.name == "frag4").value.cigar.toString shouldBe "30M1D70M"
    }
  }

  Seq(true, false).foreach { empty =>
    Seq("missing", "wrong length").foreach { problem =>
      it should s"validate an unused $problem contig with empty=$empty" in {
        val dict = SequenceDictionary(SequenceMetadata("chr1", 10), SequenceMetadata("unused", 10))
        val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
        if (!empty) builder.addFrag(name="unchangeable", start=1, bases="A" * 10)
        val contigs = Seq("chr1" -> ("A" * 10)) ++ (if (problem == "missing") Seq.empty else Seq("unused" -> "AAAA"))
        val output = makeTempFile("test.", ".bam")
        val ex = intercept[FailureException] {
          executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(fasta(contigs))))
        }
        ex.message.value should include (if (problem == "missing") "Contig unused is not in the reference FASTA" else "Contig unused has length")
      }
    }
    Seq("chr1", "chr1 another description").foreach { duplicate =>
      it should s"reject duplicate FASTA names after whitespace truncation with empty=$empty and second name '$duplicate'" in {
        val builder = new SamBuilder(readLength=4, sd=Some(SequenceDictionary(SequenceMetadata("chr1", 4))))
        if (!empty) builder.addFrag(name="unchangeable", start=1, bases="ACGT")
        val ref = makeTempFile("test.", ".fa")
        Files.write(ref, s">chr1 first description\nACGT\n>$duplicate\nTGCA\n".getBytes)
        val output = makeTempFile("test.", ".bam")
        val ex = intercept[FailureException] {
          executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(ref)))
        }
        ex.message.value should include ("Contig chr1 occurs more than once in the reference FASTA")
      }
    }
  }

  private val PreservedHeaders = Seq(
    (SortOrder.queryname, GroupOrder.none, None),
    (SortOrder.queryname, GroupOrder.none, Some("queryname:natural")),
    (SortOrder.queryname, GroupOrder.query, Some("queryname:lexicographical")),
    (SortOrder.queryname, GroupOrder.query, Some("unrecognized")),
    (SortOrder.unsorted, GroupOrder.query, None),
    (SortOrder.unsorted, GroupOrder.query, Some("random-query")),
    (SortOrder.unsorted, GroupOrder.query, Some("unsorted:random-query")),
    (SortOrder.unsorted, GroupOrder.query, Some("unsorted:unrecognized"))
  )

  PreservedHeaders.foreach { case (so, go, ss) =>
    it should s"preserve group order, within-group order, and the header for SO=$so GO=$go SS=$ss" in {
      val builder = new SamBuilder(readLength=10)
      val read2 = builder.addPair(name="read2", start1=101, start2=201)
      val read10 = builder.addPair(name="read10", start1=11, start2=51)
      val records = Seq(read2(1), read2(0), read10(1), read10(0))
      val header = builder.header.clone()
      header.setSortOrder(so)
      header.setGroupOrder(go)
      header.setAttribute("SS", ss.orNull)
      header.addComment("Preserve unrelated metadata")
      val input = writeRecords(header, records)
      val originalHeader = Using.resource(SamSource(input))(_.header.getSAMString)
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=LeadingIndelMode.Trim))
      identities(readBamRecs(output)) shouldBe identities(records)
      Using.resource(SamSource(output))(_.header.getSAMString) shouldBe originalHeader
    }
  }

  private val RestoredHeaders = Seq(
    (SortOrder.coordinate, GroupOrder.none, None, SamOrder.Coordinate),
    (SortOrder.coordinate, GroupOrder.query, Some("coordinate:unrecognized"), SamOrder.Coordinate),
    (SortOrder.unsorted, GroupOrder.query, Some("template-coordinate"), SamOrder.TemplateCoordinate),
    (SortOrder.unsorted, GroupOrder.query, Some("unsorted:template-coordinate"), SamOrder.TemplateCoordinate),
    (SortOrder.unsorted, GroupOrder.query, Some("unsorted:template-coordinate:extra"), SamOrder.TemplateCoordinate),
    (SortOrder.unsorted, GroupOrder.query, Some("template-coordinate:extra"), SamOrder.TemplateCoordinate),
    (SortOrder.unsorted, GroupOrder.none, Some("random"), SamOrder.Random),
    (SortOrder.unsorted, GroupOrder.none, Some("unsorted:random"), SamOrder.Random),
    (SortOrder.unsorted, GroupOrder.none, Some("unsorted:random:extra"), SamOrder.Random),
    (SortOrder.unsorted, GroupOrder.none, Some("random:extra"), SamOrder.Random),
    (SortOrder.unsorted, GroupOrder.none, None, SamOrder.Queryname),
    (SortOrder.unknown, GroupOrder.query, None, SamOrder.Queryname),
    (SortOrder.unsorted, GroupOrder.none, Some("unrecognized"), SamOrder.Queryname)
  )

  RestoredHeaders.foreach { case (so, go, ss, order) =>
    it should s"group complete templates and restore $order for SO=$so GO=$go SS=$ss" in {
      val builder = new SamBuilder(readLength=10)
      val z = builder.addPair(name="z", start1=11, start2=201)
      val a = builder.addPair(name="a", start1=101, start2=301)
      val header = builder.header.clone()
      header.setSortOrder(so)
      header.setGroupOrder(go)
      header.setAttribute("SS", ss.orNull)
      header.addComment("Retain this comment when normalizing order")
      // Missing and stale mate CIGARs show that mate information is repaired
      z(0).remove("MC")
      a(1)("MC") = "5M"
      // Even GO:query does not make coordinate/template-coordinate input safe for passthrough.
      val input = writeRecords(header, Seq(z(0), a(0), z(1), a(1)))
      val originalHeader = header.getSAMString
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=LeadingIndelMode.Keep))
      val records = readBamRecs(output)
      records.size shouldBe 4
      assertSorted(records, order)
      records.foreach(r => r.mateCigar.value.toString shouldBe "10M")
      Using.resource(SamSource(output)) { source =>
        val expected = order.applyTo(header.clone())
        source.header.getSAMString shouldBe expected.getSAMString
      }
      header.getSAMString shouldBe originalHeader
    }
  }

  it should "treat GO:query input without SO as being of unknown order, so group it and write it in query name order" in {
    val builder = new SamBuilder(readLength=10)
    val z = builder.addPair(name="z", start1=11, start2=201)
    val a = builder.addPair(name="a", start1=101, start2=301)
    val header = builder.header.clone()
    header.setGroupOrder(GroupOrder.query)
    // htsjdk always writes SO, so write the SAM text without it, and with the templates' records not adjacent
    val headerText = header.getSAMString.replaceAll("\tSO:[^\t\n]*", "")
    headerText should not include ("SO:")
    val input = makeTempFile("test.", ".sam")
    Files.write(input, (headerText + Seq(z(0), a(0), z(1), a(1)).map(_.asSam.getSAMString.trim).mkString("", "\n", "\n")).getBytes)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output))
    val records = readBamRecs(output)
    records.map(_.id) shouldBe Seq("a/1", "a/2", "z/1", "z/2")
    Using.resource(SamSource(output))(source => SamOrder(source.header)) shouldBe Some(SamOrder.Queryname)
  }

  Seq("template-coordinate", "unsorted:template-coordinate", "unsorted:template-coordinate:extra").foreach { ss =>
    LeadingIndelMode.values.foreach { mode =>
      it should s"repair interleaved primary, supplementary, and secondary records with SS=$ss in $mode" in {
        val in = leadingIndelInputs()
        val primary = in.sam.iterator.toIndexedSeq
        val secondary = in.sam.addPair(name="pair2", start1=201, start2=551)(1)
        secondary.secondary = true
        secondary("MC") = "4M2D96M"
        val bridgePair = in.sam.addPair(name="bridge-pair", start1=201, start2=451)
        val bridgeFrag = in.sam.addFrag(name="bridge-frag", start=301).value
        val header = SamOrder.TemplateCoordinate.applyTo(in.sam.header.clone())
        header.setAttribute("SS", ss)
        val records = (primary ++ bridgePair ++ Seq(secondary, bridgeFrag)).sortBy(SamOrder.TemplateCoordinate.sortkey)
        Seq("pair2", "frag5").foreach { name =>
          val indices = records.indices.filter(i => records(i).name == name)
          indices.last - indices.head should be > (indices.size - 1)
        }
        val input = writeRecords(header, records)
        val output = makeTempFile("test.", ".bam")
        executeFgbioTool(new LeftAlignIndels(input=input, output=output, ref=Some(in.ref), leadingIndelMode=mode))
        val result = readBamRecs(output)
        assertSorted(result, SamOrder.TemplateCoordinate)
        Using.resource(SamSource(output))(_.header.getAttribute("SS")) shouldBe "unsorted:template-coordinate"
        result.find(r => r.name == "frag5" && !r.supplementary).value[String]("SA") shouldBe "chr1,181,+,50S20M2D30M,60,2;"
        // Only pair2 has a leading indel, so with Filter only its records are removed, and otherwise none are
        val expected = records.filterNot(r => mode == LeadingIndelMode.Filter && r.name == "pair2")
        result.map(r => (r.name, r.flags)) should contain theSameElementsAs expected.map(r => (r.name, r.flags))
        if (mode != LeadingIndelMode.Filter) {
          val sec = result.find(_.secondary).value
          sec.mateStart shouldBe (if (mode == LeadingIndelMode.Trim) 203 else 201)
          sec[String]("MC") shouldBe (if (mode == LeadingIndelMode.Trim) "100M" else "2D100M")
        }
      }
    }
  }

  Seq(LeadingIndelMode.Keep, LeadingIndelMode.Filter).foreach { mode =>
    it should s"restore template-coordinate order from repaired stale MC keys without trimming in $mode" in {
      val builder = new SamBuilder(readLength=10)
      val z = builder.addPair(name="z", start1=101, start2=301, strand1=Minus, strand2=Plus)
      val a = builder.addPair(name="a", start1=101, start2=301, strand1=Minus, strand2=Plus)
      // A stale short mate CIGAR puts z's R2 first; repair makes the templates tie until their names.
      z(1)("MC") = "5M"
      val header = SamOrder.TemplateCoordinate.applyTo(builder.header.clone())
      val input = writeRecords(header, (z ++ a).sortBy(SamOrder.TemplateCoordinate.sortkey))
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=mode))
      val records = readBamRecs(output)
      records.foreach(r => r[String]("MC") shouldBe "10M")
      assertSorted(records, SamOrder.TemplateCoordinate)
      records.map(_.name) shouldBe Seq("a", "a", "z", "z")
    }
  }

  /** A retained record whose mapped mate is unavailable to repair its missing MC tag. */
  private def missingMcRecord(source: String): SamRecord = {
    val builder = new SamBuilder(readLength=10)
    val rec = builder.addPair(name=if (source == "unavailable name") "*" else "missing-mc", start1=11, start2=101,
      unmapped1=source == "unmapped read with a mapped mate")(0)
    rec.secondary = source == "secondary"
    rec.remove("MC")
    rec
  }

  Seq("secondary", "unavailable name", "allowed missing mate", "unmapped read with a mapped mate").foreach { source =>
    it should s"fail contextually for retained missing MC from $source in template-coordinate output" in {
      val rec = missingMcRecord(source)
      val input = writeRecords(SamOrder.TemplateCoordinate.applyTo(rec.header.clone()), Seq(rec))
      val output = makeTempFile("test.", ".bam")
      val ex = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true, leadingIndelMode=LeadingIndelMode.Keep))
      }
      ex.message.value should include (rec.name)
      ex.message.value should include ("MC")
      ex.message.value should include ("template-coordinate")
      ex.message.value should include ("samtools fixmate -m")
      ex.message.value should include ("SetMateInformation")
      ex.message.value should include ("query-name")
    }
    it should s"allow missing MC from $source when output is query-name sorted" in {
      val rec = missingMcRecord(source)
      val input = writeRecords(SamOrder.Queryname.applyTo(rec.header.clone()), Seq(rec))
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true))
      readBamRecs(output).size shouldBe 1
    }
  }

  it should "reject an unmatched secondary without MC even when both primary mates are present" in {
    val builder = new SamBuilder(readLength=10)
    val pair = builder.addPair(name="unmatched", start1=11, start2=101)
    val secondary = builder.addPair(name="unmatched", start1=201, start2=301).head
    secondary.secondary = true
    secondary.remove("MC")
    // Its mate start refers to a different alignment, so repair must not assign the primary mate's CIGAR.
    secondary.mateStart should not be pair(1).start
    val input = writeRecords(SamOrder.TemplateCoordinate.applyTo(builder.header.clone()), pair :+ secondary)
    val output = makeTempFile("test.", ".bam")
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=LeadingIndelMode.Keep))
    }
    ex.message.value should include ("unmatched")
    ex.message.value should include ("MC")
    ex.message.value should include ("template-coordinate")
  }

  it should "allow unpaired and mate-unmapped records without MC in template-coordinate output" in {
    val builder = new SamBuilder(readLength=10)
    val unpaired = builder.addFrag(name="*", start=11).value
    val mateUnmapped = builder.addPair(name="mate-unmapped", start1=21, unmapped2=true)(0)
    mateUnmapped.remove("MC")
    val input = writeRecords(SamOrder.TemplateCoordinate.applyTo(builder.header.clone()), Seq(unpaired, mateUnmapped))
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true, leadingIndelMode=LeadingIndelMode.Keep))
    val records = readBamRecs(output)
    records.size shouldBe 2
    records.foreach(_.contains("MC") shouldBe false)
    assertSorted(records, SamOrder.TemplateCoordinate)
  }

  it should "repair missing MC before checking template-coordinate output" in {
    val builder = new SamBuilder(readLength=10)
    val pair = builder.addPair(name="repaired", start1=11, start2=101)
    pair.foreach(_.remove("MC"))
    val input = writeRecords(SamOrder.TemplateCoordinate.applyTo(builder.header.clone()), pair)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=LeadingIndelMode.Keep))
    val records = readBamRecs(output)
    records.map(_.id) should contain theSameElementsAs pair.map(_.id)
    records.foreach(r => r[String]("MC") shouldBe "10M")
  }
  Seq("secondary", "unavailable name", "allowed missing mate").foreach { source =>
    it should s"filter a shifted leading indel before checking missing MC from $source" in {
      val rec = missingMcRecord(source)
      rec.bases = "A" * 10
      rec.cigar = "5M1D5M"
      rec("MD") = "5^A5"
      val input = writeRecords(SamOrder.TemplateCoordinate.applyTo(rec.header.clone()), Seq(rec))
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true, leadingIndelMode=LeadingIndelMode.Filter))
      readBamRecs(output) shouldBe empty
    }
  }

  Seq("template", "unavailable name", "duplicate unpaired").foreach { branch =>
    it should s"align before filtering indels shifted to the start in the $branch branch" in {
      val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname))
      val name = if (branch == "unavailable name") "*" else "reads"
      val changed = if (branch == "template") {
        builder.addPair(name=name, bases1="A" * 10, start1=11, start2=101, cigar1="5M1D5M")(0)
      }
      else builder.addFrag(name=name, bases="A" * 10, start=11, cigar="5M1D5M").value
      changed("MD") = "5^A5"
      if (branch != "template") builder.addFrag(name=name, start=201, attrs=Map("SA" -> "chr1,11,+,5M1D5M,60,1;"))
      val output = makeTempFile("test.", ".bam")
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, leadingIndelMode=LeadingIndelMode.Filter))
      val records = readBamRecs(output)
      if (branch == "template") records shouldBe empty
      else {
        records.map(_.start) shouldBe Seq(201)
        records.head[String]("SA") shouldBe "chr1,11,+,5M1D5M,60,1;"
      }
    }
  }

  it should "align every record of a valid template before filtering the template" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname))
    val pair = builder.addPair(name="align-all", bases1="A" * 10, start1=11, start2=101, cigar1="1D10M", cigar2="5M1D5M")
    pair.head("MD") = "0^A10"
    val output = makeTempFile("test.", ".bam")
    val ex = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, leadingIndelMode=LeadingIndelMode.Filter))
    }
    ex.message.value should include ("has no MD tag")
  }

  Seq(false, true).foreach { secondary =>
    it should s"report a missing mate before missing MD, including before filtering, with secondary=$secondary" in {
      val builder = new SamBuilder(readLength=10)
      val rec = builder.addPair(name="precedence", start1=11, start2=101, cigar1="5M1D5M")(0)
      rec.secondary = secondary
      val input = writeRecords(SamOrder.Queryname.applyTo(builder.header.clone()), Seq(rec))
      val output = makeTempFile("test.", ".bam")
      val ex = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=input, output=output, leadingIndelMode=LeadingIndelMode.Filter))
      }
      ex.message.value should include ("is paired but its mate was not found")
      val allowed = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true, leadingIndelMode=LeadingIndelMode.Filter))
      }
      allowed.message.value should include ("has no MD tag")
    }
  }

}
