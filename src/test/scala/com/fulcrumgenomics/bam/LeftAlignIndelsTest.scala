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
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import htsjdk.samtools.util.SequenceUtil
import org.scalatest.OptionValues

import java.nio.file.Files

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
                    expectedStart: Int = 100): SamRecord = {
    val builder = new SamBuilder(readLength=bases.length)
    val rec     = builder.addFrag(bases=bases, cigar=cigar, start=100, attrs=attrs).value
    withClue(f"Hint: bases: $bases cigar: $cigar ::") {
      aligner(useEqualsAndX=useEqualsAndX, trim=trim).align(rec, ref).changed shouldBe changed
      rec.cigar.toString shouldBe expected
      rec.start shouldBe expectedStart
    }
    rec
  }

  /** Fixes the template with the given method, returning the number of SA tag entries that were not updated. */
  private def unmatched(fix: (Template, String => Unit) => Unit, template: Template): Int = {
    var count = 0
    fix(template, _ => count += 1)
    count
  }

  /** Builds an aligner that does not trim leading indels unless `trim` is true. */
  private def aligner(useEqualsAndX: Boolean = false, trim: Boolean = false): IndelLeftAligner =
    new IndelLeftAligner(useEqualsAndX=useEqualsAndX, trimLeadingIndels=trim)

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

  "IndelLeftAligner.align" should "convert M to = and X when useEqualsAndX is true" in {
    align("A"*10, "10M", "A"*10, changed=true, expected="10=", useEqualsAndX=true)
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
    aligner(trim=true).align(rec, "CACG") shouldBe LeftAlignResult(changed=true, trimmed=true)
    rec.cigar.toString shouldBe "1S3M"
    rec.start shouldBe 101
    rec[String]("MD") shouldBe "3"
    rec[Int]("NM") shouldBe 0
    rec[Int]("UQ") shouldBe 0
  }

  it should "not trim a leading indel when no aligned bases follow it in the first block" in {
    align("AAAA", "2I10N2M", "N"*10 + "AA", changed=false, expected="2I10N2M", trim=true)
  }

  it should "count a record as trimmed only when an indel was trimmed" in {
    val builder = new SamBuilder(readLength=4)
    val rec     = builder.addFrag(bases="CAAG", cigar="3M1D1M", start=100).value
    aligner(trim=true).align(rec, "CAAAG") shouldBe LeftAlignResult(changed=true, trimmed=false)
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
    val rec = align("TTTTT" + "CAAG", "5M100N3M1D1M", "TTTTT" + "N"*100 + "CAAAG", changed=true, expected="5M100N1M1D3M",
      attrs=Map("MD" -> "8^A1", "NM" -> 1))
    rec[String]("MD") shouldBe "6^A3"
    rec[Int]("NM") shouldBe 1
  }

  it should "not shift an indel across a reference skip" in {
    align("AAAA" + "AAAA", "4M10N1M1I2M", "AAAA" + "N"*10 + "AAA", changed=true, expected="4M10N1I3M")
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

  it should "count an N in the read or reference as a mismatch in the CIGAR, MD, and NM, as samtools does" in {
    val rec = align("NAC", "2M1D1M", "NAAC", changed=true, expected="1X1D2=", useEqualsAndX=true,
      attrs=Map("MD" -> "0N1^A1", "NM" -> 2))
    rec[String]("MD") shouldBe "0N0^A2"
    rec[Int]("NM") shouldBe 2
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
    sec("MC") = "10M"
    r1("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    supp("SA") = s"${r1.refName},100,+,10M,60,0;"
    // pretend r2 and supp were left-aligned
    r2.cigar   = "2I8M"
    supp.cigar = "5S5="

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
    // pretend supp was left-aligned
    supp.cigar = "5S1I4M"

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
    supp.cigar = "5H1I4M" // pretend supp was left-aligned

    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, supp)))
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S1I4M,60,1;"
  }

  it should "not update an SA entry that does not match exactly one alignment" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100, cigar="10M").value
    val supp1   = builder.addFrag(name="q1", start=1000, cigar="5S1I4M").value
    val supp2   = builder.addFrag(name="q1", start=1000, cigar="5S1I4M").value
    Seq(supp1, supp2).foreach(_.supplementary = true)
    // two alignments match the first entry, and no alignment matches the CIGAR of the second entry
    r1("SA") = s"${supp1.refName},1000,+,5S2M1I2M,60,1;${supp1.refName},1000,+,10M,60,0;"

    unmatched(IndelLeftAligner.fixMateAndSupplementaryInfo(_, _), Template(Iterator(r1, supp1, supp2))) shouldBe 2
    r1[String]("SA") shouldBe s"${supp1.refName},1000,+,5S2M1I2M,60,1;${supp1.refName},1000,+,10M,60,0;"
  }

  it should "not update an SA entry from an alignment of the other read of the pair" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(name="q1", start1=100, start2=300)
    val supp = builder.addFrag(name="q1", start=1000, cigar="5S5M").value
    supp.paired = true; supp.secondOfPair = true; supp.supplementary = true
    r1("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    supp.cigar = "5S5=" // pretend the R2 supplementary was left-aligned
    unmatched(IndelLeftAligner.fixMateAndSupplementaryInfo(_, _), Template(Iterator(r1, r2, supp))) shouldBe 1
    r1[String]("SA") shouldBe s"${supp.refName},1000,+,5S5M,60,0;"
  }

  it should "not update an SA entry from an alignment on the other strand" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100).value
    val supp    = builder.addFrag(name="q1", start=1000, cigar="5S5M", strand=Plus).value
    supp.supplementary = true
    r1("SA") = s"${supp.refName},1000,-,5S5M,60,0;"
    supp.cigar = "5S5="
    unmatched(IndelLeftAligner.fixMateAndSupplementaryInfo(_, _), Template(Iterator(r1, supp))) shouldBe 1
    r1[String]("SA") shouldBe s"${supp.refName},1000,-,5S5M,60,0;"
  }

  it should "not update an SA entry from the record itself" in {
    val builder = new SamBuilder(readLength=10)
    val r1      = builder.addFrag(name="q1", start=100).value
    val supp    = builder.addFrag(name="q1", start=1000, cigar="5S5M").value
    supp.supplementary = true
    supp("SA") = s"${supp.refName},1000,+,5S5M,60,0;"
    supp.cigar = "5S5="
    unmatched(IndelLeftAligner.fixMateAndSupplementaryInfo(_, _), Template(Iterator(r1, supp))) shouldBe 1
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
    supp1.start = 1002; supp1.cigar = "5S5M"
    supp2.cigar = "6S4M"
    unmatched(IndelLeftAligner.fixMateAndSupplementaryInfo(_, _), Template(Iterator(r1, supp1, supp2))) shouldBe 0
    r1[String]("SA") shouldBe s"${supp1.refName},1002,+,5S5M,60,2;${supp2.refName},2000,+,6S4M,60,1;"
  }

  it should "update the mate start and cigar of a secondary record whose mate's leading deletion was trimmed" in {
    val builder = new SamBuilder(readLength=10)
    val Seq(r1, r2) = builder.addPair(name="q1", start1=100, start2=300, cigar1="10M", cigar2="2D10M")
    def secondary(mateStart: Int, mateCigar: String): SamRecord = {
      val sec = builder.addFrag(name="q1", start=2000, cigar="10M").value
      sec.paired = true; sec.firstOfPair = true; sec.secondary = true
      sec.mateRefIndex = r2.refIndex; sec.mateStart = mateStart; sec.mateNegativeStrand = r2.negativeStrand
      sec("MC") = mateCigar
      sec
    }
    val sec   = secondary(mateStart=300, mateCigar="2D10M")
    val decoy = secondary(mateStart=310, mateCigar="10M") // neither its mate start nor its mate end match R2
    // pretend R2's leading deletion was trimmed
    r2.start = 302; r2.cigar = "10M"
    IndelLeftAligner.fixMateAndSupplementaryInfo(Template(Iterator(r1, r2, sec, decoy)))
    sec.mateStart shouldBe 302
    sec.insertSize shouldBe htsjdk.samtools.SamPairUtil.computeInsertSize(sec.asSam, r2.asSam)
    sec.insertSize should not be 0
    sec[String]("MC") shouldBe "10M"
    decoy.mateStart shouldBe 310
    decoy[String]("MC") shouldBe "10M"
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

  private case class Inputs(sam: SamBuilder, ref: java.nio.file.Path)

  /** Builds a reference FASTA and reads: a fragment with no indels, a fragment with a deletion in the repeat, a
    * fragment with a deletion in the soft-masked homopolymer, a pair with an insertion in the repeat, and an unmapped
    * fragment. */
  private def inputs(samOrder: SamOrder, md: Boolean = true, indexed: Boolean = true): Inputs = {
    val referenceBuilder = new ReferenceSetBuilder()
    referenceBuilder.add("chr1").add(RefBases)
    val ref = referenceBuilder.toTempFile()
    if (!indexed) Files.delete(ref.resolveSibling(ref.getFileName.toString + ".fai"))

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

  it should "use an indexed reference with a query name sorted input without MD tags" in {
    val in     = inputs(SamOrder.Queryname, md=false)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
  }

  it should "use a reference without an index with a coordinate sorted input without MD tags" in {
    val in     = inputs(SamOrder.Coordinate, md=false, indexed=false)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    readBamRecs(output).find(_.name == "frag2").value.cigar.toString shouldBe "50M2D50M"
  }

  it should "fail with a helpful message when the MD tag is needed but missing" in {
    val in     = inputs(SamOrder.Queryname, md=false)
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output))
    }
    ex.message.value should include ("has no MD tag")
  }

  it should "fail when given a reference FASTA without an index for an input that is not coordinate sorted" in {
    val in     = inputs(SamOrder.Queryname, indexed=false)
    val output = makeTempFile("test.", ".bam")
    val ex     = intercept[FailureException] {
      executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(in.ref)))
    }
    ex.message.value should include ("samtools faidx")
  }

  /** Writes a reference FASTA with the given contigs, with or without an index. */
  private def fasta(contigs: Seq[(String, String)], indexed: Boolean): java.nio.file.Path = {
    val referenceBuilder = new ReferenceSetBuilder()
    contigs.foreach { case (name, bases) => referenceBuilder.add(name).add(bases) }
    val ref = referenceBuilder.toTempFile()
    if (!indexed) Files.delete(ref.resolveSibling(ref.getFileName.toString + ".fai"))
    ref
  }

  Seq(true, false).foreach { indexed =>
    val (samOrder, desc) = if (indexed) (SamOrder.Queryname, "an indexed") else (SamOrder.Coordinate, "a non-indexed")
    it should f"fail when $desc reference FASTA is missing a contig in the BAM's sequence dictionary" in {
      val in     = inputs(samOrder)
      val output = makeTempFile("test.", ".bam")
      val ref    = fasta(Seq("chr2" -> RefBases), indexed=indexed)
      val ex     = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
      }
      ex.message.value should include (if (indexed) "Contig chr1 is not in the reference FASTA" else "Could not retrieve contig chr1")
    }

    it should f"fail when $desc reference FASTA has a contig length that differs from the BAM's sequence dictionary" in {
      val in     = inputs(samOrder)
      val output = makeTempFile("test.", ".bam")
      val ref    = fasta(Seq("chr1" -> (RefBases + "ACGT")), indexed=indexed)
      val ex     = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=in.sam.toTempFile(), output=output, ref=Some(ref)))
      }
      ex.message.value should include (if (indexed) "Contig chr1 has length" else "Could not retrieve contig chr1")
    }

    it should f"fail with a helpful message when a read extends past the end of its contig with $desc reference FASTA" in {
      val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
      val builder = new SamBuilder(readLength=100, sort=Some(samOrder), sd=Some(dict))
      builder.addFrag(name="frag1", bases="A"*100, start=RefBases.length - 50, cigar="50M1D50M")
      val output = makeTempFile("test.", ".bam")
      val ex     = intercept[FailureException] {
        executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output, ref=Some(fasta(Seq("chr1" -> RefBases), indexed=indexed))))
      }
      ex.message.value should include ("past the end of contig chr1")
    }
  }

  it should "use an indexed reference FASTA whose contigs are in a different order than the BAM's sequence dictionary" in {
    val in     = inputs(SamOrder.Coordinate, md=false)
    val output = makeTempFile("test.", ".bam")
    val ref    = fasta(Seq("chr0" -> "ACGT", "chr1" -> RefBases), indexed=true)
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
    val input  = makeTempFile("test.", ".bam")
    val writer = SamWriter(input, header)
    writer ++= in.sam.toIndexedSeq.filterNot(r => r.name.startsWith("read")).sortBy(SamOrder.Queryname.sortkey)
    writer ++= Seq(read2, read10)
    writer.close()
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
    val input  = makeTempFile("test.", ".bam")
    val writer = SamWriter(input, header)
    writer ++= Seq(r2, r1)
    writer.close()

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

  it should "fail on a paired read without its mate, unless missing mates are allowed" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Queryname), sd=Some(dict))
    val Seq(r1, _) = builder.addPair(name="pair1", start1=11, start2=101)
    val input  = makeTempFile("test.", ".bam")
    val writer = SamWriter(input, builder.header)
    writer += r1
    writer.close()

    val output = makeTempFile("test.", ".bam")
    val ex = intercept[FailureException] { executeFgbioTool(new LeftAlignIndels(input=input, output=output)) }
    ex.message.value should include ("is paired but its mate was not found")
    // the partially written output is removed
    Files.exists(output) shouldBe false
    executeFgbioTool(new LeftAlignIndels(input=input, output=output, allowMissingMates=true))
    readBamRecs(output).map(_.name) shouldBe Seq("pair1")
  }

  /** Builds a reference FASTA and coordinate sorted reads: a pair whose R1 has a deletion in the repeat that is shifted
    * to the start of its alignment, and a fragment with a supplementary alignment whose deletion is shifted. */
  private def leadingIndelInputs(samOrder: SamOrder = SamOrder.Coordinate): Inputs = {
    val ref     = fasta(Seq("chr1" -> RefBases), indexed=true)
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
    val header = SamOrder.TemplateCoordinate.applyTo(in.sam.header.clone())
    val input  = makeTempFile("test.", ".bam")
    val writer = SamWriter(input, header)
    writer ++= in.sam.toIndexedSeq.sortBy(SamOrder.TemplateCoordinate.sortkey)
    writer.close()
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=input, output=output, ref=Some(in.ref)))
    val records = readBamRecs(output)
    records.map(_.name) shouldBe records.sortBy(SamOrder.TemplateCoordinate.sortkey).map(_.name)
    records.find(r => r.name == "pair2" && r.firstOfPair).value.start shouldBe 203
  }

  it should "write unpaired records that share a query name without updating their SA tags" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=RefBases.length))
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.Coordinate), sd=Some(dict))
    builder.addFrag(name="dup", start=11)
    builder.addFrag(name="dup", start=101)
    val output = makeTempFile("test.", ".bam")
    executeFgbioTool(new LeftAlignIndels(input=builder.toTempFile(), output=output))
    readBamRecs(output).map(_.name) shouldBe Seq("dup", "dup")
  }

  "LeftAlignIndels.deleteOutput" should "delete the output and its index, but not a file that is not a regular file" in {
    val output = makeTempFile("test.", ".bam")
    val name   = output.getFileName.toString
    val bai    = output.resolveSibling(name.dropRight(4) + ".bai")
    val bamBai = output.resolveSibling(name + ".bai")
    Seq(bai, bamBai).foreach { p => Files.write(p, Array[Byte](1)); p.toFile.deleteOnExit() }
    LeftAlignIndels.deleteOutput(output)
    Seq(output, bai, bamBai).foreach(p => Files.exists(p) shouldBe false)
    // a missing output and a device are left alone without failing
    LeftAlignIndels.deleteOutput(output)
    LeftAlignIndels.deleteOutput(java.nio.file.Paths.get("/dev/null"))
    Files.exists(java.nio.file.Paths.get("/dev/null")) shouldBe true
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
}
