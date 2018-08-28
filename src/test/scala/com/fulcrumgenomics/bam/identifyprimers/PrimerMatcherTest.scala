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

import com.fulcrumgenomics.alignment.{Aligner, Alignment, Mode}
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.util.{Interval, SequenceUtil}
import org.scalatest.OptionValues

object PrimerMatcherTest {
  val primers = IndexedSeq(
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

  def fragFromPrimer(primer: Primer, builder: SamBuilder = new SamBuilder(), startAdjust: Int = 0): SamRecord = {
    require(primer.length <= builder.readLength)
    val contig = builder.dict.getSequenceIndex(primer.ref_name)
    val start  = if (primer.positiveStrand) primer.start else primer.end - builder.readLength + 1
    val bases  = {
      val remaining = builder.readLength - primer.length
      val primerBases = new String(primer.bases)
      if (primer.positiveStrand) primerBases + ("N" * remaining)
      else ("N" * remaining) + primerBases
    }
    val strand = if (primer.positiveStrand) SamBuilder.Plus else SamBuilder.Minus
    builder.addFrag(contig=contig, start=start+startAdjust, bases=bases, unmapped = false, strand=strand).get
  }
}

/** Tests methods in [[PrimerMatcherWithKmerFilterTest]]. */
class PrimerMatcherWithKmerFilterTest extends UnitSpec with OptionValues {
  private val kmerLength: Int = 5

  private val matcher = new UngappedAlignmentBasedPrimerMatcher(
    primers           = PrimerMatcherTest.primers,
    slop              = 5,
    ignorePrimerStrand = true,
    maxMismatchRate   = 0.1,
    kmerLength        = Some(kmerLength)
  )

  private val matcherNoCache = new UngappedAlignmentBasedPrimerMatcher(
    primers           = PrimerMatcherTest.primers,
    slop              = 5,
    ignorePrimerStrand = true,
    maxMismatchRate   = 0.1,
    kmerLength        = None
  )

  private implicit class SeqPrimerToUniqueSortedPrimerId(primers: Seq[Primer]) {
    /** Maps each primer to their primer id, then returns them sorted and distinct. */
    def uniqSortedId: Seq[String] = primers.map(_.primer_id).distinct.sorted
  }

  "PrimerMatcherWithKmerFilter.getPrimersFrom" should "return all primers that share a kmer with the given read" in {
    implicit def toBytes(s: String): Array[Byte] = s.getBytes

    def toSeqOfMatches(bases: String, startOffset: Int, endOffset: Int): Seq[Seq[Primer]] = {
      Seq(
        matcher.getPrimersFrom(bases, true, kmerLength),
        matcher.getPrimersFrom(bases, false, kmerLength),
        matcher.getPrimersFrom(bases, startOffset, endOffset, kmerLength)
      ).map(_.distinct)
    }

    an[Exception] should be thrownBy matcher.getPrimersFrom("AAAAA", 0, 1, kmerLength) // endOffset is too far
    an[Exception] should be thrownBy matcher.getPrimersFrom("AAAAA", 1, 1, kmerLength) // startOffset is too far
    matcher.getPrimersFrom("AAAAA", 1, 0, kmerLength) shouldBe 'empty

    toSeqOfMatches("AAAAA", 0, 0).foreach { matches =>
      matches.length shouldBe 3
      matches.foreach { p => p.sequence.contains("AAAAA") shouldBe true }
    }

    toSeqOfMatches("AAAAAAAAAA", 0, 5).foreach { matches =>
      matches.length shouldBe 3
      matches.foreach { p => p.sequence.contains("AAAAA") shouldBe true }
    }

    toSeqOfMatches("GATTA", 0, 0).foreach { matches =>
      matches.length shouldBe 2
      matches.foreach { p => p.sequence.contains("GATTA") shouldBe true }
    }

    {
      val read = matcher.primers.map(_.sequence).mkString("")
      val matches = matcher.getPrimersFrom(read, 0, read.length - kmerLength, kmerLength).distinct
      matches.length shouldBe matcher.primers.length
    }
  }

  "PrimerMatcherWithKmerFilter.getPrimersForAlignmentTasks" should "return alignment tasks for an unmapped read" in {
    val frag = new SamBuilder(readLength=20).addFrag(bases="A"*10+"G"*10, strand=SamBuilder.Plus, unmapped=true).value

    // AAAAA* matches (1+) [must be + and have AAAAA], while the reverse complement TTTTT* matches (1-) [must be - and have TTTTT].
    matcher.getPrimersForAlignmentTasks(frag, false).uniqSortedId should contain theSameElementsInOrderAs Seq("1+", "1-")
    // AAAAA* matches (1+, 2-, 3-) [must have AAAAA], while the reverse complement TTTTT* matches (1-, 2+, 3+) [must have TTTTT].
    matcher.getPrimersForAlignmentTasks(frag, true).uniqSortedId should contain theSameElementsInOrderAs Seq("1+", "1-", "2+", "2-", "3+", "3-")
    // all primers!
    matcherNoCache.getPrimersForAlignmentTasks(frag, false).uniqSortedId should contain theSameElementsInOrderAs matcher.primers.uniqSortedId
    // all primers, and all reverse complements!
    val matches = matcherNoCache.getPrimersForAlignmentTasks(frag, true)
    matches.uniqSortedId should contain theSameElementsInOrderAs matcher.primers.uniqSortedId
    matches.map(_.primer_id).sorted should contain theSameElementsInOrderAs (matcher.primers ++ matcher.primers).map(_.primer_id).sorted
  }

  it should "return alignment tasks for mapped reads" in {
    val fragPlus  = new SamBuilder(readLength=20).addFrag(bases="A"*20, strand=SamBuilder.Plus).value
    val fragMinus = new SamBuilder(readLength=20).addFrag(bases="A"*20, strand=SamBuilder.Minus).value

    // frag mapped to the + strand
    {
      // AAAAA* matches (1+) [must be + and have AAAAA]
      matcher.getPrimersForAlignmentTasks(fragPlus, false).uniqSortedId should contain theSameElementsInOrderAs Seq("1+")
      // AAAAA* matches (1+, 2-, 3-) [must have AAAAA], while the reverse complement TTTTT* matches (1-, 2+, 3+) [must have TTTTT].
      matcher.getPrimersForAlignmentTasks(fragPlus, true).uniqSortedId should contain theSameElementsInOrderAs Seq("1+", "1-", "2+", "2-", "3+", "3-")
      // all positive strand primers
      matcherNoCache.getPrimersForAlignmentTasks(fragPlus, false).uniqSortedId should contain theSameElementsInOrderAs matcher.primers.filter(_.positiveStrand).uniqSortedId
      // all primers, and all reverse complements!
      val matches = matcherNoCache.getPrimersForAlignmentTasks(fragPlus, true)
      matches.uniqSortedId should contain theSameElementsInOrderAs matcher.primers.uniqSortedId
      matches.map(_.primer_id).sorted should contain theSameElementsInOrderAs (matcher.primers ++ matcher.primers).map(_.primer_id).sorted
    }

    // frag mapped to the - strand
    {
      // AAAAA* matches (2-, 3-) [must be - and have AAAAA]
      matcher.getPrimersForAlignmentTasks(fragMinus, false).uniqSortedId should contain theSameElementsInOrderAs Seq("2-", "3-")
      // AAAAA* matches (1+, 2-, 3-) [must have AAAAA], while the reverse complement TTTTT* matches (1-, 2+, 3+) [must have TTTTT].
      matcher.getPrimersForAlignmentTasks(fragMinus, true).uniqSortedId should contain theSameElementsInOrderAs Seq("1+", "1-", "2+", "2-", "3+", "3-")
      // all negative strand primers
      matcherNoCache.getPrimersForAlignmentTasks(fragMinus, false).uniqSortedId should contain theSameElementsInOrderAs matcher.primers.filter(_.negativeStrand).uniqSortedId
      // all primers, and all reverse complements!
      val matches = matcherNoCache.getPrimersForAlignmentTasks(fragMinus, true)
      matches.uniqSortedId should contain theSameElementsInOrderAs matcher.primers.uniqSortedId
      matches.map(_.primer_id).sorted should contain theSameElementsInOrderAs (matcher.primers ++ matcher.primers).map(_.primer_id).sorted
    }
  }
}

class UngappedAlignmentTest extends UnitSpec with OptionValues {
  private val tool = new UngappedAlignment {
    override def maxMismatchRate: Double = 0.25 // 1 in 4!
  }

  "UngappedAlignment.numMismatches" should "count the number of mismatches (with ambiguity)" in {
    // NB: We have inspected the test cases for SequenceUtil.readBaseMatchesRefBaseWithAmbiguity in htjskd and found them
    // to be comprehensive, so we do not duplicate them here.

    case class TestCase(left: String, right: String, expectedNumMismatches: Int, maxMismatches: Double) {
      def bases: Array[Byte] = left.getBytes
      def primer: Array[Byte] = right.getBytes
    }
    Seq(
      TestCase("",          "",            0, 10d),
      TestCase("AAAAA",      "AAAAA",      0, 10d),
      TestCase("TAAAA",      "AAAAA",      1, 10d),
      TestCase("AACAA",      "AAAAA",      1, 10d),
      TestCase("AAAAG",      "AAAAA",      1, 10d),
      TestCase("TTTTT",      "AAAAA",      5, 10d),
      TestCase("AAAAATTTTT", "AAAAA",      0, 10d),
      TestCase("AAAAA",      "AAAAATTTTT", 0, 10d),
      TestCase("A",          "N",          0, 10d),
      TestCase("A",          "M",          0, 10d),
      TestCase("A",          "S",          1, 10d),
      TestCase("M",          "V",          0, 10d),
      TestCase("V",          "M",          1, 10d), // NB: the bases of V are not a sub-set of the bases of M
      TestCase("TTTTT",      "AAAAA",      2, 1d), // NB: the maximum # of mismatches is 1, so we stop at 2
      TestCase("TTTTT",      "AAAAA",      2, 1.99) // NB: the maximum # of mismatches is 1.99, so we stop at 2
    ).foreach { testCase =>
      tool.numMismatches(testCase.bases, testCase.primer, 0, testCase.maxMismatches) shouldBe testCase.expectedNumMismatches
    }

    val testCase = TestCase("TTAAA", "AAAAA", 0, 10d)
    tool.numMismatches(testCase.bases, testCase.primer, 0, 10d) shouldBe 2
    tool.numMismatches(testCase.bases, testCase.primer, 1, 10d) shouldBe 1
    tool.numMismatches(testCase.bases, testCase.primer, 2, 10d) shouldBe 0
  }

  "UngappedAlignment.mismatchAlign" should "align a SamRecord to a Primer" in {
    val unmapped    = new SamBuilder(readLength=10).addFrag(bases="A"*10, unmapped=true).value
    val mappedPlus  = new SamBuilder(readLength=10).addFrag(bases="A"*10, strand=SamBuilder.Plus).value
    val mappedMinus = new SamBuilder(readLength=10).addFrag(bases="A"*10, strand=SamBuilder.Minus).value

    val primer      = Primer("1", "1+", "AAAAAAAAAA", "chr1", 1, 10, true)
    val mmPrimer    = Primer("1", "1+", "AAAAATAAAA", "chr1", 1, 10, true)
    val shortPrimer = Primer("1", "1+", "AAAAA", "chr1", 1, 5, true)
    val longPrimer  = Primer("1", "1+", "AAAAAAAAAATTTTT", "chr1", 1, 15, true)
    val mmShort     = Primer("1", "1+", "AAAAT", "chr1", 1, 5, true)
    val mmLong      = Primer("1", "1+", "TTTTTAAAAAAAAAA", "chr1", 1, 15, true)

    // both unmapped and positive strand mapped reads match at the star
    Seq(unmapped, mappedPlus).foreach { rec =>
      tool.mismatchAlign(rec, primer).value shouldBe UngappedAlignmentPrimerMatch(primer, 0, Int.MaxValue)
      tool.mismatchAlign(rec, mmPrimer).value shouldBe UngappedAlignmentPrimerMatch(mmPrimer, 1, Int.MaxValue)
      tool.mismatchAlign(rec, shortPrimer).value shouldBe UngappedAlignmentPrimerMatch(shortPrimer, 0, Int.MaxValue)
      tool.mismatchAlign(rec, longPrimer).value shouldBe UngappedAlignmentPrimerMatch(longPrimer, 0, Int.MaxValue) // matches at the start of the primer
      tool.mismatchAlign(rec, mmShort).value shouldBe UngappedAlignmentPrimerMatch(mmShort, 1, Int.MaxValue) // matches at the start of the read! OK 1/5 < 1/4 mm rate
      tool.mismatchAlign(rec, mmLong) shouldBe 'empty // matches at the start of the read! NOK 4/5 > 1/4 mm rate

    }

    // negative matches at the end
    {
      tool.mismatchAlign(mappedMinus, primer).value shouldBe UngappedAlignmentPrimerMatch(primer, 0, Int.MaxValue)
      tool.mismatchAlign(mappedMinus, mmPrimer).value shouldBe UngappedAlignmentPrimerMatch(mmPrimer, 1, Int.MaxValue) // mismatch rate ok
      tool.mismatchAlign(mappedMinus, shortPrimer).value shouldBe UngappedAlignmentPrimerMatch(shortPrimer, 0, Int.MaxValue)
      tool.mismatchAlign(mappedMinus, longPrimer) shouldBe 'empty // matches at the end of the primer! NOK 4/10 > 1/4 mm rate
      tool.mismatchAlign(mappedMinus, mmShort).value shouldBe UngappedAlignmentPrimerMatch(mmShort, 1, Int.MaxValue) // matches at the end of the read! OK 1/5 < 1/4 mm rate
      tool.mismatchAlign(mappedMinus, mmLong).value shouldBe UngappedAlignmentPrimerMatch(mmLong, 0, Int.MaxValue) // matches at the end of the primer! OK 1/5 < 1/4 mm rate
    }
  }
}

class LocationBasedPrimerMatcherTest extends UnitSpec with OptionValues {
  import PrimerMatcherTest._

  private def newMatcher(slop: Int, ignorePrimerStrand: Boolean, maxMismatchRate: Double = 1.0) = {
    new LocationBasedPrimerMatcher(
      primers            = primers,
      slop               = slop,
      ignorePrimerStrand = ignorePrimerStrand,
      maxMismatchRate    = maxMismatchRate
    )
  }

  private val builder = new SamBuilder(20)

  private def toMatch(primer: Primer, numMismatches: Int) = LocationBasedPrimerMatch(primer, numMismatches)

  "LocationBasedPrimerMatcher.getOverlaps" should "find overlaps when considering primer strand" in {
    val matcher = newMatcher(0, false)
    primers.foreach { primer =>
      val frag = fragFromPrimer(primer, builder)
      matcher.getOverlaps(frag).toList should contain theSameElementsInOrderAs Seq(primer)
    }
  }

  it should "find overlaps without considering primer strand" in {
    val matcher = newMatcher(0, true)

    primers.foreach { primer =>
      val newPrimer = primer.copy(positive_strand=primer.negativeStrand)
      val frag      = fragFromPrimer(newPrimer, builder)
      val interval  = new Interval(frag.refName, frag.start, frag.end)
      val expected  = primers
        .filter { p => p.overlaps(interval) }
        .map { p =>
          // since the strand was reversed, the bases will be reversed for mismatch alignment
          if (p.positiveStrand == primer.positiveStrand) p.copy(reverseComplementBases = true) else p
        }
      matcher.getOverlaps(frag).toList should contain theSameElementsInOrderAs expected
    }
  }

  it should "find overlaps with some slop" in {
    val matcher = newMatcher(5, false)
    primers.foreach { primer =>
      val frag = fragFromPrimer(primer, builder, startAdjust=5)
      matcher.getOverlaps(frag).toList should contain theSameElementsInOrderAs Seq(primer)
    }
  }


  it should "not find overlaps when start is beyond slop" in {
    val matcher = newMatcher(5, false)
    primers.foreach { primer =>
      val frag = fragFromPrimer(primer, builder, startAdjust=6)
      matcher.getOverlaps(frag).toList shouldBe 'empty
    }
  }

  "LocationBasedPrimerMatcher.find" should "not find primer matches when the mismatch rate is too high" in {
    val matcher = newMatcher(0, false, 0.0)
    primers.foreach { primer =>
      val frag = fragFromPrimer(primer, builder)
      frag.bases = "N"*builder.readLength // So many mismatches
      matcher.find(frag) shouldBe 'empty
    }
  }

  it should "find primer matches when considering primer strand" in {
    val matcher = newMatcher(0, false)
    primers.foreach { primer =>
      val frag = fragFromPrimer(primer, builder)
      matcher.find(frag).value shouldBe toMatch(primer, 0)
    }
  }

  it should "find primer matches without considering primer strand" in {
    val matcher = newMatcher(0, true)
    // test only the first four as the primers in each pair don't overlap each other
    primers.take(4).foreach { primer =>
      val newPrimer = primer.copy(positive_strand=primer.negativeStrand, reverseComplementBases=true)
      val frag      = fragFromPrimer(newPrimer, builder)
      matcher.find(frag).value shouldBe toMatch(primer.copy(reverseComplementBases=true), 0)
    }
  }

  it should "find prioritize primers on the same strand" in {
    val matcher = newMatcher(0, true)
    // this primer pair has primers map to the same location, so priorize the match based on strand (negative)
    primers.slice(4, 6).foreach { primer =>
      val frag        = fragFromPrimer(primer, builder)
      val primerMatch = matcher.mismatchAlign(frag, primer).value
      matcher.find(frag).value shouldBe LocationBasedPrimerMatch(primerMatch.primer, primerMatch.numMismatches)
    }
  }

  it should "fail if there are more than one overlapping primers" in {
    // take the first primer pair and shift it over by a base.  Voila!
    val newPrimers = primers.take(2) ++ primers.take(2).map(p => p.copy(start=p.start+1, end=p.end+1))
    val matcher = new LocationBasedPrimerMatcher(
      primers            = newPrimers,
      slop               = 5,
      ignorePrimerStrand = false,
      maxMismatchRate    = 1d
    )
    primers.take(2).foreach { primer =>
      val frag = fragFromPrimer(primer, builder)
      an[Exception] should be thrownBy matcher.find(frag)
    }
  }
}

class UngappedAlignmentBasedPrimerMatcherTest extends UnitSpec with OptionValues {

  import PrimerMatcherTest._

  private def newMatcher(slop: Int, ignorePrimerStrand: Boolean, maxMismatchRate: Double = 1.0) = {
    new UngappedAlignmentBasedPrimerMatcher(
      primers            = primers,
      slop               = slop,
      ignorePrimerStrand = ignorePrimerStrand,
      maxMismatchRate    = maxMismatchRate
    )
  }

  private val builder = new SamBuilder(20)

  private def toMatch(numMismatches: Int) = UngappedAlignmentPrimerMatch(null, numMismatches, Int.MaxValue)

  "UngappedAlignmentBasedPrimerMatcher.getBestUngappedAlignment" should "return None if no alignments were given" in {
    val matcher = newMatcher(0, false)
    matcher.getBestUngappedAlignment(Seq.empty) shouldBe 'empty
  }

  it should "handle a single alignment" in {
    val matcher = newMatcher(0, false)
    matcher.getBestUngappedAlignment(Seq(toMatch(0))).value shouldBe toMatch(0)
  }

  it should "handle two equally likely alignments" in {
    val matcher = newMatcher(0, false)
    matcher.getBestUngappedAlignment(Seq(toMatch(0), toMatch(0))).value shouldBe toMatch(0).copy(nextNumMismatches=0)
  }

  it should "handle multiple alignments with one with the fewest mismatches" in {
    val matcher = newMatcher(0, false)
    matcher.getBestUngappedAlignment(Seq(toMatch(0), toMatch(1), toMatch(3), toMatch(1))).value shouldBe toMatch(0).copy(nextNumMismatches=1)
  }

  "UngappedAlignmentBasedPrimerMatcher.find" should "find primer matches when considering primer strand" in {
    val matcher = newMatcher(0, false)
    primers.foreach { primer =>
      val frag        = fragFromPrimer(primer, builder)
      val primerMatch = matcher.find(frag).value
      primerMatch.primer shouldBe primer
      primerMatch.numMismatches shouldBe 0
    }
  }

  it should "find primer matches without considering primer strand" in {
    val matcher = newMatcher(0, true)
    // test only the first four as the primers in each pair don't overlap each other
    primers.take(4).foreach { primer =>
      val frag      = fragFromPrimer(primer, builder)
      val primerMatch = matcher.find(frag).value
      primerMatch.primer.sequence shouldBe primer.sequence
      primerMatch.numMismatches shouldBe 0
    }
  }

  it should "find primer matches with unmapped reads" in {
    val matcher = newMatcher(0, false)
    primers.foreach { primer =>
      val frag        = fragFromPrimer(primer.copy(positive_strand=true), builder) // use + strand so we set start right
      frag.unmapped   = true
      if (primer.negativeStrand) {
        // since we used + strand above, we still need to reverse complement complement..
        val revcomp = SequenceUtil.reverseComplement(primer.sequence)
        frag.bases = revcomp + ("N" * (builder.readLength - primer.length))
        val primerMatch = matcher.find(frag).value
        // NB: the primer returned may have been matched using the reverse complement of the primer, since unmapped reads
        // could match the negative genomic strand.
        new String(primerMatch.primer.bases) shouldBe revcomp
        primerMatch.numMismatches shouldBe 0
      }
      else {
        val primerMatch = matcher.find(frag).value
        primerMatch.primer.sequence shouldBe primer.sequence
        primerMatch.numMismatches shouldBe 0
      }
    }
  }
}

class GappedAlignmentBasedPrimerMatcherTest extends UnitSpec with OptionValues {

  import PrimerMatcherTest._

  private val matchScore    = 1
  private val mismatchScore = -3
  private val gapOpen       = -6
  private val gapExtend     = -1
  private val aligner       = Aligner(matchScore, mismatchScore, gapOpen, gapExtend, mode=Mode.Glocal)

  private def newMatcher(slop: Int, ignorePrimerStrand: Boolean, maxMismatchRate: Double = 1.0) = {
    new GappedAlignmentBasedPrimerMatcher(
      primers               = primers,
      slop                  = slop,
      ignorePrimerStrand    = ignorePrimerStrand,
      aligner               = aligner,
      minAlignmentScoreRate = 0d
    )
  }

  private def toMatch(score: Int) = GappedAlignmentPrimerMatch(primers(0), score, 0)

  {
    val matcher = newMatcher(0, false)

    def toAlignmentAndPrimer(score: Int) = {
      matcher.AlignmentAndPrimer(
        alignment = new Alignment(null, null, 0, 0, null, score),
        primer    = primers(0)
      )
    }

    "GappedAlignmentBasedPrimerMatcher.getBestGappedAlignment" should "return None if no alignments were given" in {
      matcher.getBestGappedAlignment(Seq.empty) shouldBe 'empty
    }

    it should "handle a single alignment" in {
      matcher.getBestGappedAlignment(Seq(toAlignmentAndPrimer(0))).value shouldBe toMatch(0)
    }

    it should "handle two equally likely alignments" in {
      matcher.getBestGappedAlignment(Seq(toAlignmentAndPrimer(0), toAlignmentAndPrimer(0))).value shouldBe toMatch(0).copy(secondBestScore=0)
    }

    it should "handle multiple alignments with one with the fewest mismatches" in {
      matcher.getBestGappedAlignment(Seq(toAlignmentAndPrimer(0), toAlignmentAndPrimer(1), toAlignmentAndPrimer(3), toAlignmentAndPrimer(1))).value shouldBe toMatch(0).copy(secondBestScore=1)
    }
  }
}


