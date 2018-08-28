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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
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
  //"LocationBasedPrimerMatcher.getOverlaps" should
  // - ignorePrimerStrand or enforce strand
  // - slop or no slop
  // - + strand, - strand

  //"LocationBasedPrimerMatcher.find" should
  // - find no primers
  // - find a single primer
  // - find multiple primers, prioritize on the forward strand, ok
  // - find multiple primers, prioritize on the forward strand, fail
  // - find multiple primers, fail
  // - ignorePrimerStrand, + strand read, matches a reverse strand read, ok with num mismatches
}

class UngappedAlignmentBasedPrimerMatcherTest extends UnitSpec with OptionValues {
  // getBestUngappedAlignment
  // - none if no alignments
  // - (numMM, Int.MaxValue) if one alignment
  // - (numMM, numMM) if two equally likely alignments
  // - (numMM, nextMM) if not equal

  // find
  // - + strand, - strand, unmapped
}

class GappedAlignmentBasedPrimerMatcherTest extends UnitSpec with OptionValues {
  // getBestGappedAlignment
  // - none if no alignments
  // - (score, alignmentscorerate * primer length) if one alignment
  // - (score, nextBestScore)
}


