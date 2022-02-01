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

package com.fulcrumgenomics.bam.pileup

import com.fulcrumgenomics.bam.api.SamOrder.{Coordinate, Unknown}
import com.fulcrumgenomics.bam.pileup.PileupBuilder.BamAccessPattern
import com.fulcrumgenomics.bam.pileup.PileupBuilder.BamAccessPattern.{RandomAccess, Streaming}
import com.fulcrumgenomics.bam.pileup.PileupBuilderTest._
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

/** Companion object for [[PileupBuilderTest]]. */
object PileupBuilderTest {

  /** The name for chromosome 1. */
  private val Chr1: String = "chr1"

  /** The reference sequence index for chromosome 1. */
  private val Chr1Index: Int = 0

  /** The name for chromosome 2. */
  private val Chr2: String = "chr2"

  /** The reference sequence index for chromosome 2. */
  private val Chr2Index: Int = 1

  /** The read length to use for most tests. */
  private val ReadLength: Int = 50

  /** Small contig to use for tests. */
  private val SmallContigSequence: String = "ACGTGCAACGTGGCCTCAGTGATTATGCGC" * 100

  /** The custom SAM/BAM sequence dictionary to use for all tests. */
  private val TestSequenceDictionary: SequenceDictionary = {
    SequenceDictionary(
      SequenceMetadata(Chr1, SmallContigSequence.length, index = Chr1Index),
      SequenceMetadata(Chr2, SmallContigSequence.length, index = Chr2Index),
    )
  }
}

/** Unit tests for [[PileupBuilder]]. */
class PileupBuilderTest extends UnitSpec {

  BamAccessPattern.values.foreach { accessPattern =>
    s"PileupBuilder.apply with $accessPattern BAM access" should "build a simple pileup of only matched/mismatched fragment reads" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      1 to 5 foreach { _ => builder.addFrag(start = 5).foreach(_.bases = "A" * ReadLength) }
      1 to 4 foreach { _ => builder.addFrag(start = 5).foreach(_.bases = "C" * ReadLength) }
      1 to 3 foreach { _ => builder.addFrag(start = 5).foreach(_.bases = "G" * ReadLength) }
      1 to 2 foreach { _ => builder.addFrag(start = 5).foreach(_.bases = "T" * ReadLength) }
      1 to 1 foreach { _ => builder.addFrag(start = 5).foreach(_.bases = "N" * ReadLength) }

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)
      val pile   = piler.pileup(Chr1, 5).withoutIndels

      pile.depth shouldBe 5 + 4 + 3 + 2 + 1
      pile.count(_.base == 'A') shouldBe 5
      pile.count(_.base == 'C') shouldBe 4
      pile.count(_.base == 'G') shouldBe 3
      pile.count(_.base == 'T') shouldBe 2
      pile.count(_.base == 'N') shouldBe 1
      source.safelyClose()
      piler.safelyClose()
    }

    it should "pileup only matched/mismatched fragment reads while ignoring any that occur on a previous locus" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      // Occurs on a previous reference sequence
      1 to 5 foreach { _ => builder.addFrag(start = 55, contig = Chr1Index).foreach(_.bases = "A" * ReadLength) }

      // Occurs in a previous position and does not overlap requested pileup locus (they are abutting).
      1 to 5 foreach { _ => builder.addFrag(start = 5, contig = Chr2Index).foreach(_.bases = "A" * ReadLength) }

      1 to 5 foreach { _ => builder.addFrag(start = 55, contig = Chr2Index).foreach(_.bases = "A" * ReadLength) }
      1 to 4 foreach { _ => builder.addFrag(start = 55, contig = Chr2Index).foreach(_.bases = "C" * ReadLength) }
      1 to 3 foreach { _ => builder.addFrag(start = 55, contig = Chr2Index).foreach(_.bases = "G" * ReadLength) }
      1 to 2 foreach { _ => builder.addFrag(start = 55, contig = Chr2Index).foreach(_.bases = "T" * ReadLength) }
      1 to 1 foreach { _ => builder.addFrag(start = 55, contig = Chr2Index).foreach(_.bases = "N" * ReadLength) }

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)
      val pile   = piler.pileup(Chr2, 55).withoutIndels

      pile.depth shouldBe 5 + 4 + 3 + 2 + 1
      pile.count(_.base == 'A') shouldBe 5
      pile.count(_.base == 'C') shouldBe 4
      pile.count(_.base == 'G') shouldBe 3
      pile.count(_.base == 'T') shouldBe 2
      pile.count(_.base == 'N') shouldBe 1
      source.safelyClose()
      piler.safelyClose()
    }

    it should "build a pileup that contains all edge cases of indels" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addFrag(name = "q1", start = 101, cigar = "10M2D40M").foreach(_.bases = "A" * ReadLength)
      builder.addFrag(name = "q2", start = 101, cigar = "10M2I38M").foreach(_.bases = "C" * ReadLength)
      builder.addFrag(name = "q3", start = 101, cigar = "31M9I10M").foreach(_.bases = "G" * ReadLength)
      builder.addFrag(name = "q4", start = 101, cigar = "30M9D20M").foreach(_.bases = "T" * ReadLength)
      builder.addFrag(name = "q5", start = 141, cigar = "10I40M"  ).foreach(_.bases = "N" * ReadLength)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)

      // Before any of the indels, we should have 4 _bases_
      piler.pileup(Chr1, 105).baseIterator.size shouldBe 4

      // The locus before the first insertion and deletion (should include the insertion)
      val p1 = piler.pileup(Chr1, 110)
      p1.depth shouldBe 4
      p1.iterator.size shouldBe 5
      p1.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACGT"
      p1.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q2"
      p1.baseIterator.toSeq should contain theSameElementsAs p1.withoutIndels.iterator.toSeq

      // The locus of the first deleted base
      val p2 = piler.pileup(Chr1, 111)
      p2.depth shouldBe 4
      p2.iterator.size shouldBe 4
      p2.baseIterator.size shouldBe 3
      p2.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "CGT"
      p2.iterator.collect{ case x: DeletionEntry => x }.map(_.rec.name).next() shouldBe "q1"
      p2.baseIterator.toSeq should contain theSameElementsAs p2.withoutIndels.iterator.toSeq

      // The locus of the bigger indels
      val p3 = piler.pileup(Chr1, 131)
      p3.depth shouldBe 4
      p3.iterator.size shouldBe 5
      p3.baseIterator.size shouldBe 3
      p3.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACG"
      p3.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q3"
      p3.iterator.collect{ case x: DeletionEntry  => x }.map(_.rec.name).next() shouldBe "q4"
      p3.baseIterator.toSeq should contain theSameElementsAs p3.withoutIndels.iterator.toSeq

      // Locus where there is a read with a leading indel
      val p4 = piler.pileup(Chr1, 140)
      p4.depth shouldBe 4  // shouldn't count the leading indel
      p4.iterator.size shouldBe 5
      p4.baseIterator.size shouldBe 4
      p4.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACGT"
      p4.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q5"
      p4.baseIterator.toSeq should contain theSameElementsAs p4.withoutIndels.iterator.toSeq

      source.safelyClose()
      piler.safelyClose()
    }

    it should "allow a repeated locus query for a coordinate-maintaining call" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addFrag(start = 5).foreach(_.bases = "N" * ReadLength)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)

      val pile1 = piler.pileup(Chr1, 5).withoutIndels
      pile1.depth shouldBe 1
      pile1.count(_.base == 'N') shouldBe 1

      val pile2 = piler.pileup(Chr1, 5).withoutIndels
      pile2.depth shouldBe 1
      pile2.count(_.base == 'N') shouldBe 1

      source.safelyClose()
      piler.safelyClose()
    }

    it should "filter out reads below the minimum mapping quality" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addFrag(name = "q1", start = 101, mapq = 9)
      builder.addFrag(name = "q2", start = 101, mapq = 10)
      builder.addFrag(name = "q3", start = 101, mapq = 11)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, minMapQ = 10, mappedPairsOnly = false)
      val pile   = piler.pileup(Chr1, 105)

      pile.depth shouldBe 2
      pile.exists(_.rec.name == "q1") shouldBe false

      source.safelyClose()
      piler.safelyClose()
    }

    it should "filter out base entries below the minimum base quality" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), baseQuality = 19, sort = Some(Coordinate))

      builder.addFrag(name = "q1", start = 101)
      builder ++= new SamBuilder(
        readLength  = ReadLength,
        sd          = Some(TestSequenceDictionary),
        baseQuality = 20
      ).addFrag(name = "q2", start = 101)
      builder ++= new SamBuilder(
        readLength  = ReadLength,
        sd          = Some(TestSequenceDictionary),
        baseQuality = 21
      ).addFrag(name = "q3", start = 101)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)
      val pile   = piler.pileup(Chr1, 105)

      pile.depth shouldBe 2
      pile.exists(_.rec.name == "q1") shouldBe false

      source.safelyClose()
      piler.safelyClose()
    }

    it should "filter out reads that are not a part of a mapped pair" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addFrag(name = "q1", start = 101)
      builder.addPair(name = "q2", start1 = 101, start2 = 101, unmapped2 = true)
      builder.addPair(name = "q3", start1 = 101, start2 = 300)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern)
      val pile   = piler.pileup(Chr1, 105)

      pile.depth shouldBe 1
      pile.exists(r => r.rec.name == "q3" && r.rec.firstOfPair) shouldBe true

      source.safelyClose()
      piler.safelyClose()
    }

    it should "remove one half of each overlapping pair" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addPair(name = "q1", start1 = 100, start2 = 110)
      builder.addPair(name = "q2", start1 = 110, start2 = 100)
      builder.addPair(name = "q3", start1 = 50, start2 = 100)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern)
      val pile   = piler.pileup(Chr1, 125)

      pile.depth shouldBe 5
      pile.withoutOverlaps.depth shouldBe 3
      pile.withoutOverlaps.groupBy(_.rec.name).values.foreach(xs => xs.size shouldBe 1)

      source.safelyClose()
      piler.safelyClose()
    }

    it should "filter out records where a position is outside the insert for an FR pair" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addPair(name = "q2", start1 = 101, start2 = 100)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern)

      piler.pileup(Chr1, 100).depth shouldBe 0
      piler.pileup(Chr1, 101).depth shouldBe 2
      piler.pileup(Chr1, 100 + ReadLength - 1).depth shouldBe 2
      piler.pileup(Chr1, 101 + ReadLength - 1).depth shouldBe 0

      source.safelyClose()
      piler.safelyClose()
    }

    it should "include records where a position is outside the insert for an FR pair if asked to" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addPair(name = "q2", start1 = 101, start2 = 100)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, includeMapPositionsOutsideFrInsert = true)

      piler.pileup(Chr1, 100).depth shouldBe 1
      piler.pileup(Chr1, 101).depth shouldBe 2
      piler.pileup(Chr1, 100 + ReadLength - 1).depth shouldBe 2
      piler.pileup(Chr1, 101 + ReadLength - 1).depth shouldBe 1

      source.safelyClose()
      piler.safelyClose()
    }

    it should "not filter out records where a position is outside what might look like an 'insert' for a non-FR pair" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addPair(name = "q2", start1 = 101, start2 = 100, strand1 = Minus, strand2 = Plus)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern)

      piler.pileup(Chr1, 100).depth shouldBe 1
      piler.pileup(Chr1, 101).depth shouldBe 2
      piler.pileup(Chr1, 100 + ReadLength - 1).depth shouldBe 2
      piler.pileup(Chr1, 101 + ReadLength - 1).depth shouldBe 1

      source.safelyClose()
      piler.safelyClose()
    }

    it should "filter out records and entries with custom filters" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

      builder.addFrag(name = "q1", start = 100)
      builder.addFrag(name = "x2", start = 104)
      builder.addFrag(name = "q3", start = 108)
      builder.addFrag(name = "q4", start = 112)

      val source = builder.toSource
      val piler  = PileupBuilder(source, accessPattern = accessPattern, mappedPairsOnly = false)
        .withReadFilter(r => !r.name.startsWith("x"))
        .withEntryFilter(_.offset > 5)

      val pile = piler.pileup(Chr1, 115)

      pile.depth shouldBe 2
      pile.map(_.rec.name) should contain theSameElementsAs Seq("q1", "q3")

      source.safelyClose()
      piler.safelyClose()
    }

    it should "report the correct positions and offsets from pileup entries" in {
      val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), baseQuality = 35, sort = Some(Coordinate))

      builder.addPair(name = "q1", start1 = 101, start2 = 201, bases1 = "A" * ReadLength, bases2 = "C" * ReadLength)

      val source = builder.toSource

      val piler1 = PileupBuilder(source, accessPattern = accessPattern)
      val pile1  = piler1.pileup(Chr1, 105)
      val p1     = pile1.pile.head.asInstanceOf[BaseEntry]

      pile1.depth shouldBe 1
      p1.base shouldBe 'A'
      p1.baseInReadOrientation shouldBe 'A'
      p1.qual shouldBe  35
      p1.offset shouldBe 4
      p1.positionInRead shouldBe 5
      p1.positionInReadInReadOrder shouldBe 5

      val piler2 = PileupBuilder(source, accessPattern = accessPattern)
      val pile2  = piler2.pileup(Chr1, 205)
      val p2     = pile2.pile.head.asInstanceOf[BaseEntry]

      pile2.depth shouldBe 1
      p2.base shouldBe 'C'
      p2.baseInReadOrientation shouldBe 'G'
      p2.qual shouldBe  35
      p2.offset shouldBe 4
      p2.positionInRead shouldBe 5
      p2.positionInReadInReadOrder shouldBe 46

      source.safelyClose()
      piler1.safelyClose()
      piler2.safelyClose()
    }
  }

  it should "raise an exception if the input SAM source is not coordinate sorted when BAM access is `RandomAccess`" in {
    val builder  = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Unknown))
    val source   = builder.toSource
    require(!source.indexed)
    an[IllegalArgumentException] shouldBe thrownBy { PileupBuilder(source, accessPattern = RandomAccess) }
    source.safelyClose()
  }

  it should "raise an exception if the input SAM source is not coordinate sorted when BAM access is `Streaming`" in {
    val builder  = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Unknown))
    val source   = builder.toSource
    an[IllegalArgumentException] shouldBe thrownBy { PileupBuilder(source, accessPattern = Streaming) }
    source.safelyClose()
  }
}
