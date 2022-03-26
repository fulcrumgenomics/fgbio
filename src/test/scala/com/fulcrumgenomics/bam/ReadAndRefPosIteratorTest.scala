/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.ReadAndRefPosIterator.ReadAndRefPos
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class ReadAndRefPosIteratorTest extends UnitSpec {

  private def test(start: Int, cigar: String, readPositions: Seq[Int], refPositions: Seq[Int]): Unit = {
    val builder   = new SamBuilder(Cigar(cigar).lengthOnQuery)
    val frag      = builder.addFrag(start=start, cigar=cigar).value
    val positions = new ReadAndRefPosIterator(rec=frag).toIndexedSeq
    positions.map(_.readPos) should contain theSameElementsInOrderAs readPositions
    positions.map(_.refPos) should contain theSameElementsInOrderAs refPositions
  }

  private def ri(start: Int, end: Int): Range = Range.inclusive(start=start, end=end)

  "ReadAndRefPosIterator" should "return the read and reference positions" in {
    test(start=10, cigar="10M",    readPositions=ri(1, 10),             refPositions=ri(10, 19))
    test(start=10, cigar="1S9M",   readPositions=ri(2, 10),             refPositions=ri(10, 18))
    test(start=10, cigar="9M1S",   readPositions=ri(1, 9),              refPositions=ri(10, 18))
    test(start=10, cigar="5M2D5M", readPositions=ri(1, 10),             refPositions=ri(10, 14) ++ ri(17, 21))
    test(start=1, cigar="4M2I4M",  readPositions=ri(1, 4) ++ ri(7, 10), refPositions=ri(1, 8))
    test(start=6, cigar="10S15M5D5M5I65M5H", readPositions=ri(11, 30) ++ ri(36, 100), refPositions=ri(6, 20) ++ ri(26, 95))
  }

  it should "return the read and reference positions for read windows that overlap the alignment" in {
    val builder   = new SamBuilder(100)
    val frag      = builder.addFrag(start=100, cigar="10S15M5D5M5I65M").value

    // read window is soft-clipped, nothing returned
    new ReadAndRefPosIterator(rec=frag, minReadPos=1,  maxReadPos=10) shouldBe empty

    // read window contains the first mapped base after leading soft-clips
    new ReadAndRefPosIterator(rec=frag, minReadPos=1,  maxReadPos=11).toIndexedSeq should contain theSameElementsInOrderAs Seq(ReadAndRefPos(11, 100))

    // read window includes the base before and after the deletion
    new ReadAndRefPosIterator(rec=frag, minReadPos=25, maxReadPos=26).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(25, 114),
      ReadAndRefPos(26, 120)
    )

    // read window that includes the base before and after the insertion
    new ReadAndRefPosIterator(rec=frag, minReadPos=30, maxReadPos=36).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(30, 124),
      ReadAndRefPos(36, 125)
    )

    // read window that goes beyond the read
    new ReadAndRefPosIterator(rec=frag, minReadPos=100, maxReadPos=101).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(100, 189),
    )
  }

  it should "return the read and reference positions for reference windows that overlap the alignment" in {
    val builder   = new SamBuilder(100)
    val frag      = builder.addFrag(start=100, cigar="10S15M5D5M5I65M").value

    // before the alignment and the first mapped base
    new ReadAndRefPosIterator(rec=frag, minRefPos=1,   maxRefPos=100).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(11, 100),
    )

    // reference window includes the base before and after the deletion
    new ReadAndRefPosIterator(rec=frag, minRefPos=114, maxRefPos=120).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(25, 114),
      ReadAndRefPos(26, 120)
    )

    // reference window that includes the base before and after the insertion
    new ReadAndRefPosIterator(rec=frag, minRefPos=124, maxRefPos=125).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(30, 124),
      ReadAndRefPos(36, 125)
    )

    // reference window that goes beyond the read
    new ReadAndRefPosIterator(rec=frag, minRefPos=189, maxRefPos=200).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(100, 189),
    )
  }

  it should "return the read and reference positions for both a read and reference window" in {
    val builder   = new SamBuilder(100)
    val frag      = builder.addFrag(start=100, cigar="10S15M5D5M5I65M").value

    // The read window includes 6b before the deletion, and 1bp after.  The reference window includes 1bp before the
    // deletion, and 6bp after.  The intersection of the two is just 1bp before and after the deletion.
    new ReadAndRefPosIterator(
      rec=frag, minReadPos=20, maxReadPos=26, minRefPos=114, maxRefPos=125
    ).toIndexedSeq should contain theSameElementsInOrderAs Seq(
      ReadAndRefPos(25, 114),
      ReadAndRefPos(26, 120)
    )
  }

  it should "return an empty iterator for read or reference windows that do not overlap the alignment" in {
    val builder   = new SamBuilder(100)
    val frag      = builder.addFrag(start=100, cigar="100M").value
    // ref window OOB
    new ReadAndRefPosIterator(rec=frag,  minRefPos=1,   maxRefPos=99) shouldBe empty
    new ReadAndRefPosIterator(rec=frag,  minRefPos=200, maxRefPos=299) shouldBe empty
    new ReadAndRefPosIterator(rec=frag,  minRefPos=101, maxRefPos=100) shouldBe empty
    // read window OOB
    new ReadAndRefPosIterator(rec=frag, minReadPos= -1, maxReadPos=0) shouldBe empty
    new ReadAndRefPosIterator(rec=frag, minReadPos=201, maxReadPos=299) shouldBe empty
    new ReadAndRefPosIterator(rec=frag, minReadPos=101, maxReadPos=100) shouldBe empty
  }

  it should "return the read and reference positions for the first or last base" in {
    val builder   = new SamBuilder(100)
    val frag      = builder.addFrag(start=100, cigar="50M10D50M").value

    new ReadAndRefPosIterator(rec=frag,
      minReadPos=1, minRefPos=1, // before the alignment read/ref start
      maxReadPos=1, maxRefPos=100 // first mapped base in the alignments
    ).toIndexedSeq should contain theSameElementsInOrderAs Seq(ReadAndRefPos(1, 100))

    new ReadAndRefPosIterator(rec=frag,
      minReadPos=100, minRefPos=209, // last mapped base in the alignment
      maxReadPos=10000, maxRefPos=10000 // after the alignment read/ref end
    ).toIndexedSeq should contain theSameElementsInOrderAs Seq(ReadAndRefPos(100, 209))

  }

  "ReadAndRefPosIterator.apply" should "return read and reference positions for a perfectly overlapping read pair" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=1, cigar1="100M", start2=1, cigar2="100M")
    val positions1   = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2   = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=100)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=100)

    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=100)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=100)
  }

  it should "return an empty iterator for non-overlapping mates" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=1, cigar1="100M", start2=101, cigar2="100M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Seq()
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Seq()
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Seq()
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Seq()
  }

  it should "return read and reference positions for a partially overlapping read pair" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=1, cigar1="100M", start2=50, cigar2="100M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=50, end=100)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=50, end=100)
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=51)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=50, end=100)
  }

  it should "return read and reference positions for a read pair overlapping 1bp" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=1, cigar1="100M", start2=100, cigar2="100M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Seq(100)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Seq(100)
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Seq(1)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Seq(100)
  }

  it should "return read and reference positions for a read pair with soft-clipping in the overlap" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=100, cigar1="20S80M", start2=150, cigar2="20S80M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=71, end=100)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=150, end=179)
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=21, end=50)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=150, end=179)
  }

  it should "return read and reference positions for a read pair where each read extends past the other's end" in {
    val builder     = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=100, cigar1="10S90M", start2=90, cigar2="90M10S")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    // R1 read: first 10bp of read is soft-clipped, and last 10bp extend past mate's end
    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=11, end=90)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=100, end=179)
    // R2 read: first 10bp of read is before the mate's start, and last 10bp is soft-clipped
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=11, end=90)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=100, end=179)
  }

  it should "not return a value for an insertion" in {
    val builder     = new SamBuilder(10)
    val Seq(r1, r2) = builder.addPair(start1=10, cigar1="4M2I4M", start2=10, cigar2="10M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    positions1.map(_.readPos) should contain theSameElementsInOrderAs Seq(1, 2, 3, 4, 7, 8, 9, 10)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Seq(10, 11, 12, 13, 14, 15, 16, 17)
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=8)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Seq(10, 11, 12, 13, 14, 15, 16, 17)
  }

  it should "not return a value for a deletion" in {
    val builder     = new SamBuilder(10)
    val Seq(r1, r2) = builder.addPair(start1=10, cigar1="5M2D5M", start2=11, cigar2="10M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    // R1 read: last 2bp are past the end of it's mate
    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=2, end=9)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Seq(11, 12, 13, 14, 17, 18, 19, 20)
    // R2 read: all mapped bases
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=10)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=11, end=20)
  }

  it should "return values for a complicated cigar in an overlapping read pair" in {
    val builder     = new SamBuilder(100)
    val Seq(r1, r2) = builder.addPair(start1=1, cigar1="10S45M10I15M5D10M3I2M5S", start2=50, cigar2="100M")
    val positions1  = ReadAndRefPosIterator(r1, r2).toIndexedSeq
    val positions2  = ReadAndRefPosIterator(r2, r1).toIndexedSeq

    // R1 read: first 10b soft-clip (skip), next 45bp before mate start (skip), 10bp insertion (skip), 5M (out of 15M)
    //          before mate start (skip).  Next 10M (out of 15M) aligned (output), 5D skipped (advance ref), 10M used
    //          (output), 3I (not output), 2M used (output)
    positions1.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=70, end=90) ++ Seq(94, 95)
    positions1.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=50, end=60) ++ Range.inclusive(start=66, end=77)
    positions2.map(_.readPos) should contain theSameElementsInOrderAs Range.inclusive(start=1, end=28)
    positions2.map(_.refPos) should contain theSameElementsInOrderAs Range.inclusive(start=50, end=77)
  }
}
