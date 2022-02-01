/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.coord

import com.fulcrumgenomics.coord.LocatableOrderingTest._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.{Interval, Locatable}

/** Companion object for [[LocatableOrderingTest]]. */
object LocatableOrderingTest {

  /** The reference sequence name for chromosome 1. */
  private val Chr1: String = "chr1"

  /** The reference sequence name for chromosome 2. */
  private val Chr2: String = "chr2"

  /** The sequence dictionary for the ordering. */
  private val Dict: SequenceDictionary = SequenceDictionary(SequenceMetadata(Chr1), SequenceMetadata(Chr2))

  /** The ordering of the given <Dict>. */
  private val Ordering: Ordering[Locatable] = LocatableOrdering(Dict)
}

/** Unit tests for [[LocatableOrdering]]. */
class LocatableOrderingTest extends UnitSpec {

  "LocatableOrdering" should "know when two locatables are equivalent" in {
    val interval1 = new Interval(Chr1, 1, 1)
    val interval2 = new Interval(Chr1, 1, 1)
    Ordering.compare(interval1, interval2) shouldBe 0
  }

  it should "know when one Locatable is more 'left' than another Locatable on the same contig" in {
    val interval1 = new Interval(Chr1, 1, 1)
    val interval2 = new Interval(Chr1, 2, 2)
    Ordering.compare(interval1, interval2) should be < 0
  }

  it should "know when one Locatable is more 'right' than another Locatable on the same contig" in {
    val interval1 = new Interval(Chr1, 2, 2)
    val interval2 = new Interval(Chr1, 1, 2)
    Ordering.compare(interval1, interval2) should be > 0
  }

  it should "know when one Locatable is more 'left' than another Locatable on a further 'right' contig" in {
    val interval1 = new Interval(Chr1, 1, 1)
    val interval2 = new Interval(Chr2, 1, 1)
    Ordering.compare(interval1, interval2) should be < 0
  }

  it should "know when one Locatable is more 'right' than another Locatable on a further 'left' contig" in {
    val interval1 = new Interval(Chr2, 1, 1)
    val interval2 = new Interval(Chr1, 1, 1)
    Ordering.compare(interval1, interval2) should be > 0
  }

  it should "raise an exception when any of the locatables are aligned to contigs that don't exist" in {
    val interval1 = new Interval(Chr1, 1, 1)
    val interval2 = new Interval("ChrDoesNotExist", 1, 1)
    a[NoSuchElementException] shouldBe thrownBy { Ordering.compare(interval1, interval2) }
  }

  it should "order genomic locatables first by contig, then start, then end" in {
    val intervals = Seq(
      new Interval(Chr2, 2, 5),
      new Interval(Chr2, 2, 9),
      new Interval(Chr1, 2, 5),
      new Interval(Chr1, 2, 5),
      new Interval(Chr1, 4, 5)
    )
    val expected  = Seq(
      new Interval(Chr1, 2, 5),
      new Interval(Chr1, 2, 5),
      new Interval(Chr1, 4, 5),
      new Interval(Chr2, 2, 5),
      new Interval(Chr2, 2, 9)
    )
    intervals.min(Ordering) shouldBe new Interval(Chr1, 2, 5)
    intervals.max(Ordering) shouldBe new Interval(Chr2, 2, 9)
    intervals.sorted(Ordering) should contain theSameElementsInOrderAs expected
  }
}
