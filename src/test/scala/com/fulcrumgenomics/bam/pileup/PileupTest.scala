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

import com.fulcrumgenomics.FgBioDef.SafelyClosable
import com.fulcrumgenomics.bam.api.SamOrder.Coordinate
import com.fulcrumgenomics.bam.pileup.PileupTest._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

/** Companion object for [[PileupTest]]. */
object PileupTest {

  /** The name for chromosome 1. */
  private val Chr1: String = "chr1"

  /** The reference sequence index for chromosome 1. */
  private val Chr1Index: Int = 0

  /** The read length to use for most tests. */
  private val ReadLength: Int = 50

  /** Small contig to use for tests. */
  private val SmallContigSequence: String = "ACGTGCAACGTGGCCTCAGTGATTATGCGC" * 100

  /** The custom SAM/BAM sequence dictionary to use for all tests. */
  private val TestSequenceDictionary: SequenceDictionary = {
    SequenceDictionary(SequenceMetadata(Chr1, SmallContigSequence.length, index = Chr1Index))
  }
}

/** Unit tests for [[Pileup]]. */
class PileupTest extends UnitSpec {

  "BaseEntry" should "report the correct offsets/positions/bases" in {
    val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), baseQuality = 35, sort = Some(Coordinate))
    builder.addPair(name = "q1", start1 = 101, start2 = 201, bases1 = "A" * ReadLength, bases2 = "C" * ReadLength)

    val source1 = builder.toSource
    val piler1  = PileupBuilder(source1, mappedPairsOnly = true)
    val pile1   = piler1.pileup(Chr1, 105)
    pile1.depth shouldBe 1
    val p1 = pile1.pile.head.asInstanceOf[BaseEntry]
    p1.base shouldBe 'A'
    p1.baseInReadOrientation shouldBe 'A'
    p1.qual shouldBe  35
    p1.offset shouldBe 4
    p1.positionInRead shouldBe 5
    p1.positionInReadInReadOrder shouldBe 5
    source1.safelyClose()
    piler1.safelyClose()

    val source2 = builder.toSource
    val piler2  = PileupBuilder(source2, mappedPairsOnly = true)
    val pile2   = piler2.pileup(Chr1, 205)
    pile2.depth shouldBe 1
    val p2 = pile2.pile.head.asInstanceOf[BaseEntry]
    p2.base shouldBe 'C'
    p2.baseInReadOrientation shouldBe 'G'
    p2.qual shouldBe  35
    p2.offset shouldBe 4
    p2.positionInRead shouldBe 5
    p2.positionInReadInReadOrder shouldBe 46
    source2.safelyClose()
    piler2.safelyClose()
  }
}
