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

package com.fulcrumgenomics.bam.pileup

import com.fulcrumgenomics.bam.api.SamOrder.Coordinate
import com.fulcrumgenomics.bam.pileup.StreamingPileupBuilderTest._
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

/** Companion object for [[StreamingPileupBuilderTest]]. */
object StreamingPileupBuilderTest {

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

/** Unit tests for [[StreamingPileupBuilder]]. */
class StreamingPileupBuilderTest extends UnitSpec {

  "StreamingPileupBuilder" should "raise an exception if asked to query a locus prior to the last locus" in {
    val builder = new SamBuilder(readLength = ReadLength, sd = Some(TestSequenceDictionary), sort = Some(Coordinate))

    builder.addFrag(start = 5).foreach(_.bases = "N" * ReadLength)

    val source   = builder.toSource
    val piler    = StreamingPileupBuilder(source, mappedPairsOnly = false)

    val pile1 = piler.pileup(Chr1, 6).withoutIndels
    pile1.depth shouldBe 1
    pile1.count(_.base == 'N') shouldBe 1

    an[IllegalArgumentException] shouldBe thrownBy { piler.pileup(Chr1, 5) }

    source.safelyClose()
    piler.safelyClose()
  }
}
