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

package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.coord.LocatableOrderingTest._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.util.{Interval, Locatable}

/** Companion object for [[LocatableOrderingTest]]. */
object LocatableOrderingTest {

  /** The reference sequence name for chromosome 1. */
  private val Chr1: String = "chr1"

  /** The reference sequence name for chromosome 2. */
  private val Chr2: String = "chr2"

  /** The sequence dictionary for the ordering. */
  private val Dict: SequenceDictionary = SequenceDictionary(SequenceMetadata(Chr1, length = 10_000), SequenceMetadata(Chr2, length = 100_000))

  private val Builder:SamBuilder = new SamBuilder(sd=Some(Dict))
  Builder.addPair()
}

/** Unit tests for [[NormalizeCoverage]]. */

class NormalizeCoverageTest extends UnitSpec {





}
