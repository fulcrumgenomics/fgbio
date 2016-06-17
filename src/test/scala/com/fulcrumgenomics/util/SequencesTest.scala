/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

/**
  * Tests for the Sequences utility object.
  */
class SequencesTest extends UnitSpec {
  "Sequences.countMismatches" should "return 0 when comparing a string with itself" in {
    val seq = "ACGTAGCTGACTGCA"
    Sequences.countMismatches(seq, seq) shouldBe 0
  }

  it should "throw an exception when comparing sequences of different lengths" in {
    an[IllegalArgumentException] should be thrownBy Sequences.countMismatches("ACGT", "ACGTGG")
  }

  it should "correct count mismatches" in {
    Sequences.countMismatches("ACGACTATATGCAT", "acgactatatgcat") shouldBe 0
    Sequences.countMismatches("ACGACTATATGCAT", "acgactaNatgcat") shouldBe 1
    Sequences.countMismatches("ACGACTATATGCAT", "acgactataGTcat") shouldBe 2
    Sequences.countMismatches("ACGACTATATGCAT", "NNNNNNNNNNNNNN") shouldBe 14
  }
}
