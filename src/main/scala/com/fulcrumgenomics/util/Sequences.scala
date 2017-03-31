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

import com.fulcrumgenomics.FgBioDef._
/**
  * Utility methods for working with DNA or RNA sequences
  */
object Sequences {
  /** Common contig/chrom names for non-autosomal sequences in mammals. */
  val CommonNonAutosomalContigNames = Seq("M", "chrM", "MT", "X", "chrX", "Y", "chrY")

  /** Counts the number of mismatches between two sequences of the same length. */
  def countMismatches(s1: String, s2: String): Int = {
    require(s1.length == s2.length, s"Cannot count mismatches in strings of differing lengths: $s1 $s2")

    var count = 0
    forloop (from=0, until=s1.length) { i =>
      val a = Character.toUpperCase(s1.charAt(i))
      val b = Character.toUpperCase(s2.charAt(i))
      if (a != b) count += 1
    }

    count
  }

  /**
    * Returns the 0-based index of, and the length of, the longest homopolymer in a non-empty string.
    * If there are multiple homopolymers of the same length the _first_ one is returned.
    *
    * @param s a sequence
    * @return a tuple of (0-based index, length) of the longest homopolymer in the string
    */
  def longestHomopolymer(s: String): (Int, Int) = {
    val xs = s.toUpperCase
    val len = xs.length
    assert(len > 0, "Cannot compute longest homopolymer of a zero length string")

    var foundIndex = 0
    var foundLength = 1

    forloop (from=0, until=len) { i =>
      val x = xs.charAt(i)
      var j = i+1
      while (j < len && xs.charAt(j) == x) j += 1
      val newLength = j - i
      if (newLength > foundLength) {
        foundIndex = i
        foundLength = newLength
      }
    }

    (foundIndex, foundLength)
  }
}
