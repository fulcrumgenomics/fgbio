/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

class AmpliconDetectorTest extends UnitSpec {

  // maxPrimerLength
  // find(refName: String, start: Int, end: Int, positiveStrand: Boolean)
  // - positive strand
  // - negative strand
  // - does not overlap amplicon
  // - overlaps amplicon, but not correct primer
  // - overlaps amplicon, but not near correct primer
  // - overlaps amplicon, but just within slop
  // - overlaps amplicon, exactly primer
  // def find(rec: SamRecord): Option[Amplicon]
  // - with clipped/without clipped (without clipped within slop, with clipped outside slop)
  // def findMate(rec: SamRecord): Option[Amplicon]
  // - with clipped/without clipped (without clipped within slop, with clipped outside slop)
  // def find(r1: SamRecord, r2: SamRecord): Option[Amplicon]
  // - not paired -> None
  // - read unmapped -> None
  // - mate unmapped -> None
  // - different contigs -> None
  // - not FR -> None
  // - r1(F)/r2(R)
  //   - with clipped/without clipped (without clipped within slop, with clipped outside slop)
  // - r1(R)/r2(F)
  //   - with clipped/without clipped (without clipped within slop, with clipped outside slop)
}
