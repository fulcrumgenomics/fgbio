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

class AmpliconTest extends UnitSpec {

  "Amplicon" should "store an amplicon" in {
    val amplicon = Amplicon("chr1", 1, 2, 4, 8, id=Some("id"))
    amplicon.leftPrimerLength shouldBe 2
    amplicon.rightPrimerLength shouldBe 5
    amplicon.longestPrimerLength shouldBe 5
    amplicon.leftPrimerLocation shouldBe "chr1:1-2"
    amplicon.rightPrimerLocation shouldBe "chr1:4-8"
    amplicon.ampliconLocation shouldBe "chr1:1-8"
    amplicon.identifier shouldBe "id"
  }

  it should "round trip to disk amplicons with ids" in {
    val amplicons = Range.inclusive(start=1, end=100).map { i =>
      Amplicon(f"chr${1 + (i%10)}", 1, 2+i, 4+i, 8+i, id=Some(f"id-$i"))
    }
    val path = makeTempFile("some.", "metrics.txt")
    Metric.write(path=path, amplicons)
    Metric.read[Amplicon](path=path) should contain theSameElementsInOrderAs amplicons
  }

  it should "round trip to disk amplicons without ids" in {
    val amplicons = Range.inclusive(start=1, end=100).map { i =>
      Amplicon(f"chr${1 + (i%10)}", 1, 2+i, 4+i, 8+i, id=None)
    }
    val path = makeTempFile("some.", "metrics.txt")
    Metric.write(path=path, amplicons)
    Metric.read[Amplicon](path=path) should contain theSameElementsInOrderAs amplicons
  }
}
