/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec

class BetterBufferedIteratorTest extends UnitSpec {
  "BetterBufferedIterator" should "takeWhile and dropWhile without losing elements" in {
    val list = List(1,2,3,4,5,6,7,8,9)
    val xs = list.iterator.bufferBetter
    xs.takeWhile(_ < 5).toSeq shouldBe Seq(1,2,3,4)
    xs.toSeq shouldBe Seq(5,6,7,8,9)

    val ys = list.iterator.bufferBetter
    ys.dropWhile(_ <= 3)
    ys.takeWhile(_ <  7).toSeq shouldBe Seq(4,5,6)
    ys.dropWhile(_ <= 8)
    ys.takeWhile(_ => true).toSeq shouldBe Seq(9)
  }
}
