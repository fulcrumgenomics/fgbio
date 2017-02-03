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
 *
 */

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

/**
  * Tests for SimpleCounter.
  */
class SimpleCounterTest extends UnitSpec {

  "SimpleCounter" should "count a comparable type" in {
    val counter = new SimpleCounter[Int]()
    counter.count(10, 5)
    counter.count(9, 4)
    counter.countOf(10) shouldBe 5
    counter.countOf(9) shouldBe 4
    counter.iterator.toList should contain theSameElementsAs List((9, 4), (10, 5))
    counter.total shouldBe 9
  }

  it should "count a type that does not implement comparable" in {
    val counter = new SimpleCounter[Any]()
    counter.count(10, 5)
    counter.count("10.5", 4)
    counter.countOf(10) shouldBe 5
    counter.countOf("10.5") shouldBe 4
    counter.iterator.toList should contain theSameElementsAs List(("10.5", 4), (10, 5))
    counter.total shouldBe 9
  }
}
