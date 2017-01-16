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
  * Tests for NumericCounter.
  */
class NumericCounterTest extends UnitSpec {

  "NumericCounter" should "operate on Ints" in {
    val counter = new NumericCounter[Int]()
    Stream.range(1, 10, 1).foreach { i => counter.count(i, 1) }
    counter.count(1)
    counter.iterator.sliding(2).foreach { case Seq(lower, upper) => lower should be < upper }
    counter.mean() shouldBe 4.6 +- 0.0001
    counter.stddev() shouldBe  2.875181 +- 0.0001
    counter.mode() shouldBe 1
    counter.median() shouldBe 4.5
    counter.mad() shouldBe 2.5
    counter.totalMass shouldBe 46
  }

  it should "operate on Doubles" in {
    val counter = new NumericCounter[Double]()
    Stream.range(1, 10, 1).foreach { i => counter.count(i / 10.0, 1) }
    counter.count(1 / 10.0)
    counter.iterator.sliding(2).foreach { case Seq(lower, upper) => lower should be < upper }
    counter.mean() shouldBe 0.46 +- 0.0001
    counter.stddev() shouldBe  0.2875181 +- 0.0001
    counter.mode() shouldBe 0.1
    counter.median() shouldBe 0.45
    counter.mad() shouldBe 0.25 +- 0.0001
    counter.totalMass shouldBe 4.6
  }

  it should "handle when no objects are counted" in {
    val counter = new NumericCounter[Double]()
    counter.mean() shouldBe 0.0
    counter.stddev() shouldBe 0.0
    counter.mode() shouldBe 0.0
    counter.median() shouldBe 0.0
    counter.mad() shouldBe 0.0
    counter.totalMass shouldBe 0
  }

  it should "handle when a single object is counted" in {
    val counter = new NumericCounter[Double]()
    counter.count(10.0)
    counter.mean() shouldBe 10.0
    counter.stddev() shouldBe 0.0
    counter.mode() shouldBe 10.0
    counter.median() shouldBe 10.0
    counter.mad() shouldBe 0.0
    counter.totalMass shouldBe 10.0
  }

  it should "handle when two objects are counted" in {
    val counter = new NumericCounter[Double]()
    counter.count(10.0)
    counter.count(20.0)
    counter.mean() shouldBe 15.0
    counter.stddev() shouldBe 7.071 +- 0.0001
    counter.mode() shouldBe 10.0
    counter.median() shouldBe 15.0
    counter.mad() shouldBe 5.0
    counter.totalMass shouldBe 30.0
  }

  it should "handle updating the count and re-computing statistics" in {
    val counter = new NumericCounter[Int]()
    counter.count(1)
    counter.mean() shouldBe 1
    counter.stddev() shouldBe 0.0
    counter.mode() shouldBe 1
    counter.median() shouldBe 1
    counter.count(9)
    counter.mean() shouldBe 5
    counter.stddev() shouldBe 5.656854 +- 0.00001
    counter.mode() shouldBe 1
    counter.median() shouldBe 5.0
    counter.mad() shouldBe 4.0
    counter.totalMass shouldBe 10
  }

  it should "compute mean and standard deviation together and separately" in {
    val counter = new NumericCounter[Int]()
    Range(0, 100).foreach(counter.count)
    val (mean, stddev) = counter.meanAndStddev()
    mean shouldBe counter.mean()
    stddev shouldBe counter.stddev(m=mean)
  }
}
