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

import com.fulcrumgenomics.testing.UnitSpec

/** Tests for the delimited data parser. */
class DelimitedDataParserTest extends UnitSpec {
  def csv(lines: Seq[String], trim: Boolean=true) = new DelimitedDataParser(lines, delimiter=',', trimFields=trim)

  "DelimitedDataParser" should "parse some basic fields" in {
    val parser = csv(Seq("zero,one,two", "1,foo,5.0", "2,bar,10.0"))
    parser.hasNext shouldBe true
    val row1 = parser.next()
    row1[String](0) shouldBe "1"
    row1[Int](0) shouldBe 1
    row1[Int]("zero") shouldBe 1
    row1[String](1) shouldBe "foo"
    row1[String]("one") shouldBe "foo"
    row1[Double](2) shouldBe 5.0
    row1[Double]("two") shouldBe 5.0

    parser.hasNext shouldBe true
    val row2 = parser.next()
    row2[Int]("zero") shouldBe 2
    row2[String]("one") shouldBe "bar"
    row2[Double]("two") shouldBe 10.0

    parser.hasNext shouldBe false
  }

  it should "handle consecutive delimiters correctly" in {
    val parser = csv(Seq("zero,one,two", "0,,2"))
    val row = parser.next()
    row[Int]("zero") shouldBe 0
    row[String]("one") shouldBe ""
    an[Exception] shouldBe thrownBy { row[Int]("one") }
    row.get[String]("one") shouldBe None
    row[Int]("two") shouldBe 2
  }

  it should "parse an empty set of lines" in {
    csv(Seq.empty[String]).hasNext shouldBe false
    csv(Seq("", "")).hasNext shouldBe false // blank lines are suppressed by default
  }

  it should "report the set of header fields" in {
    csv(Seq("foo,bar,splat ")).headers should contain theSameElementsInOrderAs Seq("foo", "bar", "splat")
  }

  it should "trim values when requested" in {
    val lines = Seq("one,two", " foo , bar")
    val row = csv(lines).next()
    row[String]("one") shouldBe "foo"
    row[String]("two") shouldBe "bar"

    val row2 = csv(lines, trim=false).next()
    row2[String]("one") shouldBe " foo "
    row2[String]("two") shouldBe " bar"
  }

  it should "throw an exception when lines have more or less values than the header" in {
    val p1 = csv(Seq("one,two,three", "one,two,three,four"))
    an[Exception] shouldBe thrownBy { p1.next() }

    val p2 = csv(Seq("one,two,three", "one,two"))
    an[Exception] shouldBe thrownBy { p2.next() }
  }

  it should "work with alternative delimiters" in {
    val parser = DelimitedDataParser(lines=Seq("one|two|three", "a|bee|c"), delimiter='|')
    val row = parser.next()
    row[String]("one") shouldBe "a"
    row[String]("two") shouldBe "bee"
    row[String]("three") shouldBe "c"
  }

  it should "handle blank lines anywhere in the file" in {
    val parser = csv(Seq("", "one,two", "", "", "1,2", "3,4", "", "5,6", ""))
    parser.headers should contain theSameElementsInOrderAs Seq("one", "two")
    val row1 = parser.next()
    val row2 = parser.next()
    val row3 = parser.next()
    parser.hasNext shouldBe false
    row1[Int]("one") shouldBe 1
    row1[Int]("two") shouldBe 2
    row2[Int]("one") shouldBe 3
    row2[Int]("two") shouldBe 4
    row3[Int]("one") shouldBe 5
    row3[Int]("two") shouldBe 6
  }
}
