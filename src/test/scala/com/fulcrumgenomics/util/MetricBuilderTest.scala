/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.testing.UnitSpec


case class MetricBuilderTestMetric(name: String, count: Long = 1) extends Metric

case class MetricBuilderNameMetric(first: String, second: String, age: Int) extends Metric

class MetricBuilderTest extends UnitSpec {
  private val builder = new MetricBuilder[MetricBuilderTestMetric]()

  "MetricBuilder.fromValues" should "build a metric from a list of values" in {
    builder.fromValues(Seq("foo", "2")) shouldBe MetricBuilderTestMetric(name="foo", count=2)
  }

  it should "fail if the number of values is different than expected" in {
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq())
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq("Foo"))
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq("1", "2", "3"))
  }

  it should "fail if the # of values are incorrect or wrong type" in {
    val builder = new MetricBuilder[MetricBuilderNameMetric]()
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq("one", "two"))
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq("one", "two", "3", "four"), lineNumber=Some(1))
    a[MetricBuilderException] should be thrownBy builder.fromValues(Seq("one", "two", "three"))
  }

  "MetricBuilder.fromRow" should "build a metric from a Row" in {
    val path = makeTempFile("data", ".tab")
    Io.writeLines(path, Seq("name\tcount", "Foo Bar\t42", "Car Dog\t32"))
    val parser = DelimitedDataParser(path=path, delimiter='\t')
    val metrics = parser.zipWithIndex.map { case (row, lineNumber) =>
      builder.fromRow(row=row, headers=parser.headers, lineNumber=Some(lineNumber))
    }.toIndexedSeq
    metrics.length shouldBe 2
    metrics.head shouldBe MetricBuilderTestMetric(name="Foo Bar", count=42)
    metrics.last shouldBe MetricBuilderTestMetric(name="Car Dog", count=32)
  }

  it should "build a metric if the row is missing an optional field" in {
    val path = makeTempFile("data", ".tab")
    Io.writeLines(path, Seq("name\tnumber", "Foo Bar\t42", "Car Dog\t32"))
    val parser = DelimitedDataParser(path=path, delimiter='\t')
    val metrics = parser.zipWithIndex.map { case (row, lineNumber) =>
      builder.fromRow(row=row, headers=parser.headers, lineNumber=Some(lineNumber))
    }.toIndexedSeq
    metrics.length shouldBe 2
    metrics.head shouldBe MetricBuilderTestMetric(name="Foo Bar")
    metrics.last shouldBe MetricBuilderTestMetric(name="Car Dog")
  }

  it should "fail if the row is missing a required field" in {
    val path = makeTempFile("data", ".tab")
    Io.writeLines(path, Seq("desc\tcount", "Foo Bar\t42", "Car Dog\t32"))
    val parser = DelimitedDataParser(path=path, delimiter='\t')
    parser.zipWithIndex.foreach { case (row, lineNumber) =>
      a[MetricBuilderException] should be thrownBy builder.fromRow(row=row, headers=parser.headers, lineNumber=Some(lineNumber))
    }
  }
}
