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

import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.testing.UnitSpec


case class MetricBuilderTestMetric(name: String, count: Long = 1) extends Metric

case class MetricBuilderNameMetric(first: String, second: String, age: Int) extends Metric

class MetricBuilderTest extends UnitSpec {
  private val builder = new MetricBuilder[MetricBuilderTestMetric]()

  "MetricBuilder.fromArgMap" should "build a metric from an argmap with all value specified" in {
    builder.fromArgMap(Map("name" -> "foo", "count" -> "2")) shouldBe MetricBuilderTestMetric(name="foo", count=2)
  }

  it should "build a metric from an argmap with only required values specified" in {
    builder.fromArgMap(Map("name" -> "foo")) shouldBe MetricBuilderTestMetric(name="foo")
  }

  it should "build a metric from a delimited line" in {
    builder.fromLine(line = "Foo Bar\t42", delim = "\t") shouldBe MetricBuilderTestMetric(name="Foo Bar", count=42)
    builder.fromLine(line = "Foo Bar,42", delim = ",") shouldBe MetricBuilderTestMetric(name="Foo Bar", count=42)

    val nameBuilder = new MetricBuilder[MetricBuilderNameMetric]()
    nameBuilder.fromLine(line = "Foo Bar 42", delim = " ") shouldBe MetricBuilderNameMetric(first="Foo", second="Bar", age=42)
  }

  it should "throw an FailureException when the # of values are incorrect or wrong type" in {
    val builder = new MetricBuilder[MetricBuilderNameMetric]()
    an[FailureException] should be thrownBy builder.fromLine(line="one\ttwo")
    an[FailureException] should be thrownBy builder.fromLine(line="one\ttwo\t3\tfour", lineNumber=Some(1))
    an[FailureException] should be thrownBy builder.fromLine(line="one\ttwo\tthree")
  }
}
