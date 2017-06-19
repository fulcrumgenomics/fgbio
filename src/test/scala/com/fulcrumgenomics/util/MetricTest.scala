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

import java.io.StringWriter
import java.nio.file.Path

import com.fulcrumgenomics.testing.UnitSpec
import org.scalatest.OptionValues

private case class TestMetric(foo: String, bar: Int, car: String = "default") extends Metric

private case class TestMetricWithOption(foo: String, bar: Int, option: Option[String]) extends Metric

private case class TestMetricWithIntOption(foo: String, bar: Int, option: Option[Int]) extends Metric


/**
  * Tests for Metric.
  */
class MetricTest extends UnitSpec with OptionValues {

  "Metric.header" should "return the header names in order" in {
    val testMetric = TestMetric(foo="fooValue", bar=1)
    testMetric.names should contain theSameElementsInOrderAs Seq("foo", "bar", "car")
    Metric.names[TestMetric] should contain theSameElementsInOrderAs Seq("foo", "bar", "car")
  }

  "Metric.values" should "return the values in order" in {
    val testMetric = TestMetric(foo="fooValue", bar=1)
    testMetric.values should contain theSameElementsInOrderAs Seq("fooValue", "1", "default")
  }

  "Metric.build" should "build a metric from the given arguments" in {
    val testMetric = Metric.build[TestMetric](Seq(("foo", "fooValue"), ("bar", "1"), ("car", "2")))
    testMetric.foo shouldBe "fooValue"
    testMetric.bar shouldBe 1
    testMetric.car shouldBe "2"

  }

  "Metric.build" should "build a metric when omitting a field with a default value" in {
    val testMetric = Metric.build[TestMetric](Seq(("foo", "fooValue"), ("bar", "1")))
    testMetric.foo shouldBe "fooValue"
    testMetric.bar shouldBe 1
    testMetric.car shouldBe "default"
  }

  it should "fail when an argument is not given and has no default value" in {
    an[Exception] should be thrownBy Metric.build[TestMetric](Seq(("foo", "fooValue"))) // bar has no default
  }

  it should "fail when an unknown argument is given" in {
    an[Exception] should be thrownBy Metric.build[TestMetric](Seq(("foo", "fooValue"), ("doh", "1"), ("bar", "1")))
  }

  it should "fail when an argument cannot be built from the given string" in {
    an[Exception] should be thrownBy Metric.build[TestMetric](Seq(("foo", "fooValue"), ("bar", "bar"))) // bar is an int
  }

  "Metrics.read" should "read metrics from a sequence of lines" in {
    val lines   = Seq("foo\tbar", "fooValue1\t1", "fooValue2\t2")
    val metrics = Metric.read[TestMetric](lines.iterator)
    metrics.zipWithIndex.foreach { case (metric, idx) =>
        metric.foo shouldBe s"fooValue${idx+1}"
        metric.bar shouldBe (idx+1)
    }
  }

  it should "read metrics from a sequence of lines with the header in any order" in {
    val lines   = Seq("bar\tfoo", "1\tfoo", "2\tbar")
    val metrics = Metric.read[TestMetric](lines.iterator)
    metrics.head.foo shouldBe "foo"
    metrics.head.bar shouldBe 1
    metrics.head.car shouldBe "default"
    metrics.last.foo shouldBe "bar"
    metrics.last.bar shouldBe 2
    metrics.last.car shouldBe "default"
  }

  it should "throw an exception when a metric line has a different number of values than the header" in {
    val lines   = Seq("foo\tbar", "fooValue1\t1", "fooValue2")
    an[Exception] should be thrownBy Metric.read[TestMetric](lines.iterator).toList
  }

  // Various write methods to test with a writer
  private val writeWithWriters = Seq(
    ((writer: StringWriter, metrics: Seq[TestMetric]) => Metric.write[TestMetric](writer, metrics:_*), "write[T <: Metric](writer: Writer, metric: T*)"),
    ((writer: StringWriter, metrics: Seq[TestMetric]) => Metric.write[TestMetric](writer, metrics), "write[T <: Metric](writer: Writer, metrics: TraversableOnce[T])")
  )

  writeWithWriters.foreach { case (writeMethod, desc) =>
    "Metrics.write" should s"write metrics in a sequence of lines [$desc]" in {
      val metrics = Seq(TestMetric(foo = "fooValue1", bar = 1), TestMetric(foo = "fooValue2", bar = 2))
      val writer  = new StringWriter()
      writeMethod(writer, metrics)
      val lines   = writer.toString
      lines shouldBe "foo\tbar\tcar\nfooValue1\t1\tdefault\nfooValue2\t2\tdefault\n"
    }
  }

  // Various write methods to test with a path
  private val writeWithPaths = Seq(
    ((path: Path, metrics: Seq[TestMetric]) => Metric.write[TestMetric](path, metrics:_*), "write[T <: Metric](path: Path, metric: T*)"),
    ((path: Path, metrics: Seq[TestMetric]) => Metric.write[TestMetric](path, metrics), "write[T <: Metric](path: Path, metric: TraversableOnce[T])")
  )

  writeWithPaths.foreach { case (writeMethod, desc) =>
    "Metrics.write" should s"write metrics in a sequence of lines [$desc]" in {
      val output      = makeTempFile("MetricTest", ".txt")
      val metrics = Seq(TestMetric(foo="fooValue1", bar=1), TestMetric(foo="fooValue2", bar=2))
      writeMethod(output, metrics)
      val lines   = Io.readLines(output).mkString("\n")
      lines shouldBe "foo\tbar\tcar\nfooValue1\t1\tdefault\nfooValue2\t2\tdefault"
    }
  }

  "Metrics" should "write to and read from a file" in {
    val metrics     = Seq(TestMetric(foo="fooValue1", bar=1), TestMetric(foo="fooValue2", bar=2))
    val output      = makeTempFile("MetricTest", ".txt")
    Metric.write[TestMetric](output, metrics:_*)
    val readMetrics = Metric.read[TestMetric](output)
    readMetrics.length shouldBe 2
    readMetrics.zip(metrics).foreach { case (in, out) => in.toString shouldBe out.toString }
  }

  it should "write and read metrics with an option" in {
    val metricWithSome = TestMetricWithOption(foo="fooValue1", bar=1, option=Some("option1"))
    val metricWithNone = TestMetricWithOption(foo="fooValue2", bar=2, option=None)
    val metrics = Seq(metricWithSome, metricWithNone)
    val output      = makeTempFile("MetricTest", ".txt")
    Metric.write(path=output, metrics:_*)
    val readMetrics = Metric.read[TestMetricWithOption](path=output)
    readMetrics.length shouldBe 2
    readMetrics.head shouldBe metricWithSome
    readMetrics.last shouldBe metricWithNone
  }

  it should "write and read metrics with an integer option" in {
    val metricWithSome = TestMetricWithIntOption(foo="fooValue1", bar=1, option=Some(1))
    val metricWithNone = TestMetricWithIntOption(foo="fooValue2", bar=2, option=None)
    val metrics = Seq(metricWithSome, metricWithNone)
    val output      = makeTempFile("MetricTest", ".txt")
    Metric.write(path=output, metrics:_*)
    val readMetrics = Metric.read[TestMetricWithIntOption](path=output)
    readMetrics.length shouldBe 2
    readMetrics.head shouldBe metricWithSome
    readMetrics.last shouldBe metricWithNone
  }

  it should "get the field by name" in {
    val testMetric = TestMetric(foo="fooValue", bar=1)

    testMetric("foo") shouldBe "fooValue"
    testMetric("bar") shouldBe "1"
    testMetric("car") shouldBe "default"
    an[NoSuchElementException] should be thrownBy testMetric("dar")

    testMetric.get("foo").value shouldBe "fooValue"
    testMetric.get("bar").value shouldBe "1"
    testMetric.get("car").value shouldBe "default"
    testMetric.get("dar") shouldBe 'empty
  }

  it should "be iterable over tuples of names and values" in {
    val testMetric = TestMetric(foo="fooValue", bar=1)
    testMetric.iterator.toSeq should contain theSameElementsInOrderAs Seq(("foo", "fooValue"), ("bar", "1"), ("car", "default"))
  }
}
