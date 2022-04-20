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

import com.fulcrumgenomics.testing.UnitSpec

case class MetricSorterTestMetric(name: String, count: Long) extends Metric with Ordered[MetricSorterTestMetric] {
  override def compare(that: MetricSorterTestMetric): Int = {
    var retval = this.count.compare(that.count)
    if (retval == 0) retval = this.name.compare(that.name)
    retval
  }
}

class MetricSorterTest extends UnitSpec {

  private val metrics = IndexedSeq(
    MetricSorterTestMetric(name="foo", count=10),
    MetricSorterTestMetric(name="foo", count=1),
    MetricSorterTestMetric(name="bar", count=1),
    MetricSorterTestMetric(name="foo", count=5),
    MetricSorterTestMetric(name="roger", count=2),
    MetricSorterTestMetric(name="nadal", count=2),
  )

  private val metricsSorted = metrics.sortBy(m => (m.count, m.name))

  private case class Key()

  "MetricSorter" should "sort metrics in memory" in {
    val sorter = new MetricSorter[MetricSorterTestMetric, MetricSorterTestMetric](
      maxObjectsInRam = 10000,
      keyfunc         = identity
    )
    sorter ++= metrics
    sorter.iterator.toSeq should contain theSameElementsInOrderAs metricsSorted
  }

  it should "sort metrics after spilling to disk" in {
    val sorter = new MetricSorter[MetricSorterTestMetric, MetricSorterTestMetric](
      maxObjectsInRam = 2,
      keyfunc         = identity
    )
    sorter ++= metrics
    sorter.iterator.toSeq should contain theSameElementsInOrderAs metricsSorted
  }
}
