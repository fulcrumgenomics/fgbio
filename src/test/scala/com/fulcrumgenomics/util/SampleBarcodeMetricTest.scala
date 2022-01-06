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
import scala.collection.mutable

/**
  * Tests for SampleBarcodeMetric.
  */
class SampleBarcodeMetricTest extends UnitSpec {

  "SampleBarcodeMetric.finalizeMetrics" should "fill out all metrics that require summary counts" in {
    val metricOne = new SampleBarcodeMetric("Alice", "Alice.lib", "AAAAAAA", 10, 5,   9,  4, 1, 1, total_number_of_bases = 10, q30_bases = 5, q20_bases = 10)
    val metricTwo = new SampleBarcodeMetric("Bob",   "Bob.lib",   "CCCCCCC", 20, 10, 19,  9, 1, 1)
    val noMatch   = new SampleBarcodeMetric("Eve",   "Eve.lib",   "NNNNNNN", 30, 20, 29, 19, 1, 1)
    val barcodeToMetrics = Seq(metricOne, metricTwo, noMatch).map { metric => (metric.barcode, metric) }.toMap
    SampleBarcodeMetric.finalizeMetrics(barcodeToMetrics, noMatch.barcode)

    metricOne.fraction_matches shouldBe 10 / 60d
    metricOne.ratio_this_barcode_to_best_barcode shouldBe 5 / 15d
    metricOne.pf_fraction_matches shouldBe 5 / 35d
    metricOne.pf_ratio_this_barcode_to_best_barcode shouldBe 5 / 20d
    metricOne.pf_normalized_matches shouldBe 3 / 7d +- 0.00001
    metricOne.frac_q30_bases shouldBe 1 / 2d

    metricTwo.fraction_matches shouldBe 20 / 60d
    metricTwo.ratio_this_barcode_to_best_barcode shouldBe 10 / 15d
    metricTwo.pf_fraction_matches shouldBe 10 / 35d
    metricTwo.pf_ratio_this_barcode_to_best_barcode shouldBe 10 / 20d
    metricTwo.pf_normalized_matches shouldBe 6 / 7d +- 0.00001

    noMatch.fraction_matches shouldBe 30 / 60d
    noMatch.ratio_this_barcode_to_best_barcode shouldBe 1d
    noMatch.pf_fraction_matches shouldBe 20 / 35d
    noMatch.pf_ratio_this_barcode_to_best_barcode shouldBe 1d
    noMatch.pf_normalized_matches shouldBe 12 / 7d +- 0.00001
  }

  "SampleBarcodeMetric.increment" should "not increment `pf_*` variables when the read doesn't pass QC" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 1, isPf = false, basesToAdd = 1, q30Bases = 1, q20Bases = 1, omitFailing = true)
    metric.templates shouldEqual 1
    metric.pf_templates shouldEqual 0
    metric.one_mismatch_matches shouldEqual 1
    metric.pf_one_mismatch_matches shouldEqual 0
  }

  it should "increment `pf_*` variables as well when the read passes QC" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 1, isPf = true, basesToAdd = 1, q30Bases =  1, q20Bases = 1, omitFailing = true)
    metric.templates shouldEqual 1
    metric.pf_templates shouldEqual 1
    metric.one_mismatch_matches shouldEqual 1
    metric.pf_one_mismatch_matches shouldEqual 1
  }

  it should "increment the total number of bases if the read passes QC, and --omit-failing-reads=false" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 0, isPf=true, basesToAdd = 1, q30Bases = 1, q20Bases = 1, omitFailing = false)
    metric.total_number_of_bases shouldBe 1
    metric.q20_bases shouldBe 1
    metric.q30_bases shouldBe 1
  }

  it should "increment the total number of bases if the read passes QC, and --omit-failing-reads=true" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 0, isPf=true, basesToAdd = 1, q30Bases = 1, q20Bases = 1, omitFailing = true)
    metric.total_number_of_bases shouldBe 1
    metric.q20_bases shouldBe 1
    metric.q30_bases shouldBe 1
  }

  it should "increment the total number of bases if the read does not pass QC and --omit-failing-reads=false" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 0, isPf=false, basesToAdd = 1, q30Bases = 1, q20Bases = 1, omitFailing = false)
    metric.total_number_of_bases shouldBe 1
    metric.q20_bases shouldBe 1
    metric.q30_bases shouldBe 1
  }

  it should "not increment the total number of bases if the read does not pass QC and --omit-failing-reads=true" in {
    val metric = new SampleBarcodeMetric()

    metric.increment(numMismatches = 0, isPf=false, basesToAdd = 1, q30Bases = 1, q20Bases = 1, omitFailing = true)
    metric.total_number_of_bases shouldBe 0
    metric.q20_bases shouldBe 0
    metric.q30_bases shouldBe 0
  }
}