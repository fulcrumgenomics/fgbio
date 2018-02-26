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


object SampleBarcodeMetric {
  def apply(barcodeName: String, libraryName: String, barcode: String): SampleBarcodeMetric = {
    new SampleBarcodeMetric(barcode_name=barcodeName, library_name=libraryName, barcode=barcode)
  }

  /** Computes values that are require the summary counts across multiple barcode metrics, such as certain fractions.
    *
    * @param barcodeToMetrics the map in which barcode metrics per sample are stored.
    * @param noMatchBarcode the barcode the unmatched reads.  This should stored in `barcodeToMetrics`.
    */
  def finalizeMetrics(barcodeToMetrics: Map[String, SampleBarcodeMetric],
                      noMatchBarcode: String): Unit = {
    val noMatchMetric = barcodeToMetrics(noMatchBarcode)

    var totalReads: Long           = 0
    var totalPfReads: Long         = 0
    var totalPfReadsAssigned: Long = 0

    barcodeToMetrics.foreach { case (barcode, metric) =>
      totalReads           += metric.reads
      totalPfReads         += metric.pf_reads
      totalPfReadsAssigned += metric.pf_reads
    }

    if (totalReads > 0) {
      noMatchMetric.pct_matches = noMatchMetric.reads / totalReads.toDouble
      var bestPctOfAllBarcodeMatches: Double = 0
      barcodeToMetrics.foreach { case (_, metric) =>
        val pctMatches =  metric.reads / totalReads.toDouble
        if (pctMatches > bestPctOfAllBarcodeMatches) {
          bestPctOfAllBarcodeMatches = pctMatches
        }
        metric.pct_matches = pctMatches
      }
      if (bestPctOfAllBarcodeMatches > 0) {
        noMatchMetric.ratio_this_barcode_to_best_barcode_pct = noMatchMetric.pct_matches / bestPctOfAllBarcodeMatches
        barcodeToMetrics.foreach { case (_, metric) =>
          metric.ratio_this_barcode_to_best_barcode_pct = metric.pct_matches / bestPctOfAllBarcodeMatches
        }
      }
    }
    if (totalPfReads > 0) {
      var bestPfPctOfAllBarcodeMatches: Double = 0
      barcodeToMetrics.foreach { case (_, metric) =>
        val pctPfMatches = metric.pf_reads / totalPfReads.toDouble
        if (pctPfMatches > bestPfPctOfAllBarcodeMatches) {
          bestPfPctOfAllBarcodeMatches = pctPfMatches
        }
        metric.pf_pct_matches = pctPfMatches
      }
      if (bestPfPctOfAllBarcodeMatches > 0) {
        noMatchMetric.pf_ratio_this_barcode_to_best_barcode_pct = noMatchMetric.pf_pct_matches / bestPfPctOfAllBarcodeMatches
        barcodeToMetrics.foreach { case (_, metric) =>
          metric.pf_ratio_this_barcode_to_best_barcode_pct = metric.pf_pct_matches / bestPfPctOfAllBarcodeMatches
        }
      }

    }
    if (totalPfReadsAssigned > 0) {
      val mean: Double = totalPfReadsAssigned.toDouble / barcodeToMetrics.values.size.toDouble
      barcodeToMetrics.foreach { case (barcode, metric) =>
        metric.pf_normalized_matches = metric.pf_reads / mean
      }
    }
  }
}

/**
  * Metrics for matching reads to sample barcodes.
  */
case class SampleBarcodeMetric(
  var barcode_name: String                                         = "",
  var library_name: String                                         = "",
  var barcode: String                                              = "",
  var reads: Metric.Count                                          = 0,
  var pf_reads: Metric.Count                                       = 0,
  var perfect_matches: Metric.Count                                = 0,
  var pf_perfect_matches: Metric.Count                             = 0,
  var one_mismatch_matches: Metric.Count                           = 0,
  var pf_one_mismatch_matches: Metric.Count                        = 0,
  var pct_matches: Metric.Proportion                               = 0d,
  var ratio_this_barcode_to_best_barcode_pct: Metric.Proportion    = 0d,
  var pf_pct_matches: Metric.Proportion                            = 0d,
  var pf_ratio_this_barcode_to_best_barcode_pct: Metric.Proportion = 0d,
  var pf_normalized_matches: Metric.Proportion                     = 0d
) extends Metric {

  /** Increments the counts for the metric. */
  def increment(numMismatches: Int, isPf: Boolean = true): Unit = {
    this.reads += 1
    if (isPf) this.pf_reads += 1
    if (numMismatches == 0) {
      this.perfect_matches += 1
      if (isPf) this.pf_perfect_matches += 1
    }
    else if (numMismatches == 1) {
      this.one_mismatch_matches += 1
      if (isPf) this.pf_one_mismatch_matches += 1
    }
  }
}
