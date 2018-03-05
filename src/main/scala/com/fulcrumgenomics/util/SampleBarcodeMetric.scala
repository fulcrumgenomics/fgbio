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
    * @param noMatchBarcode the barcode for the unmatched templates.  This should stored in `barcodeToMetrics`.
    */
  def finalizeMetrics(barcodeToMetrics: Map[String, SampleBarcodeMetric],
                      noMatchBarcode: String): Unit = {
    val noMatchMetric = barcodeToMetrics(noMatchBarcode)

    var totalReads: Long           = 0
    var totalPfReads: Long         = 0
    var totalPfReadsAssigned: Long = 0

    barcodeToMetrics.foreach { case (barcode, metric) =>
      totalReads           += metric.templates
      totalPfReads         += metric.pf_templates
      totalPfReadsAssigned += metric.pf_templates
    }

    if (totalReads > 0) {
      noMatchMetric.fraction_matches = noMatchMetric.templates / totalReads.toDouble
      var bestPctOfAllBarcodeMatches: Double = 0
      barcodeToMetrics.foreach { case (_, metric) =>
        val fracMatches =  metric.templates / totalReads.toDouble
        if (fracMatches > bestPctOfAllBarcodeMatches) {
          bestPctOfAllBarcodeMatches = fracMatches
        }
        metric.fraction_matches = fracMatches
      }
      if (bestPctOfAllBarcodeMatches > 0) {
        noMatchMetric.ratio_this_barcode_to_best_barcode = noMatchMetric.fraction_matches / bestPctOfAllBarcodeMatches
        barcodeToMetrics.foreach { case (_, metric) =>
          metric.ratio_this_barcode_to_best_barcode = metric.fraction_matches / bestPctOfAllBarcodeMatches
        }
      }
    }
    if (totalPfReads > 0) {
      var bestPfPctOfAllBarcodeMatches: Double = 0
      barcodeToMetrics.foreach { case (_, metric) =>
        val fracPfMatches = metric.pf_templates / totalPfReads.toDouble
        if (fracPfMatches > bestPfPctOfAllBarcodeMatches) {
          bestPfPctOfAllBarcodeMatches = fracPfMatches
        }
        metric.pf_fraction_matches = fracPfMatches
      }
      if (bestPfPctOfAllBarcodeMatches > 0) {
        noMatchMetric.pf_ratio_this_barcode_to_best_barcode = noMatchMetric.pf_fraction_matches / bestPfPctOfAllBarcodeMatches
        barcodeToMetrics.foreach { case (_, metric) =>
          metric.pf_ratio_this_barcode_to_best_barcode = metric.pf_fraction_matches / bestPfPctOfAllBarcodeMatches
        }
      }

    }
    if (totalPfReadsAssigned > 0) {
      val mean: Double = totalPfReadsAssigned.toDouble / barcodeToMetrics.values.size.toDouble
      barcodeToMetrics.foreach { case (barcode, metric) =>
        metric.pf_normalized_matches = metric.pf_templates / mean
      }
    }
  }
}

/**
  * Metrics for matching templates to sample barcodes primarily used in [[com.fulcrumgenomics.fastq.DemuxFastqs]].
  *
  * The number of templates will match the number of reads for an Illumina single-end sequencing run, while the number
  * of templates will be half the number of reads for an Illumina paired-end sequencing run (i.e. R1 & R2 observe the
  * same template).
  *
  * @param barcode_name the name for the sample barcode, typically the sample name from the SampleSheet.
  * @param library_name the name of the library, typically the library identifier from the SampleSheet.
  * @param barcode the sample barcode bases.  Dual index barcodes will have two sample barcode sequences delimited by a
  *                dash.
  * @param templates the total number of templates matching the given barcode.
  * @param pf_templates the total number of pass-filter templates matching the given barcode.
  * @param perfect_matches the number of templates that match perfectly the given barcode.
  * @param pf_perfect_matches the number of pass-filter templates that match perfectly the given barcode.
  * @param one_mismatch_matches the number of pass-filter templates that match the given barcode with exactly one
  *                             mismatch.
  * @param pf_one_mismatch_matches the number of pass-filter templates that match the given barcode with exactly
  *                                one mismatch.
  * @param fraction_matches the fraction of all templates that match the given barcode.
  * @param ratio_this_barcode_to_best_barcode the rate of all templates matching this barcode to all template
  *                                               reads matching the most prevalent barcode. For the most prevalent
  *                                               barcode this will be 1, for all others it will be less than 1 (except
  *                                               for the possible exception of when there are more unmatched templates
  *                                               than for any other barcode, in which case the value may be arbitrarily
  *                                               large).  One over the lowest number in this column gives you the
  *                                               fold-difference in representation between barcodes.
  * @param pf_fraction_matches the fraction of all pass-filter templates that match the given barcode.
  * @param pf_ratio_this_barcode_to_best_barcode the rate of all pass-filter templates matching this barcode to
  *                                                  all templates matching the most prevalent barcode. For the
  *                                                  most prevalent barcode this will be 1, for all others it will be
  *                                                  less than 1 (except for the possible exception of when there are
  *                                                  more unmatched templates than for any other barcode, in which
  *                                                  case the value may be arbitrarily large).  One over the lowest
  *                                                  number in this column gives you the fold-difference in
  *                                                  representation between barcodes.
  * @param pf_normalized_matches The "normalized" matches to each barcode. This is calculated as the number of
  *                              pass-filter templates matching this barcode over the mean of all pass-filter
  *                              templates matching any barcode (excluding unmatched). If all barcodes are
  *                              represented equally this will be
  *                              1.
  */
case class SampleBarcodeMetric
( var barcode_name: String                                     = "",
  var library_name: String                                     = "",
  var barcode: String                                          = "",
  var templates: Metric.Count                                  = 0,
  var pf_templates: Metric.Count                               = 0,
  var perfect_matches: Metric.Count                            = 0,
  var pf_perfect_matches: Metric.Count                         = 0,
  var one_mismatch_matches: Metric.Count                       = 0,
  var pf_one_mismatch_matches: Metric.Count                    = 0,
  var fraction_matches: Metric.Proportion                      = 0d,
  var ratio_this_barcode_to_best_barcode: Metric.Proportion    = 0d,
  var pf_fraction_matches: Metric.Proportion                   = 0d,
  var pf_ratio_this_barcode_to_best_barcode: Metric.Proportion = 0d,
  var pf_normalized_matches: Metric.Proportion                 = 0d
) extends Metric {

  /** Increments the counts for the metric. */
  def increment(numMismatches: Int, isPf: Boolean = true): Unit = {
    this.templates += 1
    if (isPf) this.pf_templates += 1
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
