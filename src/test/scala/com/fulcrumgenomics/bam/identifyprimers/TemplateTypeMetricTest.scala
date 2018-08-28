/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.identifyprimers

import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.testing.UnitSpec

class TemplateTypeMetricTest extends UnitSpec {

  "TemplateTypeMetric" should "be built from primer matches and be collatable" in {
    // build a couple
    val primer = new Primer("pair1", "primer1+", "GATTACA", "chr1", 1, 7, true)
    val metrics = Seq(
      TemplateTypeMetric(TemplateType.MappedPair, true, PrimerPairMatchType.Canonical, Some(LocationBasedPrimerMatch(primer, 1)), None),
      TemplateTypeMetric(TemplateType.MappedPair, true, PrimerPairMatchType.Canonical, Some(LocationBasedPrimerMatch(primer, 2)), None),
      TemplateTypeMetric(TemplateType.MappedFragment, false, PrimerPairMatchType.Single, Some(LocationBasedPrimerMatch(primer, 3)), None)
    )

    // collate
    val counter = new SimpleCounter[TemplateTypeMetric]()
    metrics.foreach { m => counter.count(m) }
    val newMetrics = TemplateTypeMetric.metricsFrom(counter)
    newMetrics should have length 2
    newMetrics.head shouldBe metrics.head.copy(count = 2, frac = 2/3.toDouble)
    newMetrics.last shouldBe metrics.last.copy(count = 1, frac = 1/3.toDouble)
  }
}
