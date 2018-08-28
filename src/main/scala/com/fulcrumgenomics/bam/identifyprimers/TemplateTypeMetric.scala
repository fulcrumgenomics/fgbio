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
import com.fulcrumgenomics.util.Metric

private[identifyprimers] object TemplateTypeMetric {

  /** Builds a [[TemplateTypeMetric]]. */
  def apply(readType: TemplateType,
            is_fr: Boolean,
            primerPairMatchType: PrimerPairMatchType,
            forwardPrimerMatchType: Option[PrimerMatch],
            reversePrimerMatchType: Option[PrimerMatch]
           ): TemplateTypeMetric = {
    TemplateTypeMetric(
      template_type          = readType,
      is_fr                  = is_fr,
      primer_pair_match_type = primerPairMatchType,
      r1_primer_match_type   = PrimerMatch.toName(forwardPrimerMatchType),
      r2_primer_match_type   = PrimerMatch.toName(reversePrimerMatchType)
    )
  }

  /** Creates a sequence of [[TemplateTypeMetric]] from a counter.  The metrics stored in the counter must each have
    * `metric.count == 0`.  The returned sequence will be ordered by most frequent template type (i.e. highest count).
    */
  def metricsFrom(counter: SimpleCounter[TemplateTypeMetric]): Seq[TemplateTypeMetric] = {
    val total   = counter.total.toDouble
    counter.map { case (metric, count) =>
      require(metric.count == 0)
      metric.copy(count = count, frac = count / total)
    }.toSeq.sortBy(-_.count)
  }
}


/** Stores summary counts for how many templates of type [[template_type]] have the given same value for
  * [[primer_pair_match_type]], [[r1_primer_match_type]], and [[r2_primer_match_type]].
  *
  * @param template_type          the [[TemplateType]].
  * @param is_fr                  true if a mapped read pair is in FR orientation: `( 5' --F-->  <--R-- 5' )`  - aka. innie
  * @param primer_pair_match_type the [[PrimerPairMatchType]].
  * @param r1_primer_match_type   [[PrimerMatch]]for matches read one (or fragment reads), otherwise [[None]].
  * @param r2_primer_match_type   [[PrimerMatch]]for matches read two, otherwise [[None]].
  * @param count                  the number of templates of type [[template_type]] having the same value for [[primer_pair_match_type]],
  *                               [[r1_primer_match_type]], and [[r2_primer_match_type]].
  * @param frac                   the fraction of all templates having the same value for [[primer_pair_match_type]],
  *                               [[r1_primer_match_type]], and [[r2_primer_match_type]].
  */
private[identifyprimers] case class TemplateTypeMetric
(template_type: TemplateType,
 is_fr: Boolean,
 primer_pair_match_type: PrimerPairMatchType,
 r1_primer_match_type: String,
 r2_primer_match_type: String,
 count: Long = 0,
 frac: Double = 0d
) extends Metric

