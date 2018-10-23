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

import com.fulcrumgenomics.commons.util.{NumericCounter, SimpleCounter}
import com.fulcrumgenomics.util.Metric


private[identifyprimers] object IdentifyPrimersMetric {

  /** Creates a [[IdentifyPrimersMetric]] from [[TemplateTypeMetric]]s. */
  def apply(templateTypeMetrics: Seq[TemplateTypeMetric]): IdentifyPrimersMetric = {
    import TemplateType._

    val readTypeCounter        = new SimpleCounter[TemplateType]()
    val primerMatchCounter     = new SimpleCounter[String]()
    val endsWithMatchesCounter = new NumericCounter[Int]()

    val fr_pairs = templateTypeMetrics.map { templateTypeMetric =>
      val count = templateTypeMetric.count

      // Add the count of the type of template
      readTypeCounter.count(templateTypeMetric.template_type, count)

      // Add the count of types of primer matches
      val primer_match_types = Seq(templateTypeMetric.r1_primer_match_type, templateTypeMetric.r2_primer_match_type)
      primer_match_types.foreach { primerMatchType =>
        primerMatchCounter.count(primerMatchType, count)
      }

      // Add the count of # of primer matches per template
      primer_match_types.map(_ == PrimerMatch.NoPrimerMatchName) match {
        case Seq(true, true)   => endsWithMatchesCounter.count(0, count)
        case Seq(false, false) => endsWithMatchesCounter.count(2, count)
        case Seq(_, _)         => endsWithMatchesCounter.count(1, count)
      }

      // return this to count the # of FR pairs
      if (templateTypeMetric.is_fr) templateTypeMetric.count else 0
    }.sum

    val mapped_pairs       = readTypeCounter.countOf(MappedPair)
    val unpaired           = readTypeCounter.countOf(Unpaired)
    val unmapped_pairs     = readTypeCounter.countOf(UnmappedPair)
    val mapped_fragments   = readTypeCounter.countOf(MappedFragment)
    val unmapped_fragments = readTypeCounter.countOf(UnmappedFragment)

    new IdentifyPrimersMetric(
      templates                 = readTypeCounter.total,
      pairs                     = mapped_pairs + unpaired + unmapped_pairs,
      fragments                 = mapped_fragments + unmapped_fragments,
      mapped_pairs              = mapped_pairs,
      unpaired                  = unpaired,
      unmapped_pairs            = unmapped_pairs,
      fr_pairs                  = fr_pairs,
      mapped_fragments          = mapped_fragments,
      unmapped_fragments        = unmapped_fragments,
      paired_matches            = endsWithMatchesCounter.countOf(2),
      unpaired_matches          = endsWithMatchesCounter.countOf(1),
      no_paired_matches         = endsWithMatchesCounter.countOf(0),
      match_attempts            = primerMatchCounter.total,
      location                  = primerMatchCounter.countOf(PrimerMatch.toName[LocationBasedPrimerMatch]),
      ungapped                  = primerMatchCounter.countOf(PrimerMatch.toName[UngappedAlignmentPrimerMatch]),
      gapped          = primerMatchCounter.countOf(PrimerMatch.toName[GappedAlignmentPrimerMatch]),
      no_match                  = primerMatchCounter.countOf(PrimerMatch.NoPrimerMatchName)
    )
  }
}

/** Provides summary information on the counts of templates, [[TemplateType]]s, and [[PrimerMatch]] types.
  *
  * @param templates the number of templates examined.
  * @param pairs the number of reads pairs examined.
  * @param fragments the number of fragment raeds examined.
  * @param mapped_pairs the number of mapped pairs examined.
  * @param unpaired the number of pairs with exactly one end mapped examined.
  * @param unmapped_pairs the number of unmapped pairs examined.
  * @param fr_pairs the number of mapped pairs that are in FR orientation: `( 5' --F-->  <--R-- 5' )`  - aka. innie
  * @param mapped_fragments the number of mapped fragment reads examined.
  * @param unmapped_fragments the number of unmapped fragment reads examined.
  * @param paired_matches the number of templates with both ends of a pair with matches.
  * @param unpaired_matches tthe number of templates with only one end of a pair with a match.
  * @param no_paired_matches the number of templates with neither end of a pair with a match.
  * @param match_attempts the number of attempts at primer match (i.e. `(2 * pairs) + fragments`).
  * @param location the number of primer matches made based on mapping location.
  * @param ungapped the number of primer matches made based on ungapped (mismatch-only) alignment.
  * @param gapped the number of primer matches made based on gapped alignment.
  * @param no_match the number of reads that do no match a primer.
  */
case class IdentifyPrimersMetric
(templates: Long = 0,
 pairs: Long = 0,
 fragments: Long = 0,
 mapped_pairs: Long = 0,
 unpaired: Long = 0,
 unmapped_pairs: Long = 0,
 fr_pairs: Long = 0,
 mapped_fragments: Long = 0,
 unmapped_fragments: Long = 0,
 paired_matches: Long = 0,
 unpaired_matches: Long = 0,
 no_paired_matches: Long = 0,
 match_attempts: Long = 0,
 location: Long = 0,
 ungapped: Long = 0,
 gapped: Long = 0,
 no_match: Long = 0
) extends Metric