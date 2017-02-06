/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.filtration

import com.fulcrumgenomics.bam.{Pileup, PileupEntry}
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.vcf.{VCFFilterHeaderLine, VCFHeaderLineType, VCFInfoHeaderLine}

/**
  * Trait for classes that can compute annotations on a somatic variant call and optionally
  * apply one or more filters to the calls.
  */
trait SomaticVariantFilter {
  /** The collection of VCF INFO header lines that the filter may reference. */
  val VcfInfoLines: Traversable[VCFInfoHeaderLine]

  /** The collection of VCF Filter header lines that the filter may reference. */
  val VcfFilterLines: Traversable[VCFFilterHeaderLine]

  /** Helper method to construct an INFO header line. */
  protected def vcfInfoLine(name: String, description: String, count: Int = 1, vcfType: VCFHeaderLineType = VCFHeaderLineType.Float) = {
    new VCFInfoHeaderLine(name, count, vcfType, description)
  }

  /** Helper method to construct a filter header line. */
  protected def vcfFilterLine(name: String, description: String) = new VCFFilterHeaderLine(name, description)

  /** Returns true if the filter can be applied to a genotype, false otherwise. */
  def appliesTo(gt: Genotype): Boolean

  /** Calculate the set of annotations for the filter, and return them as a Map. */
  def annotations(pileup: Pileup[PileupEntry], gt: Genotype): Map[String,Any]

  /** Given the set of annotations calculated by [[annotations()]] determine the set of filters
    * to be applied to the VCF record.
    */
  def filters(annotations: Map[String, Any]): Traversable[String]
}
