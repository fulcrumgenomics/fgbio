/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef._

/**
  * Map used to store genotypes and ensure that they are stored in the order they were read in the VCF. This class
  * is immutable.  To avoid copying of the gts array, once passed to this class no other references to the gts
  * should be retained.
  *
  * @param gts the genotypes for the samples
  * @param sampleIndex a mapping of sample name to index within the genotypes
  */
class GenotypeMap private[api] (private val gts: Array[Genotype],
                                private val sampleIndex: Map[String, Int]) extends Map[String, Genotype] {

  override def removed(key: String): Map[String, Genotype] = throw new UnsupportedOperationException
  override def updated[V1 >: Genotype](key: String, value: V1): Map[String, V1] = throw new UnsupportedOperationException
  override def get(key: String): Option[Genotype] = sampleIndex.get(key).map(i => gts(i))
  override def apply(key: String): Genotype = gts(sampleIndex(key))
  override def iterator: Iterator[(String, Genotype)] = gts.iterator.map(gt => (gt.sample, gt))
}
