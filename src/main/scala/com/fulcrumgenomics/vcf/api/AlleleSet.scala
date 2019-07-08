/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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
import com.fulcrumgenomics.vcf.api.Allele.SimpleAllele

case class AlleleSet(ref: SimpleAllele, alts: IndexedSeq[Allele]) extends Iterable[Allele] {
  override def iterator: Iterator[Allele] = Iterator(ref) ++ alts.iterator

  def apply(index: Int): Allele = if (index == 0) ref else alts(index-1)

  def indexOf(a: Allele): Int = if (a == ref) 0 else alts.indexOf(a) + 1
}

object AlleleSet {
  def apply(ref: Allele, alts: Traversable[Allele]): AlleleSet = ref match {
    case r: SimpleAllele => AlleleSet(r, alts.toIndexedSeq)
    case _ => throw new IllegalArgumentException(s"Cannot have a non-simple ref allele: $ref")
  }
}
