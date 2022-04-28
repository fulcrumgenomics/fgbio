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
 */

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.vcf.api.Variant

object JointVariantContextIterator {
  def apply(iters: Seq[Iterator[Variant]],
            dict: SequenceDictionary
           ): JointVariantContextIterator = {
    new JointVariantContextIterator(
      iters=iters,
      dictOrComp = Left(dict)
    )
  }

  def apply(iters: Seq[Iterator[Variant]],
            comp: VariantComparator
           ): JointVariantContextIterator = {
    new JointVariantContextIterator(
      iters=iters,
      dictOrComp = Right(comp)
    )
  }
}

/**
  * Iterates over multiple variant context iterators such that we return a list of contexts for the union of sites
  * across the iterators.  If samples is given, we subset each variant context to just that sample.
  */
class JointVariantContextIterator private(iters: Seq[Iterator[Variant]],
                                          dictOrComp: Either[SequenceDictionary, VariantComparator]
                                         )
extends Iterator[Seq[Option[Variant]]] {

  if (iters.isEmpty) throw new IllegalArgumentException("No iterators given")

  private val iterators = iters.map(_.buffered)
  private val comparator = dictOrComp match {
    case Left(dict)  => VariantComparator(dict)
    case Right(comp) => comp
  }

  def hasNext: Boolean = iterators.exists(_.nonEmpty)

  def next(): Seq[Option[Variant]] = {
    val minCtx = iterators.filter(_.nonEmpty).map(_.head).sortWith {
      // this just checks that the variants are the the same position? Shouldn't be difficult to replace.
      case (left: Variant, right: Variant) => comparator.compare(left, right) < 0
    }.head
    // TODO: could use a TreeSet to store the iterators, examine the head of each iterator, then pop the iterator with the min,
    // and add that iterator back in.
    iterators.zipWithIndex.map { case(iter, idx) =>
      if (iter.isEmpty || this.comparator.compare(minCtx, iter.head) != 0) None
      else Some(iter.next())
    }
  }
}

private object VariantComparator {
  def apply(dict: SequenceDictionary): VariantComparator = {
    new VariantComparator(dict)
  }
}

/** A class for comparing Variants using a sequence dictionary */
private class VariantComparator(dict: SequenceDictionary) {
  /** Function for comparing two variants. Returns negative if left < right, and positive if right > left
    * To mimic the VariantContextComparator, throws an exception if the contig isn't found. */
  def compare(left: Variant, right: Variant): Int = {
    val idx1 = this.dict(name=left.chrom).index
    val idx2 = this.dict(name=right.chrom).index
    if (idx1 - idx2 == 0) left.pos - right.pos
    else idx1 - idx2
  }
}
