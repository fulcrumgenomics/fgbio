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
import htsjdk.variant.variantcontext.{VariantContext, VariantContextComparator}

object JointVariantContextIterator {
  def apply(iters: Seq[Iterator[VariantContext]],
            dict: SequenceDictionary
           ): JointVariantContextIterator = {
    new JointVariantContextIterator(
      iters=iters,
      dictOrComp = Left(dict)
    )
  }

  def apply(iters: Seq[Iterator[VariantContext]],
            comp: VariantContextComparator
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
class JointVariantContextIterator private(iters: Seq[Iterator[VariantContext]],
                                          dictOrComp: Either[SequenceDictionary, VariantContextComparator]
                                         )
extends Iterator[Seq[Option[VariantContext]]] {
  import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary

  if (iters.isEmpty) throw new IllegalArgumentException("No iterators given")

  private val iterators = iters.map(_.buffered)
  private val comparator = dictOrComp match {
    case Left(dict)  => new VariantContextComparator(dict.asSam)
    case Right(comp) => comp
  }

  def hasNext: Boolean = iterators.exists(_.nonEmpty)

  def next(): Seq[Option[VariantContext]] = {
    val minCtx = iterators.filter(_.nonEmpty).map(_.head).sortWith {
      case (left: VariantContext, right: VariantContext) => comparator.compare(left, right) < 0
    }.head
    // TODO: could use a TreeSet to store the iterators, examine the head of each iterator, then pop the iterator with the min,
    // and add that iterator back in.
    iterators.zipWithIndex.map { case(iter, idx) =>
      if (iter.isEmpty || this.comparator.compare(minCtx, iter.head) != 0) None
      else Some(iter.next())
    }
  }
}
