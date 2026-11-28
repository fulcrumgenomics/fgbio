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

object JointVariantIterator {
  /** Creates an iterator over the union of sites across the given iterators, ordering variants by the index of their
    * contig in the sequence dictionary and then by their start position. */
  def apply(iters: Seq[Iterator[Variant]], dict: SequenceDictionary): JointVariantIterator = {
    new JointVariantIterator(iters=iters, dict=dict)
  }
}

/**
  * Iterates over multiple variant iterators such that we return a list of variants for the union of sites
  * across the iterators.  Variants are ordered by the index of their contig in the sequence dictionary and then by
  * their start position, and variants at the same contig and start position are returned together.
  */
class JointVariantIterator private(iters: Seq[Iterator[Variant]], dict: SequenceDictionary)
extends Iterator[Seq[Option[Variant]]] {

  if (iters.isEmpty) throw new IllegalArgumentException("No iterators given")

  private val iterators = iters.map(_.buffered).toIndexedSeq

  def hasNext: Boolean = iterators.exists(_.nonEmpty)

  def next(): Seq[Option[Variant]] = {
    if (!hasNext) throw new NoSuchElementException("Called next when hasNext is false")
    // TODO: could use a TreeSet to store the iterators, examine the head of each iterator, then pop the iterator with the min,
    // and add that iterator back in.
    val minVariant     = iterators.iterator.filter(_.nonEmpty).map(_.head).reduce((min, v) => if (compare(v, min) < 0) v else min)
    val minContigIndex = dict(minVariant.chrom).index
    iterators.map { iter =>
      if (iter.nonEmpty && iter.head.pos == minVariant.pos && dict(iter.head.chrom).index == minContigIndex) Some(iter.next())
      else None
    }
  }

  /** Compares two variants by the index of their contig in the sequence dictionary and then by their start position. */
  private def compare(left: Variant, right: Variant): Int = {
    val byContig = Integer.compare(dict(left.chrom).index, dict(right.chrom).index)
    if (byContig != 0) byContig else Integer.compare(left.pos, right.pos)
  }
}
