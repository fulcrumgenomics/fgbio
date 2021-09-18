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

/**
  * Simple class to encapsulate a set of alleles where one allele represents the allele
  * that is present in the reference genome, and the remainder represent alternative alleles.
  *
  * @param ref the reference allele
  * @param alts the ordered list of alternate alleles
  */
case class AlleleSet(ref: SimpleAllele, alts: IndexedSeq[Allele]) extends Iterable[Allele] {
  /** Provides an iterator over all the alleles. */
  override def iterator: Iterator[Allele] = Iterator(ref) ++ alts.iterator

  /** Retrieves the allele at `index`.  The reference allele is given index 0, and alternate
    * alleles are indexed starting at 1.  Requesting an index  `< 1` or `> size-1` will result
    *  in a [[IndexOutOfBoundsException]].
    *
    * @param index The index of the allele to return
    * @return the allele present at the index
    */
  def apply(index: Int): Allele = if (index == 0) ref else alts(index-1)

  /** Retrieves the allele at `index`.  The reference allele is given index 0, and alternate
    * alleles are indexed starting at 1.
    *
    * @param index The index of the allele to return
    * @return the allele present at the index, none if not found (i.e. index is  `< 1` or `> size-1`)
    */
  def get(index: Int): Option[Allele] = {
    if (index < 0 || index > size - 1) None else Some(this.apply(index))
  }

  /** Returns the position of the given allele within the allele set. */
  def indexOf(a: Allele): Int = if (a == ref) 0 else alts.indexOf(a) + 1

  /** The number of alleles, including reference, in the set. */
  override def size: Int = alts.size + 1
}

object AlleleSet {
  private val NoAlts: IndexedSeq[Allele] = IndexedSeq.empty

  /**
    * Generates an AlleleSet from a reference allele and zero or more alternative alleles.
    *
    * @param ref the reference allele; if this is not a [[SimpleAllele]] an [[IllegalArgumentException]]
    *            will be thrown. The parameter type is [[Allele]] to avoid callers having to cast.
    * @param alts zero or more alternative alleles
    */
  def apply(ref: Allele, alts: Iterable[Allele]): AlleleSet = ref match {
    case r: SimpleAllele => AlleleSet(r, alts.toIndexedSeq)
    case _ => throw new IllegalArgumentException(s"Cannot have a non-simple ref allele: $ref")
  }

  /**
    * Generates an AlleleSet from a reference allele and zero or more alternative alleles.
    *
    * @param ref the reference allele; if this is not a [[SimpleAllele]] an [[IllegalArgumentException]]
    *            will be thrown. The parameter type is [[Allele]] to avoid callers having to cast.
    * @param alts zero or more alternative alleles
    */
  def apply(ref: Allele, alts: Allele*): AlleleSet = apply(ref, alts)

  /**
    * Generates an AlleleSet from one or more String alleles. The first allele will be the reference
    * allele and must be composed of regular bases.
    *
    * @param alleles one or more allele strings
    */
  def apply(alleles: String*): AlleleSet = {
    require(alleles.nonEmpty, "Must provide at least one allele.")
    AlleleSet(Allele(alleles.head), alleles.drop(1).map(s => Allele(s)))
  }
}
