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

import com.fulcrumgenomics.vcf.api.Allele.{NoCallAllele, SpannedAllele}

/** A genotype for a variant.
  *
  * @param alleles the set of alleles for the variant
  * @param sample the name of the sample for which this is a genotype
  * @param calls the called alleles for the sample
  * @param phased whether or not the calls are phased
  */
final case class Genotype(alleles: AlleleSet,
                          sample: String,
                          calls: IndexedSeq[Allele],
                          phased: Boolean = false,
                          attrs: Map[String, Any] = Variant.EmptyGtAttrs
                         ) {
  require(calls.nonEmpty, "Genotype must have ploidy of at least 1!.")

  /** The indices of the calls within the AlleleSet. If any allele is no-called, that is returned as -1. */
  private [api] val callIndices: IndexedSeq[Int] = calls.map(c => if (c == NoCallAllele) -1 else alleles.indexOf(c))

  /** Returns the separator that should be used between alleles when displaying the genotype. */
  private def separator: String = if (phased) "|" else "/"

  /** Returns the set of alleles for the genotype filtering out no-calls and spanned alleles. */
  private def calledAlleles: Seq[Allele] = calls.filter(a => a != NoCallAllele && a != SpannedAllele)

  /** The ploidy of the called genotype - equal to calls.length. */
  def ploidy: Int = calls.length

  /** True if all [[calls]] are no-call alleles. */
  def isNoCall: Boolean = calls.forall(_ == Allele.NoCallAllele)

  /** True if none of [[calls]] are no-call alleles. */
  def isFullyCalled: Boolean = !calls.contains(Allele.NoCallAllele)

  /** True if all [[calls]] are the reference allele. */
  def isHomRef: Boolean = calls.forall(_ == alleles.ref)

  /** True if all [[calls]] are the same non-reference allele. */
  def isHomVar: Boolean = isFullyCalled && calls(0) != alleles.ref && calls.forall(_ == calls(0))

  /** True if the genotype contains at least two called alleles that are different from one another. */
  def isHet: Boolean = {
    val actuallyCalled = calledAlleles
    actuallyCalled.exists(_ != actuallyCalled.head)
  }

  /** True if the genotype is fully called, the genotype is heterozygous and none of the called alleles are the reference. */
  def isHetNonRef: Boolean = isHet && !calls.contains(alleles.ref)

  /** Retrieves a value from the INFO map.  Will throw an exception if the key does not exist. */
  def apply[A](key: String): A = attrs(key).asInstanceOf[A]

  /** Retrieves an optional value from the INFO map.  Will return [[None]] if the key does not exist. */
  def get[A](key: String): Option[A] = attrs.get(key).asInstanceOf[Option[A]]

  /** Retrieves an optional attribute.  Will return `default` if the key does not exist. */
  def getOrElse[A](key: String, default: => A): A = attrs.getOrElse(key, default).asInstanceOf[A]

  /** Yields a genotype string using bases, e.g. A/C, or CTTT|C. */
  def gtWithBases: String = calls.map(_.toString).mkString(separator)

  /** Yields a genotype string using allele indices, e.g. 0/1. */
  def gtVcfStyle: String = callIndices.iterator.map(i => if (i >= 0) i.toString else ".").mkString(separator)
}
