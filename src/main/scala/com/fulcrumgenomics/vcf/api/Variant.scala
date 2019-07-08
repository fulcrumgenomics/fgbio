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

import scala.collection.immutable.ListMap

case class Variant(chrom: String,
                   pos: Int,
                   id: Option[String],
                   alleles: AlleleSet,
                   qual: Option[Double],
                   filter: Seq[String],
                   info: ListMap[String,Any],
                   genotypes: Map[String, Genotype]
                  ) {

  val end: Int = get[Int]("END").getOrElse(pos + alleles.ref.length - 1)

  def apply[A](key: String): A = info(key).asInstanceOf[A]
  def get[A](key: String): Option[A] = info.get(key).asInstanceOf[Option[A]]
  def getOrElse[A](key: String, default: => A): Option[A] = info.getOrElse(key, default).asInstanceOf[Option[A]]
}


/** A genotype for a variant.
  *
  * @param alleles the set of alleles for the variant
  * @param sample the name of the sample for which this is a genotype
  * @param calls the called alleles for the sample
  * @param phased whether or not the calls are phased
  */
case class Genotype(alleles: AlleleSet,
                    sample: String,
                    calls: IndexedSeq[Allele],
                    phased: Boolean = false,
                    attributes: Map[String, Any]
                   ) {
  require(calls.nonEmpty, "Genotype must have ploidy of at least 1!.")

  val callIndices: IndexedSeq[Int] = calls.map(alleles.indexOf)

  def ploidy: Int = calls.length
  def isNoCall   : Boolean = calls.forall(_ == Allele.NoCallAllele)
  def isHomRef   : Boolean = calls.forall(_ == alleles.ref)
  def inHomVar   : Boolean = calls(0) != alleles.ref && calls.forall(_ == calls(0))
  def isHet      : Boolean = calls.exists(_ != calls(0))
  def isHetNonRef: Boolean = isHet && !calls.contains(alleles.ref)

  def apply[A](key: String): A = attributes(key).asInstanceOf[A]
  def get[A](key: String): Option[A] = attributes.get(key).asInstanceOf[Option[A]]
  def getOrElse[A](key: String, default: => A): Option[A] = attributes.getOrElse(key, default).asInstanceOf[Option[A]]

  // Would be nice to have typed accessors for common, spec-defined attributes, but how to deal with optionality?
  def dp: Option[Int] = get[Int]("DP")
  def DP: Int = apply[Int]("DP")
}



