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

import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api.Variant.ArrayAttribute

import scala.collection.immutable.ListMap

object Variant {
  /** Value used in VCF for values that are missing. */
  val Missing: String = "."

  /** The value stored for fields of type `Flag` when stored in a Map. */
  val FlagValue: String = Missing

  /** The type used in attribute maps (INFO, genotype) to store multi-valued attributed. */
  type ArrayAttribute[A] = IndexedSeq[A]

  /** Function to convert any scala collection to an [[com.fulcrumgenomics.vcf.api.Variant.ArrayAttribute]]. */
  def toArrayAttribute[A](values: TraversableOnce[A]): ArrayAttribute[A] = values.toIndexedSeq

  /** The set of filters that are applied to passing variants. */
  val PassingFilters: Set[String] = Set("PASS")
}


/**
  * Represents a variant from a VCF or similar source.
  *
  * @param chrom the chromosome on which a variant resides.
  * @param pos the start position of the variant
  * @param id the optional ID (e.g. rs number) of the variant
  * @param alleles the set of alleles recorded for the variant
  * @param qual the optional phred-scaled quality score of the variant
  * @param filter the set of filters applied to the variant.  An empty set indicates the variant has not had
  *               filtration applied.  A single value of "PASS" is applied if the variant passes filters, and
  *               one or more non-PASS strings if the variant fails filters.
  * @param attributes a map of attributes that come from the INFO field of a VCF
  * @param genotypes a map of sample name -> genotype for the variant
  */
case class Variant(chrom: String,
                   pos: Int,
                   id: Option[String],
                   alleles: AlleleSet,
                   qual: Option[Double],
                   filter: Set[String],
                   attributes: ListMap[String,Any],
                   genotypes: Map[String, Genotype]
                  ) {

  /** The end position of the variant based on either the `END` INFO field _or_ the length of the reference allele. */
  val end: Int = get[Int]("END").getOrElse(pos + alleles.ref.length - 1)

  /** Retrieves a value from the INFO map.  Will throw an exception if the key does not exist. */
  def apply[A](key: String): A = attributes(key).asInstanceOf[A]

  /** Retrieves an optional value from the INFO map.  Will return [[None]] if the key does not exist. */
  def get[A](key: String): Option[A] = attributes.get(key).asInstanceOf[Option[A]]

  /** Retrieves an optional value from the INFO map.  Will return `default` if the key does not exist. */
  def getOrElse[A](key: String, default: => A): Option[A] = attributes.getOrElse(key, default).asInstanceOf[Option[A]]
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

  /** The indices of the calls within the AlleleSet. If any allele is no-called, that is returned as -1. */
  val callIndices: IndexedSeq[Int] = calls.map(c => if (c == NoCallAllele) -1 else alleles.indexOf(c))

  /** The ploidy of the called genotype - equal to calls.length. */
  def ploidy: Int = calls.length

  /** True if all [[calls]] are no-call alleles. */
  def isNoCall: Boolean = calls.forall(_ == Allele.NoCallAllele)

  /** True if none of [[calls]] are no-call alleles. */
  def isFullyCalled: Boolean = !calls.contains(Allele.NoCallAllele)

  /** True if all [[calls]] are the reference allele. */
  def isHomRef: Boolean = calls.forall(_ == alleles.ref)

  /** True if all [[calls]] are the same non-reference allele. */
  def inHomVar: Boolean = isFullyCalled && calls(0) != alleles.ref && calls.forall(_ == calls(0))

  /** True if the genotype is fully called and there are at least two distinct alleles called in the genotype. */
  def isHet: Boolean = isFullyCalled && calls.exists(_ != calls(0))

  /** True if the genotype is fully called, the genotype is heterozygous and none of the called alleles are the reference. */
  def isHetNonRef: Boolean = isFullyCalled && isHet && !calls.contains(alleles.ref)

  /** Retrieves a value from the INFO map.  Will throw an exception if the key does not exist. */
  def apply[A](key: String): A = attributes(key).asInstanceOf[A]

  /** Retrieves an optional value from the INFO map.  Will return [[None]] if the key does not exist. */
  def get[A](key: String): Option[A] = attributes.get(key).asInstanceOf[Option[A]]

  /** Retrieves an optional value from the INFO map.  Will return `default` if the key does not exist. */
  def getOrElse[A](key: String, default: => A): Option[A] = attributes.getOrElse(key, default).asInstanceOf[Option[A]]
}



