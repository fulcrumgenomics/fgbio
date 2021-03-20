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

object Variant {
  /** Value used in VCF for values that are missing. */
  val Missing: String = "."

  /** Value used in VCF for values that are missing. */
  val MissingChar: Char = '.'

  /** Value used in arrays of ints for missing values. */
  val MissingInt: Int = Int.MinValue

  /** Value used in arrays of floats for missing values. This is the same value used in the BCF specification
    * and is the "quiet NaN"; i.e. it is not a number, but the exact value is different than Float.NaN.
    *
    * To check for this value you must use `value.equals(MissingFloat)` NOT `value == MissingFloat` as the latter
    * will always return false.
    */
  val MissingFloat: Float = {
    import java.lang.{Float => JFloat, Long => JLong}
    val l = JLong.parseLong("7F800001", 16)
    JFloat.intBitsToFloat(l.intValue())
  }

  /** Used in comparisons in [[isMissingValue()]] since we can't compare NaNs using == or .equals(). */
  private val MissingFloatIntBits: Int = java.lang.Float.floatToRawIntBits(MissingFloat)

  /** Returns true if the passed value is one of the missing values. */
  def isMissingValue(x: Any): Boolean = x match {
    case ch: Char   => ch == '.'
    case s: String  => s == "."
    case i: Int     => i == MissingInt
    case f: Float   => java.lang.Float.floatToRawIntBits(f) == MissingFloatIntBits
    case _          => false
  }

  /** The value stored for fields of type `Flag` when stored in a Map. */
  val FlagValue: String = Missing

  /** The set of filters that are applied to passing variants. */
  val PassingFilters: Set[String] = Set("PASS")

  // Empty collections used in default constructors
  private[api] val EmptyFilters:    Set[String]           = Set.empty
  private[api] val EmptyInfo:       ListMap[String, Any]  = ListMap.empty
  private[api] val EmptyGtAttrs:    Map[String, Any]      = Map.empty
  private[api] val EmptyGenotypes:  Map[String, Genotype] = Map.empty
}


/**
  * Represents a variant from a VCF or similar source.
  *
  * Attributes from the INFO field are available via the [[apply()]], [[get()]], and [[getOrElse()]] methods.
  * Attribute that are defined but missing in the VCF (e.g. DP=.) are not stored.  All multi-valued attributes
  * are stored as instances of [[ArrayAttr]] which is a custom subclass of [[IndexedSeq]] which provides methods
  * to help in dealing with missing values in cases where a subset of values are missing (e.g. AD=10,.,20).
  *
  * @param chrom the chromosome on which a variant resides.
  * @param pos the start position of the variant
  * @param id the optional ID (e.g. rs number) of the variant
  * @param alleles the set of alleles recorded for the variant
  * @param qual the optional phred-scaled quality score of the variant
  * @param filters the set of filters applied to the variant.  An empty set indicates the variant has not had
  *               filtration applied.  A single value of "PASS" is applied if the variant passes filters, and
  *               one or more non-PASS strings if the variant fails filters.
  * @param attrs a map of attributes that come from the INFO field of a VCF
  * @param genotypes a map of sample name -> genotype for the variant
  */
final case class Variant(chrom: String,
                         pos: Int,
                         id: Option[String] = None,
                         alleles: AlleleSet,
                         qual: Option[Double] = None,
                         filters: Set[String] = Variant.EmptyFilters,
                         attrs: ListMap[String,Any] = Variant.EmptyInfo,
                         genotypes: Map[String, Genotype] = Variant.EmptyGenotypes
                        ) {

  /** The end position of the variant based on either the `END` INFO field _or_ the length of the reference allele. */
  val end: Int = get[Int]("END").getOrElse(pos + alleles.ref.length - 1)

  /** Retrieves a value from the INFO map.  Will throw an exception if the key does not exist. */
  def apply[A](key: String): A = attrs(key).asInstanceOf[A]

  /** Retrieves an optional value from the INFO map.  Will return [[None]] if the key does not exist. */
  def get[A](key: String): Option[A] = attrs.get(key).asInstanceOf[Option[A]]

  /** Retrieves an optional value from the INFO map.  Will return `default` if the key does not exist. */
  def getOrElse[A](key: String, default: => A): A = attrs.getOrElse(key, default).asInstanceOf[A]

  /** Convenience method to get a single sample's genotype. */
  def gt(sample: String): Genotype = this.genotypes(sample)

  /** Returns an iterator over the genotypes for this variant. */
  def gts: Iterator[Genotype] = this.genotypes.valuesIterator
}
