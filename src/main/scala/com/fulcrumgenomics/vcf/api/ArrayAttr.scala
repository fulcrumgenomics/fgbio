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

import scala.collection.compat._
import scala.reflect.ClassTag

object ArrayAttr {
  /** Constructs an instance with the values supplied. */
  def apply[A : ClassTag](values: IterableOnce[A]): ArrayAttr[A] = {
    new ArrayAttr[A](values.iterator.toArray)
  }

  /** Apply method that re-uses the supplied array.  Should only be used from within the `api` package
    * and only when it can be guaranteed no other references to the array exist and no other caller
    * can modify the values.
    */
  private[api] def apply[A](values: Array[A]): ArrayAttr[A] = new ArrayAttr[A](values)
}


/**
  * Class that is used to store multi-valued attributes from a VCF (e.g. `AD=10,20`) and correctly
  * handle missing values.
  *
  * It is possible for one or all values in the collection to be missing. If accessed directly, e.g. by index
  * or by iteration, missing values will return one of the following based on the type of the attribute:
  *   - [[Variant.Missing]] for Strings
  *   - [[Variant.MissingInt]] for Ints
  *   - [[Variant.MissingFloat]] for Floating point numbers
  *
  * If you need to deal with the possibility of missing values it is strongly recommended that you use
  * the [[isMissing()]] and/or [[isDefined()]] methods or use the [[get()]] which returns an option type.
  *
  * @param values the values stored in the collection
  * @tparam A the type of values stored in the collection
  */
class ArrayAttr[A] private(private val values: Array[A]) extends IndexedSeq[A] {
  /** True if there are any missing values in the collection. */
  def hasMissingValues: Boolean = values.exists(Variant.isMissingValue)

  /** True if the element at the index is missing, false otherwise. */
  def isMissing(idx: Int): Boolean = Variant.isMissingValue(values(idx))

  /** True if the element at the index is not missing, false otherwise. */
  def isDefined(idx: Int): Boolean = !isMissing(idx)

  /** The number of elements (including missing) in the collection. */
  override def length: Int = values.length

  /** Accesses the value at the specified index. May return a Missing value if no value is defined at the index.
    * To avoid dealing with Missing values use [[isMissing(idx)]] or [[isDefined(idx)]] prior to accessing the
    * element, or use [[get()]] instead.
    */
  override def apply(idx: Int): A = this.values(idx)

  def get(idx: Int): Option[A] = if (isDefined(idx)) Some(apply(idx)) else None
}

