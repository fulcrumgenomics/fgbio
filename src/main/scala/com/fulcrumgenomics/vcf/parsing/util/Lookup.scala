/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.vcf.parsing.util

import com.fulcrumgenomics.vcf.parsing.util.ParseResult.{fail, parseOrExcept, success}
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex.{Index, SI}
import fastparse.P

/** Provides lookups in a key-value map of strings, while retaining the start and end index of the values from
  * parsing.  This allows for validation to occur following parsing of key-value pairs in a VCF header entry/line.
  *
  * @param lookup map from [[String]] key to a [[String]] value wrapped with where it was parsed in the input
  * @param startIndex the start index in the input to report upon any lookup failure.
  */
case class Lookup private (private val lookup: Map[String, SI], startIndex: Index)(implicit ctx: P[_]) {
  /** Iterator over keys and values. */
  def iterator: Iterator[(String, SI)] = this.lookup.iterator

  /** Optionally returns the [[String]] associated with a key.
    *
    * @param key the key
    * @return a parser with an option value containing the value associated with key in this map, or None if none exists.
    */
  def get(key: String): P[Option[String]] = lookup.get(key) match {
    case None         => success(None)
    case Some(result) => success(value=Some(result.value))
  }

  /** Retrieves the value which is associated with the given key.
    *
    *  @param  key the key
    *  @return     a parser with a value associated with the given key, or failure if none exists.
    */
  def apply(key: String): P[String] = lookup.get(key) match {
    case Some(result) => success(value=result.value)
    case None         => fail(startIndex=startIndex, message="key `$key` not found")
  }

  /** Retrieves the [[String]] and where it was parsed in the input associated with the given key.
    *
    *  @param  key the key
    *  @return     the value and where it was parsed in the input which is associated with the given key
    */
  @throws[NoSuchElementException]
  def si(key: String): SI = this.lookup(key)

  /** Optionally returns the [[String]] and where it was parsed in the input associated with the given key.
    *
    *  @param  key the key
    *  @return     an option value and where it was parsed in the input associated with key in this map, or None if
    *              none exists.
    */
  def getSI(key: String): Option[SI] = this.lookup.get(key)

  /** Optionally returns the [[String]] associated with a key.  Fails if the
    * value associated with the key exists but the parser fails.
    *
    * @param key the key
    * @param f a method to transform a [[String]] value into the type
    * @return a parser with an option value containing the transformed value associated with key in this map, or None
    *         if none exists.
    */
  def getAndMap[T](key: String, f: String => T): P[Option[T]] = {
    lookup.get(key) match {
      case None         => success(None)
      case Some(result) => applyF(key=key, result=result, f=f).map(Some(_))
    }
  }

  /** Returns a parser with an optional value transformed from the [[String]] associated with a key.  Fails if the
    * value associated with the key exists but the parser fails.
    *
    * @param key the key
    * @param parser a parser to parse a [[String]] value into the type
    */
  def getAndParse[T](key: String, parser: P[_] => P[T]): P[Option[T]] = {
    getAndMap(
      key = key,
      f   = input => parseOrExcept[T](input=input, parser=parser)
    )
  }

  /** Returns the [[String]] associated with a key.  Fails if the
    * value associated with the key exists but the parser fails.
    *
    * @param key the key
    * @param f a method to transform a [[String]] value into the type
    * @return a parser containing the transformed value associated with key in this map, or None
    *         if none exists.
    */
  def andMap[T](key: String, f: String => T): P[T] = {
    lookup.get(key) match {
      case None         => fail(startIndex=startIndex, message="key `$key` not found")
      case Some(result) => applyF(key=key, result=result, f=f)
    }
  }

  /** Returns a parser with the value transformed from the [[String]] associated with a key.  Fails if the
    * value associated with the key exists but the parser fails.
    *
    * @param key the key
    * @param parser a parser to parse a [[String]] value into the type
    */
  def andParse[T](key: String, parser: P[_] => P[T]): P[T] = {
    andMap(
      key = key,
      f   = input => parseOrExcept[T](input=input, parser=parser)
    )
  }

  /** Return the transformed value associated with the given key.  Fail if the value could not be transformed from the
    * [[String]] value with the given method, failing at the given start index.
    * @param key the key associated with the value
    * @param result the result associated with the key
    * @param f the method to transform the [[String]] value
    * @tparam T the type of the return value
    */
  private def applyF[T](key: String, result: SI, f: String => T): P[T] = {
    try {
      success(f(result.value))
    } catch {
      case _: NumberFormatException =>
        result.fail(message=f"key `$key`: value is not a number: `${result.value}`")
      case ex: Exception            =>
        result.fail(message=f"key `$key`: value could not be parsed `${result.value}`: ${ex.getMessage}")
    }
  }
}

object Lookup {
  /** Build a [[Lookup]] based on the given the key-value pairs.  Fails if the keys are not unique, succeeds otherwise.
    *
    * @param attrs the key-value pairs
    */
  def parse(attrs: Seq[(SI, SI)])(implicit ctx: P[_]): P[Lookup] = {
    // require distinct keys
    attrs.groupBy(_._1).find(_._2.length > 1) match {
      case None =>
        val map: Map[String, SI] = attrs.map { case (key, value) => (key.value, value) }.toMap
        val startIndex: Index = attrs.head._1.startIndex
        ctx.freshSuccess(Lookup(map, startIndex=startIndex))
      case Some((key, values: Seq[(SI, SI)])) =>
        values(1)._2.fail(f"found `${values.length}` values for key `$key`: " + values.map(_._2.value).mkString(", "))
    }
  }
}
