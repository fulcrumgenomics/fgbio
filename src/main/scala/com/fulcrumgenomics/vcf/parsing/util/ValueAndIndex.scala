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

import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex.Index
import fastparse.P
import com.fulcrumgenomics.vcf.parsing.util.ParseResult
/**
  * Stores a value of a given type along with the start and end index of the source string in the original input.  This
  * is useful when downstream validation of the value fails, and the start and end index want to be known.
  *
  * @param value the parsed value
  * @param start the start index in the original input (0-based)
  * @param end the end index in the original input (0-based)
  * @tparam T the type of the value
  */
case class ValueAndIndex[T](value: T, start: Index, end: Index) {
  /** The start index in the original input (0-based)  */
  def startIndex: Index = start
  /** The end index in the original input (0-based)  */
  def endIndex: Index = start

  /** Returns a parse failure with the given message and start index in the input. */
  def fail[_: P](message: String): P[Nothing] = ParseResult.fail(startIndex=this.startIndex, message=message)
}

object ValueAndIndex {
  /** Alias for an index into a [[String]] */
  type Index = Int
  /** Succinct alias for [[ValueAndIndex]] with a given type */
  type VI[T] = ValueAndIndex[T]
  /** Succinct alias for a [[ValueAndIndex]] for a [[String]] value */
  type SI    = ValueAndIndex[String]

  /** Converts the result of capturing the start index, value, and end index after parsing. */
  def toResultAndIndex[T](values: (Int, T, Int)): VI[T] = {
    ValueAndIndex(value=values._2, start=values._1, end=values._3)
  }
  /** Implicitly extracts the value from a [[ValueAndIndex]]. */
  implicit def toResult[T](valueAndIndex: ValueAndIndex[T]): T = valueAndIndex.value
  /** Implicitly extracts the [[String]] from a [[ValueAndIndex]] with [[String]] value. */
  implicit def stringAndIndexToResult(resultAndIndex: SI): String = resultAndIndex.value

  /** Adds a method to the start index, value, and end index after parsing.
    *
    * For example: `Index ~ CharPred(_ == "A-Z").rep.! ~ Index).si`
    * */
  implicit class CaptureResultAndIndex[T](val parse0: P[(Int, T, Int)]) {
    def si: P[VI[T]] = parse0.map(ValueAndIndex.toResultAndIndex[T])
  }
}