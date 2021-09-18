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

import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex._
import fastparse.NoWhitespace._
import fastparse.{&, CharIn, CharPred, Index, P, _}

/** Methods for parsing key-value pairs.  In particular, the values in `<...>` portion of a VCF header entry/line */
object Parsing {

  /** Parses one or more `key=value` pairs, comma-delimited. */
  def csvKeyValues[_: P]: P[Seq[(SI, SI)]] = (keyValue ~ ("," ~ csvKeyValues).rep(0)).map {
    case (key: SI, value: SI, attrs: Seq[Seq[(SI, SI)]]) => (key, value) +: attrs.flatten
  }

  /** Parses a full `key=value` pair.  See [[key]] and [[value]]. */
  def keyValue[_: P]: P[(SI, SI)] = key ~ "=" ~ value ~ &(CharIn(",>"))

  /** Parses a key in a `key=value` pair.  Parses one or more non-equal characters until the next character is an
    * equals, but does not consume the equals. */
  def key[_: P]: P[SI] = (Index ~ notEquals.rep(1).! ~ Index).si ~ &("=")

  /** Parses a value in a `key=value` pair.  The value may be quoted or unquoted (former preferred).  For quoted,
    * no quotations are allowed in the value.  For unquoted, parses until a comma, greater-than, or double-quote
    * is encountered. */
  def value[_: P]: P[SI] = quotedValue | unquotedValue

  /** Parses a quoted value (eg. `"description"`) */
  def quotedValue[_: P]: P[SI] = "\"" ~ (Index ~ CharPred(c => c != '"').rep(1).! ~ Index).si ~ "\""

  /** Parses an unquoted value (eg. `description`). */
  def unquotedValue[_: P]: P[SI] = notEndOfValue ~ &(endOfValue)

  /** Parses the next character after the end of a value in a key-value pair.  Must be a comma, greater-than, or double-quote
    * is encountered. For example, `INFO=<ID=value[,]Description="description"[,]Number=2>` or `fileformat="value["]`
    * where the end character parsed is in square brackets*/
  def endOfValue[_: P]: P[_] = CharIn(",>\"")

  /** Parses a full value in a key-value pair.  The end of the value occurs when a comma, greater-than, or double-quote
    * is encountered. */
  def notEndOfValue[_: P]: P[SI] = (Index ~ CharPred(!",>\"".contains(_)).rep.! ~ Index).si

  /** Parses any non-equal character. */
  def notEquals[_: P]: P[_] = CharPred(_ != '=')
}
