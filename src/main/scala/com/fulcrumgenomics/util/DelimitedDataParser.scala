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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import dagr.commons.reflect.ReflectionUtil
import dagr.commons.util.LazyLogging

import scala.reflect.runtime.{universe => ru}

/**
  * Represents a row of parsed data.  Provides methods for accessing values in a type-safe
  * way either via apply() methods for non-optional fields or via get for optional fields.
  */
class Row private[util] (private val headerIndices: Map[String,Int], private val fields: Array[String], val trim: Boolean) {
  /* Internal method to pull a value out of the field array by index and trim it if requested. */
  private def value(index: Int) = {
    if (index > fields.length-1) throw new IndexOutOfBoundsException(s"Invalid column index supplied: ${index}")
    val string = fields(index)
    if (trim) string.trim else string
  }

  /** Fetches a value of the desired type by the 0-based column index. */
  def apply[A](index: Int)(implicit tag: ru.TypeTag[A]): A = {
    val c = ReflectionUtil.typeTagToClass(tag)
    ReflectionUtil.constructFromString(c, c, value(index)).get.asInstanceOf[A]
  }

  /** Fetches a value of the desired type by column name. */
  def apply[A](column: String)(implicit tag: ru.TypeTag[A]): A = apply(headerIndices(column))

  /** Gets a value from the column with the specified index. If the value is null or empty returns None. */
  def get[A](index: Int)(implicit tag: ru.TypeTag[A]): Option[A] = {
    val string = value(index)
    if (string == null || string.isEmpty) {
      None
    }
    else {
      val c = ReflectionUtil.typeTagToClass(tag)
      Some(ReflectionUtil.constructFromString(c, c, string).get.asInstanceOf[A])
    }
  }

  /** Gets a value from the specified column. If the value is null or empty returns None. */
  def get[A](column: String)(implicit tag: ru.TypeTag[A]): Option[A] = get(headerIndices(column))
}


object DelimitedDataParser {
  val DefaultBlankPolicy: Boolean = true
  val DefaultTrim :       Boolean = true

  /** Constructs a DelimitedDataParser for a path. */
  def apply(path: FilePath, delimiter: Char): DelimitedDataParser =
    new DelimitedDataParser(Io.toSource(path).getLines(), delimiter=delimiter)

  /** Constructs a DelimitedDataParser for a path. */
  def apply(lines: Traversable[String], delimiter: Char): DelimitedDataParser =
    new DelimitedDataParser(lines, delimiter=delimiter)
}

/**
  * A parser for files of text columns delimited by some character (e.g. tab-delimited or csv).
  *
  * @param lines the lines from the file
  * @param delimiter the delimiter between columns
  * @param ignoreBlankLines whether blank lines should be ignored
  * @param trimFields whether individual fields should have their String values trimmed
  */
class DelimitedDataParser(lines: TraversableOnce[String],
                          val delimiter: Char,
                          val ignoreBlankLines: Boolean = DelimitedDataParser.DefaultBlankPolicy,
                          val trimFields: Boolean = DelimitedDataParser.DefaultTrim) extends Iterator[Row] with LazyLogging {

  private val _lines = if (ignoreBlankLines) lines.toIterator.filter(_.nonEmpty) else lines.toIterator

  // An array of the headers
  private val _headers = {
    if (this._lines.hasNext) this._lines.next().split(delimiter).map(h => if (trimFields) h.trim else h).toIndexedSeq
    else IndexedSeq.empty
  }

  // A temporary array used for parsing each line into fields
  private val tmp = new Array[String](headers.size * 2)

  // A map of header name to column index
  private val headerIndices: Map[String,Int] = _headers.zipWithIndex.toMap

  /** Returns the headers encountered in the file, in order. */
  def headers: Seq[String] = this._headers

  /** True if there is another line available for parsing. */
  override def hasNext: Boolean = this._lines.hasNext

  /** Retrieves the next row. */
  override def next(): Row = {
    val line   = this._lines.next()
    val fields = split(line)
    if (fields.length != _headers.length) {
      logger.error(s"Line has incorrect number of fields. Header length is ${_headers.length}, row length is ${fields.length}")
      logger.error(s"Headers: ${_headers.mkString(delimiter.toString)}")
      logger.error(s"Line: ${line}")
      throw new IllegalStateException(s"Incorrect number of values in line.")
    }

    new Row(headerIndices, fields, trim=trimFields)
  }

  /** Splits the line into it's constituent fields. */
  protected[util] def split(line: String): Array[String] = {
    val cs = line.toCharArray
    val len = cs.length
    var i = 0
    var count = 0
    var start = 0
    var end = 0

    while (i <= len && count < this.tmp.length) {
      if (i == len || cs(i) == this.delimiter) {
        this.tmp(count) = new String(cs, start, end-start)
        count += 1
        start = i + 1
        end   = i + 1
      }
      else {
        end += 1
      }

      i += 1
    }

    val fields = new Array[String](count)
    System.arraycopy(this.tmp, 0, fields, 0, count)
    fields
  }
}
