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
package com.fulcrumgenomics.fastq

/**
  * Represents a record that can be read from or written to a fastq file.
  */
case class FastqRecord(name: String, bases: String, quals: String, comment: Option[String] = None, readNumber: Option[Int] = None) {
  /** Constructs the header line from the name, read number and comment. */
  def header: String = name + readNumber.map("/" + _).getOrElse("") + comment.map(" " + _).getOrElse("")

  /** Returns the length of the sequence represented by this record. */
  def length: Int = bases.length

  /** Trims the record to a given length. If the record is already shorter than the length provided,
    * returns a record at the current length.
    * @param len the length to trim to
    * @return a record with length <= len
    */
  def trimmedTo(len: Int) : FastqRecord = {
    if (len >= this.length) this
    else copy(bases=this.bases.substring(0, len), quals=this.quals.substring(0, len))
  }

  override def toString: String = s"@$header\n$bases\n+\n$quals"
}
