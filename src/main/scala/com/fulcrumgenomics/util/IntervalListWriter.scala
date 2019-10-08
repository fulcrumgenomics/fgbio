/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

package com.fulcrumgenomics.util


import java.io.{BufferedWriter, Closeable}

import com.fulcrumgenomics.FgBioDef.PathToIntervals
import com.fulcrumgenomics.commons.io.Writer
import htsjdk.samtools.{SAMFileHeader, SAMSequenceDictionary, SAMTextHeaderCodec}
import htsjdk.samtools.util.Interval


object IntervalListWriter {
  /** Constructs an [[IntervalListWriter]] that will write to the provided path. */
  def apply(path: PathToIntervals, header: SAMFileHeader): IntervalListWriter = apply(Io.toWriter(path), header)

  /** Constructs an [[IntervalListWriter]] that will write to the provided path. */
  def apply(path: PathToIntervals, dict: SAMSequenceDictionary): IntervalListWriter = apply(Io.toWriter(path), dict)

  /** Constructs an [[IntervalListWriter]] from a Writer. */
  def apply(writer: java.io.Writer, header: SAMFileHeader): IntervalListWriter = writer match {
    case bw: BufferedWriter => new IntervalListWriter(bw, header)
    case w                  => new IntervalListWriter(new BufferedWriter(w), header)
  }

  /** Constructs an [[IntervalListWriter]] from a Writer. */
  def apply(writer: java.io.Writer, dict: SAMSequenceDictionary): IntervalListWriter = {
    val header = new SAMFileHeader()
    header.setSequenceDictionary(dict)
    this.apply(writer, header)
  }
}

/**
  * Implements a writer for interval lists.
  */
class IntervalListWriter(val out: BufferedWriter, header: SAMFileHeader) extends Writer[Interval] {

  // Write out the header
  new SAMTextHeaderCodec().encode(out, header)

  private val dict = header.getSequenceDictionary

  require(!dict.isEmpty, "The header does not contain any reference sequences.")

  override def write(interval: Interval): Unit = {
    require(dict.getSequence(interval.getContig) != null,
      s"Reference sequence not found in the sequence dictionary: ${interval.getContig}")

    out.write(interval.getContig)
    out.write('\t')
    out.write(Integer.toString(interval.getStart))
    out.write('\t')
    out.write(Integer.toString(interval.getEnd))
    out.write('\t')
    out.write(interval.getStrand.encode)
    out.write('\t')
    out.write(if (interval.getName != null) interval.getName else ".")
    out.newLine()
  }

  /** Closes the underlying writer. */
  override def close(): Unit = out.close()
}

