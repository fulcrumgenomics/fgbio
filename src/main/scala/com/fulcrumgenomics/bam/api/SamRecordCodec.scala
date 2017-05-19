/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.api

import java.io.{ByteArrayInputStream, ByteArrayOutputStream}

import com.fulcrumgenomics.util.Sorter
import htsjdk.samtools.{BAMRecordCodec, SAMFileHeader}

/** Sorter.Codec implementation that wraps HTSJDK's BAMRecordCodec to read/write records to bytes. */
class SamRecordCodec(header: SAMFileHeader, maxRecordSize: Int = 128 * 1024) extends Sorter.Codec[SamRecord] {
  private val out      = new ByteArrayOutputStream(maxRecordSize)
  private val bamCodec = new BAMRecordCodec(header, SamRecord.Factory)

  /** Encode the object into an array of bytes. */
  override def encode(a: SamRecord): Array[Byte] = this.synchronized {
    out.reset()
    bamCodec.setOutputStream(out)
    bamCodec.encode(a.asSam)
    out.toByteArray
  }

  /** Decode an object from an array of bytes. */
  override def decode(bs: Array[Byte], start: Int, length: Int): SamRecord = this.synchronized {
    bamCodec.setInputStream(new ByteArrayInputStream(bs, start, length))
    bamCodec.decode().asInstanceOf[SamRecord]
  }
}
