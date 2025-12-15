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

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import htsjdk.samtools.{BAMRecordCodec, SAMFileHeader, SAMRecordIterator, SAMRecord}

import java.io.{ByteArrayInputStream, ByteArrayOutputStream}

/** An iterator over [[com.fulcrumgenomics.bam.api.SamRecord]]s that will automatically close the underlying iterator at
  * the end of iteration, and provides access to the [[htsjdk.samtools.SAMFileHeader]] from the associated source.
  */
final class SamIterator(val header: SAMFileHeader, underlying: SAMRecordIterator)
  extends SelfClosingIterator[SamRecord](SamIterator.buildIterator(underlying, header), () => underlying.close())
  with HeaderHelper

object SamIterator {
  // Lazy-initialized reusable resources for CRAM conversion (per iterator instance)
  private class CramConverter(header: SAMFileHeader) {
    private lazy val codec = new BAMRecordCodec(header, SamRecord.Factory)
    private lazy val buffer = new ByteArrayOutputStream(128 * 1024)

    def convert(plain: SAMRecord): SamRecord = {
      buffer.reset()
      codec.setOutputStream(buffer)
      codec.encode(plain)
      codec.setInputStream(new ByteArrayInputStream(buffer.toByteArray))
      codec.decode().asInstanceOf[SamRecord]
    }
  }

  private def buildIterator(underlying: SAMRecordIterator, header: SAMFileHeader): Iterator[SamRecord] = {
    // Create converter once per iterator - lazy vals inside will only initialize if CRAM records are encountered
    val converter = new CramConverter(header)

    underlying.map { rec =>
      rec match {
        case rec: SamRecord => rec  // Already enhanced (BAM/SAM) - no overhead
        case rec: SAMRecord => converter.convert(rec)  // CRAM - reuses codec and buffer
      }
    }
  }
}
