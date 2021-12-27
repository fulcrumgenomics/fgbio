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

import java.io.{BufferedInputStream, BufferedReader, FileInputStream, InputStreamReader}
import java.util.zip.GZIPInputStream

import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.BlockCompressedInputStream

/**
  * Tests for the Io utility class.
  */
class IoTest extends UnitSpec {

  Seq(".bgz", ".bgzip").foreach { ext =>
    it should s"round trip data to a bgzipped file with extension ${ext}" in {
      val text = "This is a stupid little text fragment for compression. Yay compression!"
      val data = Seq.fill(10)(text)
      val f = makeTempFile("test.", ext)
      Io.writeLines(f, data)

      val stream = new BufferedInputStream(new FileInputStream(f.toFile))
      BlockCompressedInputStream.isValidFile(stream) shouldBe true
      val reread = Io.readLines(f).toIndexedSeq

      reread shouldBe data
    }
  }
}
