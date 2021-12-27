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

import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.Interval

class IntervalListWriterTest extends UnitSpec {

  private val dict = SequenceDictionary(
    SequenceMetadata(name="chr1", length=10000),
    SequenceMetadata(name="chr2", length=50000)
  )

  "IntervalListWriter" should "write a set of intervals" in {
    val path   = makeTempFile("intervals.", ".interval_list")
    val writer = IntervalListWriter(path, dict)

    // NB: order does not matter here
    val intervals = Seq(
      new Interval("chr2", 1, 10),
      new Interval("chr1", 1, 10),
      new Interval("chr1", 1000, 2000),
      new Interval("chr1", 50, 5000)
    )
      
    writer ++= intervals
    writer.close()

    IntervalListSource(path).toSeq should contain theSameElementsInOrderAs intervals
  }

  it should "fail if no sequence dictionary is given" in {
    val path      = makeTempFile("intervals.", ".interval_list")
    val header    = new SAMFileHeader()
    val exception = intercept[Exception] { IntervalListWriter(path, header) }

    exception.getMessage should include("The header does not contain any reference sequences")
  }

  it should "fail if an interval's contig is not found in the sequence dictionary" in {
    val path   = makeTempFile("intervals.", ".interval_list")
    val writer = IntervalListWriter(path, dict)

    val exception = intercept[Exception] { writer += new Interval("chr3", 1, 10) }
    exception.getMessage should include("Reference sequence not found")

    writer.close()
  }
}
