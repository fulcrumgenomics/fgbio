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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.IntervalList

class IntervalListSourceTest extends UnitSpec {
  import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary

  private val dict = SequenceDictionary(
    SequenceMetadata(name="chr1", length=10000),
    SequenceMetadata(name="chr2", length=50000)
  )

  private val intervals = new IntervalList(this.dict.asSam)

  "IntervalListSource" should "read an interval list from a path" in {
    val path = makeTempFile("intervals.", ".interval_list")
    this.intervals.write(path)

    val source = IntervalListSource(path)
    val actual = source.toIndexedSeq
    source.close()

    this.intervals.getHeader.getSequenceDictionary.assertSameDictionary(source.dict.asSam)
    this.dict.asSam.assertSameDictionary(source.dict.asSam)
    actual should contain theSameElementsInOrderAs this.intervals.getIntervals.toIndexedSeq

    val list = IntervalListSource(path).toIntervalList
    this.intervals.getHeader.getSequenceDictionary.assertSameDictionary(list.getHeader.getSequenceDictionary)
    this.dict.asSam.assertSameDictionary(list.getHeader.getSequenceDictionary)
    list.getIntervals.toIndexedSeq should contain theSameElementsInOrderAs this.intervals.getIntervals.toIndexedSeq
  }

  it should "fail if no header is present" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("chr1\t1\t1\t+\tname")) }
    exception.getMessage should include("No header found")
  }

  it should "fail if no sequence dictionary is present" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "chr1\t1\t1\t+\tname")) }
    exception.getMessage should include("No reference sequences found")
  }

  it should "fail if a line does not have exactly five fields" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr1")).toSeq }
    exception.getMessage should include("Expected 5 fields on line 3")
  }

  it should "fail if an interval's contig is not found the in the sequence dictionary" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr2\t1\t1\t+\tname")).toSeq }
    exception.getMessage should include("key not found: chr2")
  }

  it should "fail if an interval's start is less than one" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr1\t0\t1\t+\tname")).toSeq }
    exception.getMessage should include("Start is less than 1")
  }

  it should "fail if an interval's end is beyond the contig's length in the sequence dictionary" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr1\t1\t10001\t+\tname")).toSeq }
    exception.getMessage should include("End is beyond")
  }

  it should "fail if an interval's start is greater than its end" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr1\t2\t1\t+\tname")).toSeq }
    exception.getMessage should include("Start is greater than end")
  }

  it should "fail if the strand cannot be recognized" in {
    val exception = intercept[Exception] { IntervalListSource(Seq("@HD\tVN:1.0", "@SQ\tSN:chr1\tLN:10000", "chr1\t1\t1\tN\tname")).toSeq }
    exception.getMessage should include("Unrecognized strand")
  }
}
