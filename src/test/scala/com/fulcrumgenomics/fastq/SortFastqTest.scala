/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec

class SortFastqTest extends UnitSpec {
  def sort(fq: Seq[FastqRecord], maxInMem: Int = 1000): Seq[FastqRecord] = {
    val unsorted = makeTempFile("unsorted.", ".fastq")
    val sorted   = makeTempFile("sorted.", ".fastq")
    val out = FastqWriter(unsorted)
    out ++= fq
    out.close()

    new SortFastq(input=unsorted, output=sorted, maxRecordsInRam=maxInMem).execute()
    val in = FastqSource(sorted)
    yieldAndThen(in.toIndexedSeq)(in.close())
  }


  "SortFastq" should "not fail on an empty fastq file" in {
    sort(Seq()) shouldBe empty
  }

  it should "sort several records by name" in {
    val sorted = sort(Seq(
      FastqRecord("q5", "AAA", "###"),
      FastqRecord("q3", "AAA", "###"),
      FastqRecord("q6", "AAA", "###"),
      FastqRecord("q2", "AAA", "###"),
      FastqRecord("q1", "AAA", "###"),
      FastqRecord("q4", "AAA", "###")
    ))

    sorted.map(_.name) should contain theSameElementsInOrderAs Seq("q1", "q2", "q3", "q4", "q5", "q6")
  }

  it should "sort read 1s before read 2s with the same name" in {
    val sorted = sort(Seq(
      FastqRecord("q1", "AAA", "###", readNumber=Some(2)),
      FastqRecord("q2", "AAA", "###", readNumber=Some(2)),
      FastqRecord("q2", "AAA", "###", readNumber=Some(1)),
      FastqRecord("q1", "AAA", "###", readNumber=Some(1))
    ))

    sorted shouldBe Seq(
      FastqRecord("q1", "AAA", "###", readNumber=Some(1)),
      FastqRecord("q1", "AAA", "###", readNumber=Some(2)),
      FastqRecord("q2", "AAA", "###", readNumber=Some(1)),
      FastqRecord("q2", "AAA", "###", readNumber=Some(2))
    )
  }
}
