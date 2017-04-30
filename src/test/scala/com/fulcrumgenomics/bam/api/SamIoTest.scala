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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Io

class SamIoTest extends UnitSpec {
  "SamWriter/SamSource" should "write some records to a BAM on disk and read them back" in {
    val builder = new SamBuilder()
    builder.addPair(name="q1", start1=100, start2=300)
    builder.addPair(name="q4", start1=200, start2=400)
    builder.addPair(name="q3", start1=300, start2=500)
    builder.addPair(name="q2", start1=400, start2=600)

    val bam = makeTempFile("test.", ".bam")
    val out = SamWriter(bam, builder.header, sort=Some(SamOrder.Queryname))
    out ++= builder.iterator
    out.close()

    val in = SamSource(bam)
    in.toSeq.map(_.name) should contain theSameElementsInOrderAs  Seq("q1","q1","q2","q2","q3","q3","q4","q4")
  }

  it should "write some records to a SAM text file on disk and read them back" in {
    val builder = new SamBuilder()
    builder.addPair(name="q1", start1=100, start2=300)
    builder.addPair(name="q4", start1=200, start2=400)
    builder.addPair(name="q3", start1=300, start2=500)
    builder.addPair(name="q2", start1=400, start2=600)

    val sam = makeTempFile("test.", ".sam")
    val out = SamWriter(sam, builder.header, sort=Some(SamOrder.Coordinate))
    out ++= builder.iterator
    out.close()

    val in = SamSource(sam)
    in.toSeq.map(_.start) should contain theSameElementsInOrderAs Seq(100, 200, 300, 300, 400, 400, 500, 600)

    Io.readLines(sam).next() should startWith ("@HD")
  }
}
