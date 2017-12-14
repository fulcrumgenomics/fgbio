/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics
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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class MakeUnmappedByInsertSizeTest extends UnitSpec {
  "MakeUnmappedByInsertSize" should "make fragment reads unmapped" in {
    val builder = new SamBuilder(readLength=10)
    builder.addFrag("q1", start=500)
    builder.addFrag("q2", unmapped=true)

    val out = makeTempFile("make_unmapped.", ".bam")
    new MakeUnmappedByInsertSize(input=builder.toTempFile(), output=out, minInsertSize=Some(0)).execute()
    val recs = readBamRecs(out)
    recs.filter(_.unmapped).map(_.name) should contain theSameElementsInOrderAs Seq("q1", "q2")
    recs.filter(_.mapped).map(_.name) shouldBe 'empty
  }

  it should "remove secondary and supplementary reads if they would become unmapped" in {
    val builder = new SamBuilder(readLength=10)
    builder.addPair("q1", start1=100, start2=590).foreach(_.secondary = true) // isize=500, secondary
    builder.addPair("q2", start1=100, start2=590).foreach(_.supplementary = true) // isize=500, supplementary
    builder.addPair("q3", start1=100, start2=140).foreach(_.secondary = true) // isize=50, secondary
    builder.addPair("q4", start1=100, start2=140).foreach(_.supplementary = true) // isize=50, supplementary

    val out = makeTempFile("make_unmapped.", ".bam")
    new MakeUnmappedByInsertSize(input=builder.toTempFile(), output=out, minInsertSize=Some(100)).execute()
    val recs = readBamRecs(out)
    recs should have size 4
    recs.filter(_.unmapped).map(_.name) shouldBe 'empty
    recs.filter(_.mapped).map(_.name) should contain theSameElementsInOrderAs Seq("q1", "q1", "q2", "q2")
  }

  it should "make unmapped reads below the specified insert size" in {
    val builder = new SamBuilder(readLength=10)
    builder.addPair("q1", start1=100, start2=589) // isize=499
    builder.addPair("q2", start1=100, start2=590) // isize=500
    builder.addPair("q3", start1=100, start2=591) // isize=501
    builder.addPair("q4", start1=100, start2=120) // isize=30
    builder.addPair("q5", start1=200, start2=100, strand1=Minus, strand2=Plus) // isize=110
    builder.addFrag("q6", start=500)
    builder.addPair("q7", start1=100, start2=500).foreach(r => r.refIndex = if (r.firstOfPair) 1 else 2) // isize=0/undefined

    val out = makeTempFile("make_unmapped.", ".bam")
    new MakeUnmappedByInsertSize(input=builder.toTempFile(), output=out, minInsertSize=Some(500)).execute()
    val recs = readBamRecs(out)
    recs should have size 13
    recs.filter(_.unmapped).map(_.name) should contain theSameElementsInOrderAs Seq("q1", "q1", "q4", "q4", "q5", "q5", "q6", "q7", "q7")
    recs.filter(_.mapped).map(_.name) should contain theSameElementsInOrderAs Seq("q2", "q2", "q3", "q3")
  }

  it should "make unmapped reads above the specified insert size" in {
    val builder = new SamBuilder(readLength=10)
    builder.addPair("q1", start1=100, start2=589) // isize=499
    builder.addPair("q2", start1=100, start2=590) // isize=500
    builder.addPair("q3", start1=100, start2=591) // isize=501
    builder.addPair("q4", start1=100, start2=120) // isize=30
    builder.addPair("q5", start1=200, start2=100, strand1=Minus, strand2=Plus) // isize=110
    builder.addFrag("q6", start=500)

    val out = makeTempFile("make_unmapped.", ".bam")
    new MakeUnmappedByInsertSize(input=builder.toTempFile(), output=out,  maxInsertSize=Some(499)).execute()
    val recs = readBamRecs(out)
    recs should have size 11
    recs.filter(_.unmapped).map(_.name) should contain theSameElementsInOrderAs Seq("q2", "q2", "q3", "q3", "q6")
    recs.filter(_.mapped).map(_.name) should contain theSameElementsInOrderAs Seq("q1", "q1", "q4", "q4", "q5", "q5")
  }
}
