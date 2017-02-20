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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}

class RemoveSamTagsTest extends UnitSpec {

  "RemoveSamTags" should "run end-to-end" in {
    val builder = new SamRecordSetBuilder()
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("OK", "foo"); rec.setAttribute("NO", "bar") }
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("OK", "foo"); }
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("NO", "bar") }
    builder.addFrag(unmapped=true)

    val input = builder.toTempFile()
    val output = makeTempFile("RemoveSamTags", ".bam")
    val tagsToRemove = Seq("NO")

    new RemoveSamTags(input=input, output=output, tagsToRemove=tagsToRemove).execute()

    val records = readBamRecs(output)
    records should have size 4
    records.count(_.getAttribute("OK") != null) shouldBe 2
    records.count(_.getAttribute("NO") == null) shouldBe 4
  }

  it should "run end-to-end with no tags to remove" in {
    val builder = new SamRecordSetBuilder()
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("OK", "foo"); rec.setAttribute("NO", "bar") }
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("OK", "foo"); }
    builder.addFrag(unmapped=true).foreach { rec => rec.setAttribute("NO", "bar") }
    builder.addFrag(unmapped=true)

    val input = builder.toTempFile()
    val output = makeTempFile("RemoveSamTags", ".bam")

    new RemoveSamTags(input=input, output=output).execute()

    val records = readBamRecs(output)
    records should have size 4
    records.count(_.getAttribute("OK") != null) shouldBe 2
    records.count(_.getAttribute("NO") != null) shouldBe 2
  }
}
