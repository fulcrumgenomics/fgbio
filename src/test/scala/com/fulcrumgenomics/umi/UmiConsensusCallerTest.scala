/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.GroupOrder
import org.scalatest.OptionValues

class UmiConsensusCallerTest extends UnitSpec with OptionValues {
  "UmiConsensusCaller" should "do nothing when the header has the right sort order in it" in {
    val builder = new SamBuilder(sort=Some(SamOrder.TemplateCoordinate))
    var warning : Option[String] = None
    var error   : Option[String] = None
    UmiConsensusCaller.checkSortOrder(builder.header, "foo.bam", w => warning = Some(w), e => error = Some(e))
    warning shouldBe None
    error shouldBe None
  }

  it should "fire a warning if the file is query grouped and unsorted but without sub-sort specified" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Unsorted))
    builder.header.setGroupOrder(GroupOrder.query)
    var warning : Option[String] = None
    var error   : Option[String] = None
    UmiConsensusCaller.checkSortOrder(builder.header, "foo.bam", w => warning = Some(w), e => error = Some(e))
    warning.value.indexOf("foo.bam") should be >= 0
    error shouldBe None
  }

  it should "fire an error if the file is coordinate sorted" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate))
    var warning : Option[String] = None
    var error   : Option[String] = None
    UmiConsensusCaller.checkSortOrder(builder.header, "foo.bam", w => warning = Some(w), e => error = Some(e))
    warning shouldBe None
    error.value.indexOf("foo.bam") should be >= 0
  }
}
