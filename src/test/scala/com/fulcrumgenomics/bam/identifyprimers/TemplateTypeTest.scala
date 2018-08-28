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

package com.fulcrumgenomics.bam.identifyprimers

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class TemplateTypeTest extends UnitSpec with OptionValues {
  import TemplateType._

  "TemplateType" should "be build from SamRecords" in {
    val builder = new SamBuilder()

    // frags
    TemplateType(builder.addFrag(unmapped = true).value, None) shouldBe UnmappedFragment
    TemplateType(builder.addFrag(unmapped = false).value, None) shouldBe MappedFragment

    // pairs
    Seq(
      (true, true, UnmappedPair),
      (true, false, Unpaired),
      (false, true, Unpaired),
      (false, false, MappedPair)
    ).foreach { case (unmapped1, unmapped2, expectedTemplateType) =>
      val Seq(r1, r2) = builder.addPair(unmapped1=unmapped1, unmapped2=unmapped2)
      TemplateType(r1, Some(r2)) shouldBe expectedTemplateType
    }
  }

}
