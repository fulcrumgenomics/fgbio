/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}

/** Unit tests for [[VcfUtil]]. */
class VcfUtilTest extends UnitSpec {

  /** The name for test sample number 1. */
  private val sample1 = "sample1"

  /** The name for test sample number 2. */
  private val sample2 = "sample2"

  "VcfUtil.onlySample" should "return the only sample in a VCF source" in {
    val singleSample = VcfBuilder(samples = Seq(sample1)).toSource
    val doubleSample = VcfBuilder(samples = Seq(sample1, sample2)).toSource
    VcfUtil.onlySample(singleSample) shouldBe sample1
    singleSample.safelyClose()
    an[IllegalArgumentException] shouldBe thrownBy { VcfUtil.onlySample(doubleSample) }
  }

  "VcfUtil.validateHasSingleSample" should "assert that a single sample is present in the input VCF source" in {
    val singleSample = VcfBuilder(samples = Seq(sample1)).toSource
    val doubleSample = VcfBuilder(samples = Seq(sample1, sample2)).toSource
    noException shouldBe thrownBy { VcfUtil.validateHasSingleSample(singleSample) }
    an[IllegalArgumentException] shouldBe thrownBy { VcfUtil.validateHasSingleSample(doubleSample) }
    singleSample.close()
  }
}
