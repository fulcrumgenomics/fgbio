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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import org.scalatest.OptionValues

class VariantTest extends UnitSpec with OptionValues {
  "Variant.isMissingValue" should "correctly handle missing and non-missing values" in {
    Variant.isMissingValue(".") shouldBe true
    Variant.isMissingValue('.') shouldBe true
    Variant.isMissingValue(Variant.MissingInt) shouldBe true
    Variant.isMissingValue(Variant.MissingFloat) shouldBe true

    Variant.isMissingValue(".1")      shouldBe false
    Variant.isMissingValue(Float.NaN) shouldBe false
    Range(-500, 500).foreach { i => Variant.isMissingValue(i) shouldBe false}
  }
}
