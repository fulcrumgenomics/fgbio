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
import com.fulcrumgenomics.vcf.api.Allele.{NoCallAllele, SpannedAllele}

class AlleleSetTest extends UnitSpec {
  "AlleleSet" should "perform basic operations correctly" in {
    val as = AlleleSet(ref=Allele("A"), alts=Seq(Allele("T"), Allele("G")))
    as.size shouldBe 3
    as.indexOf(Allele("A")) shouldBe 0
    as.indexOf(Allele("T")) shouldBe 1
    as.indexOf(Allele("G")) shouldBe 2
    as.toSeq shouldBe Seq(Allele("A"), Allele("T"), Allele("G"))
    as(0) shouldBe Allele("A")
    as(1) shouldBe Allele("T")
    as(2) shouldBe Allele("G")
  }

  it should "reject non-simple alleles as the reference" in {
    an[IllegalArgumentException] shouldBe thrownBy { AlleleSet(ref=NoCallAllele) }
    an[IllegalArgumentException] shouldBe thrownBy { AlleleSet(ref=SpannedAllele) }
    an[IllegalArgumentException] shouldBe thrownBy { AlleleSet(ref=Allele("<NON_REF>")) }
  }
}
