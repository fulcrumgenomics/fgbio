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
import com.fulcrumgenomics.vcf.api.Allele.{NoCallAllele, SimpleAllele, SpannedAllele, SymbolicAllele}

class AlleleTest extends UnitSpec {
  "Allele.apply()" should "generate SimpleAlleles for base strings" in {
    Allele("ACT").isInstanceOf[SimpleAllele] shouldBe true
    Allele("ACT").value shouldBe "ACT"
    Allele("G").value shouldBe "G"
    Allele("g").value shouldBe "g"
  }

  it should "return cached instances for single-base alleles" in {
    Allele("A") eq Allele("a") shouldBe false
    "ACGTNacgtn".sliding(1).foreach(b => Allele(b) eq Allele(b) shouldBe true)
  }

  it should "return the constant no-call allele for no-calls" in {
    Allele(".") shouldBe NoCallAllele
    Allele(".") eq NoCallAllele shouldBe true
  }

  it should "return the constant spanned allele for *" in {
    Allele("*") eq SpannedAllele shouldBe true
  }

  it should "construct symbolic alleles for string wrapped in <>" in {
    Allele("<NON_REF>").isInstanceOf[SymbolicAllele] shouldBe true
    Allele("<NON_REF>").value shouldBe "<NON_REF>"
  }

  it should "reject alleles with non-ACGTN bases in them" in {
    an[IllegalArgumentException] shouldBe thrownBy { Allele("ARR") }
    an[IllegalArgumentException] shouldBe thrownBy { Allele("ABC") }
    an[IllegalArgumentException] shouldBe thrownBy { Allele("Y") }
  }

  "Allele equality" should "be case insensitive" in {
    Allele("T") == Allele("A") shouldBe false

    "ACGTN".sliding(1).foreach { s =>
      val upper = Allele(s)
      val lower = Allele(s.toLowerCase)
      upper.hashCode() shouldBe lower.hashCode()
      upper == lower shouldBe true
    }
  }
}
