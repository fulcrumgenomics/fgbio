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
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import org.scalatest.OptionValues

class GenotypeTest extends UnitSpec with OptionValues {

  private def gt(ref: String, gt: String, sample: String = "s1", attrs: Map[String,Any] = Map.empty): Genotype = {
    val phased    = gt.contains('|')
    val parts     = if (phased) gt.split('|') else gt.split('/')
    val gtAlleles = parts.map(s => Allele(s)).toIndexedSeq
    val refAllele = Allele(ref)
    val alleles   = AlleleSet(ref=refAllele, alts=gtAlleles.filterNot(a => a == refAllele || a == NoCallAllele))
    Genotype(alleles, sample, gtAlleles, phased=phased)
  }

  "Genotype.ploidy" should "return the right number for various genotypes" in {
    gt(ref="A", gt="A").ploidy shouldBe 1
    gt(ref="A", gt="A/A").ploidy shouldBe 2
    gt(ref="A", gt="A/C").ploidy shouldBe 2
    gt(ref="A", gt="C/T").ploidy shouldBe 2
    gt(ref="A", gt="A|C|G").ploidy shouldBe 3
  }

  "Genotype.isNoCall" should "be true if all calls are no-calls and false otherwise" in {
    gt(ref="C", gt="C/A").isNoCall shouldBe false
    gt(ref="C", gt="C/.").isNoCall shouldBe false
    gt(ref="C", gt="C/*").isNoCall shouldBe false
    gt(ref="C", gt=".").isNoCall   shouldBe true
    gt(ref="C", gt="./.").isNoCall shouldBe true
  }

  "Genotype.isFullyCalled" should "return true only if no called allele is a no-call" in {
    gt(ref="C", gt="C/A").isFullyCalled shouldBe true
    gt(ref="C", gt="C/.").isFullyCalled shouldBe false
    gt(ref="C", gt="C/*").isFullyCalled shouldBe true
    gt(ref="C", gt=".").isFullyCalled   shouldBe false
    gt(ref="C", gt="./.").isFullyCalled shouldBe false
  }

  "Genotype.isHomRef" should "return true only if the call is homozygous reference" in {
    gt(ref="C", gt="C/A").isHomRef shouldBe false
    gt(ref="C", gt="A").isHomRef   shouldBe false
    gt(ref="C", gt="C/*").isHomRef shouldBe false
    gt(ref="C", gt="C/.").isHomRef shouldBe false
    gt(ref="C", gt="./.").isHomRef shouldBe false
    gt(ref="C", gt="A/A").isHomRef shouldBe false
    gt(ref="C", gt="C").isHomRef   shouldBe true
    gt(ref="C", gt="C/C").isHomRef shouldBe true
    gt(ref="C", gt="C/C/C/C").isHomRef shouldBe true
  }
  
  "Genotype.isHomVar" should "return true only if the call is homozygous alt" in {
    gt(ref="C", gt="C/A").isHomVar   shouldBe false
    gt(ref="C", gt="A"  ).isHomVar   shouldBe true
    gt(ref="C", gt="A/*").isHomVar   shouldBe false
    gt(ref="C", gt="A/.").isHomVar   shouldBe false
    gt(ref="C", gt="./.").isHomVar   shouldBe false
    gt(ref="C", gt="A/A").isHomVar   shouldBe true
    gt(ref="C", gt="."  ).isHomVar   shouldBe false
    gt(ref="C", gt="C/C").isHomVar   shouldBe false
    gt(ref="C", gt="A/A/A").isHomVar shouldBe true
  }

  "Genotype.isHet" should "return true if the genotype contains any heterozygosity" in {
    gt(ref="C", gt="C/C").isHet   shouldBe false
    gt(ref="C", gt="A/A").isHet   shouldBe false
    gt(ref="C", gt="A/*").isHet   shouldBe false
    gt(ref="C", gt="A/.").isHet   shouldBe false
    gt(ref="C", gt="C/A").isHet   shouldBe true
    gt(ref="C", gt="G/T").isHet   shouldBe true
    gt(ref="C", gt="A/C").isHet   shouldBe true
    gt(ref="C", gt="A/C/*").isHet shouldBe true
    gt(ref="C", gt="A/C/.").isHet shouldBe true
  }

  "Genotype.isHetNonRef" should "return true if the genotype contains any heterozygosity and doesn't contain the reference allele" in {
    gt(ref="C", gt="C/C").isHetNonRef   shouldBe false
    gt(ref="C", gt="A/A").isHetNonRef   shouldBe false
    gt(ref="C", gt="A/*").isHetNonRef   shouldBe false
    gt(ref="C", gt="A/.").isHetNonRef   shouldBe false
    gt(ref="C", gt="C/A").isHetNonRef   shouldBe false
    gt(ref="C", gt="G/T").isHetNonRef   shouldBe true
    gt(ref="C", gt="A/C").isHetNonRef   shouldBe false
    gt(ref="C", gt="A/G/*").isHetNonRef shouldBe true
    gt(ref="C", gt="A/T/.").isHetNonRef shouldBe true
  }

  "Genotype.gtWithBases" should "return a VCF-like genotype string using the actual allele strings" in {
    gt(ref="C", gt="C/C").gtWithBases       shouldBe "C/C"
    gt(ref="C", gt="A/T").gtWithBases       shouldBe "A/T"
    gt(ref="C", gt="A/*").gtWithBases       shouldBe "A/*"
    gt(ref="C", gt="T|*").gtWithBases       shouldBe "T|*"
    gt(ref="C", gt="A/C/G").gtWithBases     shouldBe "A/C/G"
    gt(ref="C", gt="A|.|G").gtWithBases     shouldBe "A|.|G"
    gt(ref="C", gt="A|<FOO>|G").gtWithBases shouldBe "A|<FOO>|G"
  }

  "Genotype.gtWithIndices" should "return a VCF-like genotype string using the indices of the alleles" in {
    gt(ref="C", gt="C/C").gtVcfStyle       shouldBe "0/0"
    gt(ref="C", gt="A/T").gtVcfStyle       shouldBe "1/2"
    gt(ref="C", gt="A/*").gtVcfStyle       shouldBe "1/2"
    gt(ref="C", gt="T|*").gtVcfStyle       shouldBe "1|2"
    gt(ref="C", gt="A/C/G").gtVcfStyle     shouldBe "1/0/2"
    gt(ref="C", gt="A|.|G").gtVcfStyle     shouldBe "1|.|2"
    gt(ref="C", gt="A|<FOO>|G").gtVcfStyle shouldBe "1|2|3"
  }
}
