/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

class GenotypeMapTest extends UnitSpec {
  private val Alleles = AlleleSet("A", "C", "G", "T")
  private val A = Alleles(0)
  private val C = Alleles(1)
  private val G = Alleles(2)
  private val T = Alleles(3)

  private val Genotypes = Array(
      Genotype(Alleles, "x", IndexedSeq(A, C)),
      Genotype(Alleles, "y", IndexedSeq(G, T)),
      Genotype(Alleles, "z", IndexedSeq(A, G)),
      Genotype(Alleles, "a", IndexedSeq(A, A)),
      Genotype(Alleles, "1", IndexedSeq(G, G))
    )

  "GenotypeMap" should "maintain ordering" in {
    val map = new GenotypeMap(Genotypes, Genotypes.zipWithIndex.map { case (g, i) => g.sample -> i }.toMap)
    map.keysIterator.toIndexedSeq   should contain theSameElementsInOrderAs Seq("x", "y", "z", "a", "1")
    map.valuesIterator.toIndexedSeq should contain theSameElementsInOrderAs Genotypes
    map.iterator.toIndexedSeq       should contain theSameElementsInOrderAs Genotypes.map(g => (g.sample, g))
  }

  it should "support generating new maps with removed and updated elements" in {
    val map = new GenotypeMap(Genotypes, Genotypes.zipWithIndex.map { case (g, i) => g.sample -> i }.toMap)

    val rem = map.removed("a")
    rem.keysIterator.toIndexedSeq should contain theSameElementsInOrderAs IndexedSeq("x", "y", "z", "1")
    rem.iterator.foreach { case (s, gt) => s shouldBe gt.sample }

    val upd = map.updated("a", Genotype(Alleles, "a", IndexedSeq(A, C)))
    upd.foreach { case (sample, gt) =>
      sample shouldBe gt.sample
      if (sample == "a") gt.gtWithBases shouldBe "A/C"
      else gt shouldBe map(sample)
    }
  }
}
