/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

import scala.util.Random

class PrefixTrieSetTest extends UnitSpec {

  "PrefixTrieSet" should "support an empty string" in {
    val set = new PrefixTrieSet()
    set.contains("") shouldBe false
    set.add("") shouldBe 1
    set.contains("") shouldBe true
    set.add("") shouldBe 2
    set.prefixes().toSeq should contain theSameElementsAs Seq("")
  }

  it should "support strings with no common prefix" in {
    val set = new PrefixTrieSet()
    "ACGT".foreach { prefix =>
      val key = prefix + "GATTACA"
      set.contains(key) shouldBe false
      set.add(key) shouldBe 1
      set.contains(key) shouldBe true
      set.add(key) shouldBe 2
      set.contains(key) shouldBe true
    }
    set.prefixes().toSeq should contain theSameElementsAs "ACGT".map(c => c + "GATTACA")
  }

  it should "collapse prefixes" in {
    val set = new PrefixTrieSet()
    val keys = Seq(
      "GATTACA",
      "GATACA",
      "CATTACA",
      "CATTACA",
      "GAT",
      "",
      "GATTACA"
    )
    keys.foreach(set.add)
    keys.forall(set.contains) shouldBe true
    set.prefixes().toSet should contain theSameElementsAs keys.toSet
    set.values().toSeq should contain theSameElementsAs Seq(
      ("GATTACA", 2),
      ("GATACA", 1),
      ("CATTACA", 2),
      ("GAT", 1),
      ("", 1),
    )
  }

  Seq(0, 1, 10, 100, 1000, 10000, 25000).foreach { numSequences =>
    Seq(1, 2, 4, 8, 16, 32, 64).foreach { length =>
      it should f"support adding in $numSequences random sequences of maximum length $length" in {
        val bases = "ACGT"
        val trieSet = new PrefixTrieSet()
        val hashSet = new scala.collection.mutable.HashSet[String]()
        val minLength = length / 2
        Range.inclusive(start=1, end=numSequences).foreach { _ =>
          val keyLength = Random.between(minLength, length)
          val key: String = Range.inclusive(start=1, end=keyLength).map { _ =>
            bases.charAt(Random.nextInt(4))
          }.toString
          trieSet.add(key)
          hashSet.add(key)
        }
        trieSet.prefixes().toSeq should contain theSameElementsAs hashSet.toSeq
        trieSet.values().map(_._2).sum shouldBe numSequences
      }
    }
  }
}
