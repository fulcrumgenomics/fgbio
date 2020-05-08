/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.testing.UnitSpec
import org.scalatest.OptionValues

class SequenceDictionaryTest extends UnitSpec with OptionValues {

  "SequenceMetadata" should "fail if not a validated sequence" in {
    an[Exception] should be thrownBy SequenceMetadata(name="(((242*(@&$)))", length=12)
  }

  it should "fail if the length is less than zero" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=(-1))
  }

  it should "fail if the sequence name tag is given in the attributes" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=1, attributes=Map("SN" -> "chr1"))
  }

  it should "fail if the sequence length tag is given in the attributes" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=1, attributes=Map("LN" -> "2"))
  }

  it should "implement apply, get, and contains" in {
    val sequence = SequenceMetadata(
      name = "chr1",
      length = 100,
      attributes = Map(
        "key1" -> "value1",
        "key2" -> "value2",
      )
    )

    // apply
    sequence("key1") shouldBe "value1"
    an[Exception] should be thrownBy sequence("key3")

    // get
    sequence.get("key1").value shouldBe "value1"
    sequence.get("key3") shouldBe 'empty

    // contains
    sequence.contains("key1") shouldBe true
    sequence.contains("key3") shouldBe false
  }

  it should "return all aliases if present in AN" in {
    val sequence = SequenceMetadata(
      name = "chr1",
      length = 100,
      attributes = Map(
        "AN" -> "chr1_1,chr1_2,chr1_3"
      )
    )
    sequence.allNames should contain theSameElementsInOrderAs Seq("chr1", "chr1_1", "chr1_2", "chr1_3")
    sequence.aliases should contain theSameElementsInOrderAs Seq("chr1_1", "chr1_2", "chr1_3")
  }

  it should "return if the sequence is an alternate locus" in {
    val noAltLocus  = SequenceMetadata(name="chr1", length=123)
    val altStar     = SequenceMetadata(name="chr1", length=123, attributes=Map("AH" -> "*"))
    val locusEquals = SequenceMetadata(name="chr1", length=123, attributes=Map("AH" -> "=:1-2"))
    val locusFull   = SequenceMetadata(name="chr1", length=123, attributes=Map("AH" -> "chr4:2-3"))

    noAltLocus.isAlternate shouldBe false
    noAltLocus.alternate shouldBe 'empty

    altStar.isAlternate shouldBe false
    altStar.alternate shouldBe 'empty

    locusEquals.isAlternate shouldBe true
    locusEquals.alternate.value shouldBe AlternateLocus(refName="chr1", start=1, end=2)

    locusFull.isAlternate shouldBe true
    locusFull.alternate.value shouldBe AlternateLocus(refName="chr4", start=2, end=3)
  }

  // md5

  // assembly

  // species

  // description

  // topology

  // same

  // SequenceDictionary

  // Converters
}
