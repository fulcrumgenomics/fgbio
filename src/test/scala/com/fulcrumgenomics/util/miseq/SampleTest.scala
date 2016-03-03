/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util.miseq

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.ReadStructure
import org.scalatest.OptionValues

/**
  * Tests for Sample.
  */
class SampleTest extends UnitSpec with OptionValues {

  "Sample.sampleBarcodeBases" should "return the sample barcodes if present" in {
    new Sample(0, "ID", "NAME").sampleBarcodeBases.flatten shouldBe 'empty
    new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA")).sampleBarcodeBases.flatten should contain theSameElementsInOrderAs Seq("GATTACA")
    new Sample(0, "ID", "NAME", i5IndexBases=Some("GATTACA")).sampleBarcodeBases.flatten should contain theSameElementsInOrderAs Seq("GATTACA")
    new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA"), i5IndexBases=Some("TGTAATC")).sampleBarcodeBases.flatten should contain theSameElementsInOrderAs Seq("GATTACA", "TGTAATC")
  }

  it should "support extended attributes" in {
    val sample = new Sample(0, "ID", "NAME", extendedAttributes=Seq(("FOO", "1"), ("BAR", "2  ")).toMap)
    sample.extendedAttribute("foo").value shouldBe "1"
    sample.extendedAttribute("Foo").value shouldBe "1"
    sample.extendedAttribute("FOO").value shouldBe "1"
    sample.extendedAttribute("bar").value shouldBe "2"
    sample.extendedAttribute("Bar").value shouldBe "2"
    sample.extendedAttribute("BAR").value shouldBe "2"
    sample.extendedAttribute("car") shouldBe 'empty
    sample.extendedAttribute("Car") shouldBe 'empty
    sample.extendedAttribute("CAR") shouldBe 'empty
  }

  it should "support throw an exception if extended attribute keys are not all uppercase" in {
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", extendedAttributes=Seq(("foo", "1"), ("bar", "2  ")).toMap)
  }

  "Sample.setSampleBarcode" should "set the sample barcode" in {
    // no sample barcodes
    new Sample(0, "ID", "NAME")
      .setSampleBarcode(Seq.empty)
      .sampleBarcode shouldBe 'empty
    // one sample barcode
    new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA"))
      .setSampleBarcode(Seq(ReadStructure("7B"))).sampleBarcode.value.toString shouldBe "GATTACA"
    // two sample barcodes
    new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA"), i5IndexBases=Some("TGTAATC"))
      .setSampleBarcode(Seq(ReadStructure("7B"), ReadStructure("7B")))
      .sampleBarcode.value.toString shouldBe "GATTACA-TGTAATC"
    // one sample barcode but two segements in the same read
    new Sample(0, "ID", "NAME", i7IndexBases=Some("GATACA"))
      .setSampleBarcode(Seq(ReadStructure("3B3B"))).sampleBarcode.value.toString shouldBe "GATACA"
  }

  it should "throw an exception if a read structure is given that contains a segment that is not a sample barcode" in {
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACAA"))
      .setSampleBarcode(Seq(ReadStructure("7B7S")))
  }

  it should "throw an exception if the # of read structures was different than the number of reads with sample barcodes" in {
    // too many
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA"))
      .setSampleBarcode(Seq(ReadStructure("7B"), ReadStructure("7B")))
    // too few
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACA"), i5IndexBases=Some("TGTAATC"))
      .setSampleBarcode(Seq(ReadStructure("7B")))
  }

  it should "throw an exception if there are a different # of sample barcode bases than in the read structure" in {
    // too many
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTACAA"))
      .setSampleBarcode(Seq(ReadStructure("7B")))
    // too few
    an[Exception] should be thrownBy new Sample(0, "ID", "NAME", i7IndexBases=Some("GATTAC"))
      .setSampleBarcode(Seq(ReadStructure("7B")))
  }
}
