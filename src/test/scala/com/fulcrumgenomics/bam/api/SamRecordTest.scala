/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.{SAMRecord, TextCigarCodec}

class SamRecordTest extends UnitSpec {
  "SamRecord.isFrPair" should "return false for reads that are not part of a mapped pair" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addFrag(start=100).exists(_.isFrPair) shouldBe false
    builder.addPair(start1=100, start2=100, unmapped2=true).exists(_.isFrPair) shouldBe false
  }

  it should "return false for pairs that are mapped to different chromosomes" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200).foreach { rec =>
      rec.mateRefIndex = rec.refIndex + 1
      rec.isFrPair shouldBe false
    }
  }

  it should "return false for mapped FF, RR and RF pairs" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200, strand1=Plus,  strand2=Plus).exists(_.isFrPair)  shouldBe false
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Minus).exists(_.isFrPair) shouldBe false
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Plus).exists(_.isFrPair)  shouldBe false
  }

  it should "return true on an actual FR pair" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addPair(start1=100, start2=200, strand1=Plus, strand2=Minus).forall(_.isFrPair)  shouldBe true
  }

  "SamRecord.cigar" should "always return the latest/correct cigar, when when SAMRecord methods are used" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10, cigar="50M").get

    rec.cigar.toString shouldBe          "50M"
    rec.asSam.getCigarString shouldBe    "50M"
    rec.asSam.getCigar.toString shouldBe "50M"

    // Test with setting via SamRecord.cigar =
    rec.cigar = Cigar("10S40M")
    rec.cigar.toString shouldBe          "10S40M"
    rec.asSam.getCigarString shouldBe    "10S40M"
    rec.asSam.getCigar.toString shouldBe "10S40M"

    // Test with setting via SAMRecord.setCigarString
    rec.asSam.setCigarString("25M2D25M")
    rec.cigar.toString shouldBe          "25M2D25M"
    rec.asSam.getCigarString shouldBe    "25M2D25M"
    rec.asSam.getCigar.toString shouldBe "25M2D25M"

    // Test with setting via SAMRecord.setCigar()
    rec.asSam.setCigar(TextCigarCodec.decode("10H50M"))
    rec.cigar.toString shouldBe          "10H50M"
    rec.asSam.getCigarString shouldBe    "10H50M"
    rec.asSam.getCigar.toString shouldBe "10H50M"
  }

  it should "work with empty cigars" in {
    val builder = new SamBuilder(readLength=50)
    val rec = builder.addFrag(start=10, unmapped=true, cigar="50M").get

    rec.cigar = Cigar("50M")
    rec.cigar = Cigar.empty
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null

    rec.cigar = Cigar("50M")
    rec.asSam.setCigar(null)
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null

    rec.cigar = Cigar("50M")
    rec.asSam.setCigarString(null)
    rec.cigar.toString shouldBe       ""
    rec.asSam.getCigarString shouldBe null
    rec.asSam.getCigar shouldBe       null
  }
}
