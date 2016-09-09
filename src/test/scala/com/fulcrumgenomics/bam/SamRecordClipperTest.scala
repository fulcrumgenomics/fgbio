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
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.bam.SamRecordClipper.ClippingMode
import com.fulcrumgenomics.testing.SamRecordSetBuilder._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.{SAMRecord, TextCigarCodec}

class SamRecordClipperTest extends UnitSpec {
  /** Returns a fragment SAM record with the start / cigar / strand requested. */
  def r(start: Int, cigar: String, strand: Strand = Plus): SAMRecord = {
    val cig = TextCigarCodec.decode(cigar)
    val len = cig.getReadLength
    new SamRecordSetBuilder(readLength=len).addFrag(start=start, cigar = cigar, strand=strand).get
  }

  "SamRecordClipper.clipStartOfAlignment" should "soft-clip 10 matched bases" in {
    val rec = r(10, "50M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "10S40M"
  }

  it should "soft-clip 10 matched/inserted bases" in {
    val rec = r(10, "4M2I44M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 18
    rec.getCigarString shouldBe "10S40M"
  }

  it should "soft-clip 10 matched/deleted bases" in {
    val rec = r(10, "6M2D44M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 22
    rec.getCigarString shouldBe "10S40M"
  }

  it should "soft-clip 10 additional bases" in {
    val rec = r(10, "10S40M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "20S30M"
  }

  it should "preserve hard-clipping while soft-clipping" in {
    val rec = r(10, "10H40M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "10H10S30M"
  }

  it should "handle a cigar with a ton of elements in it while soft-clipping" in {
    val rec = r(10, "2H4S16M10I5M5I10M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "2H14S6M10I5M5I10M"
  }

  it should "also consume any extra bases if it ends in an insertion when soft-clipping" in {
    val rec = r(10, "8M4I38M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 18
    rec.getCigarString shouldBe "12S38M"
  }

  it should "preserve an insertion that starts immediately after the soft-clipping" in {
    val rec = r(10, "10M4I36M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "10S4I36M"
  }

  it should "remove deletions that immediately follow the clipping" in {
    val rec = r(10, "10M4D40M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 24
    rec.getCigarString shouldBe "10S40M"
  }

  it should "NOT throw an exception, but do nothing, with an unmapped read" in {
    val rec = r(10, "50M")
    rec.setReadUnmappedFlag(true)
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask" in {
    val rec = r(10, "50M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.SoftWithMask)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "10S40M"
    rec.getReadBases.take(10).forall(_ == 'N'.toByte) shouldBe true
    rec.getBaseQualities.take(10).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask, including existing soft-clips" in {
    val rec = r(10, "10S40M")
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.SoftWithMask)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "20S30M"
    rec.getReadBases.take(20).forall(_ == 'N'.toByte) shouldBe true
    rec.getBaseQualities.take(20).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "hard clip 10 bases" in {
    val rec = r(10, "50M")
    val bases = rec.getReadBases.clone()
    val quals = rec.getBaseQualities.clone()
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Hard)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "10H40M"
    rec.getReadBases shouldBe bases.drop(10)
    rec.getBaseQualities shouldBe quals.drop(10)
  }

  it should "hard clip 10 extra bases and convert existing soft-clipping to hard clipping" in {
    val rec = r(10, "10S40M")
    val bases = rec.getReadBases.clone()
    val quals = rec.getBaseQualities.clone()
    SamRecordClipper.clipStartOfAlignment(rec, 10, ClippingMode.Hard)
    rec.getAlignmentStart shouldBe 20
    rec.getCigarString shouldBe "20H30M"
    rec.getReadBases shouldBe bases.drop(20)
    rec.getBaseQualities shouldBe quals.drop(20)
  }

  it should "unmap a read when more clipping is requested than there is alignment to clip" in {
    val rec = r(10, "10S40M")
    SamRecordClipper.clipStartOfAlignment(rec, 50, ClippingMode.Soft)
    rec.getReadUnmappedFlag shouldBe true
    rec.getCigar shouldBe 'empty
  }

  //////////////////////////////////////////////////////////////////////////////

  "SamRecordClipper.clipEndOfAlignment" should "soft-clip 10 matched bases" in {
    val rec = r(10, "50M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10S"
  }

  it should "soft-clip 10 matched/inserted bases" in {
    val rec = r(10, "44M2I4M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10S"
  }

  it should "soft-clip 10 matched/deleted bases" in {
    val rec = r(10, "44M2D6M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10S"
  }

  it should "soft-clip 10 additional bases" in {
    val rec = r(10, "40M10S")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "30M20S"
  }

  it should "preserve hard-clipping while soft-clipping" in {
    val rec = r(10, "40M10H")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "30M10S10H"
  }

  it should "handle a cigar with a ton of elements in it while soft-clipping" in {
    val rec = r(10, "10M5I5M10I16M4S2H")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "10M5I5M10I6M14S2H"
  }

  it should "also consume any extra bases if it ends in an insertion when soft-clipping" in {
    val rec = r(10, "38M4I8M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "38M12S"
  }

  it should "preserve an insertion that starts immediately after the soft-clipping" in {
    val rec = r(10, "36M4I10M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "36M4I10S"
  }

  it should "remove deletions that immediately follow the clipping" in {
    val rec = r(10, "40M4D10M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10S"
  }

  it should "NOT throw an exception, but do nothing, with an unmapped read" in {
    val rec = r(10, "50M")
    rec.setReadUnmappedFlag(true)
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Soft)
    rec.getAlignmentStart shouldBe 10
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask" in {
    val rec = r(10, "50M")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.SoftWithMask)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10S"
    rec.getReadBases.drop(40).forall(_ == 'N'.toByte) shouldBe true
    rec.getBaseQualities.drop(40).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask, including existing soft-clips" in {
    val rec = r(10, "40M10S")
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.SoftWithMask)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "30M20S"
    rec.getReadBases.drop(30).forall(_ == 'N'.toByte) shouldBe true
    rec.getBaseQualities.drop(30).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "hard clip 10 bases" in {
    val rec = r(10, "50M")
    val bases = rec.getReadBases.clone()
    val quals = rec.getBaseQualities.clone()
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Hard)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "40M10H"
    rec.getReadBases shouldBe bases.take(40)
    rec.getBaseQualities shouldBe quals.take(40)
  }

  it should "hard clip 10 extra bases and convert existing soft-clipping to hard clipping" in {
    val rec = r(10, "40M10S")
    val bases = rec.getReadBases.clone()
    val quals = rec.getBaseQualities.clone()
    SamRecordClipper.clipEndOfAlignment(rec, 10, ClippingMode.Hard)
    rec.getAlignmentStart shouldBe 10
    rec.getCigarString shouldBe "30M20H"
    rec.getReadBases shouldBe bases.take(30)
    rec.getBaseQualities shouldBe quals.take(30)
  }

    it should "unmap a read when more clipping is requested than there is alignment to clip" in {
    val rec = r(10, "40M10S")
    SamRecordClipper.clipEndOfAlignment(rec, 50, ClippingMode.Soft)
    rec.getReadUnmappedFlag shouldBe true
    rec.getCigar shouldBe 'empty
  }


  //////////////////////////////////////////////////////////////////////////////
  // One test each for all the derived methods.
  //////////////////////////////////////////////////////////////////////////////

  "SamRecordClipper.clipStartOfRead" should "expand existing clipping" in {
    val rec1 = r(10, "50M")
    SamRecordClipper.clipStartOfRead(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "10S40M"

    val rec2 = r(10, "5S45M")
    SamRecordClipper.clipStartOfRead(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "10S40M"

    val rec3 = r(10, "20S30M")
    SamRecordClipper.clipStartOfRead(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "20S30M"

    val rec4 = r(10, "20H30M")
    SamRecordClipper.clipStartOfRead(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "20H30M"
  }

  "SamRecordClipper.clipEndOfRead" should "expand existing clipping" in {
    val rec1 = r(10, "50M")
    SamRecordClipper.clipEndOfRead(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "40M10S"

    val rec2 = r(10, "45M5S")
    SamRecordClipper.clipEndOfRead(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "40M10S"

    val rec3 = r(10, "30M20S")
    SamRecordClipper.clipEndOfRead(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "30M20S"

    val rec4 = r(10, "30M20H")
    SamRecordClipper.clipEndOfRead(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "30M20H"
  }

  "SamRecordClipper.clip5PrimeEndOfAlignment" should "add more clipping to the 5' end" in {
    val rec1 = r(10, "50M", strand=Plus)
    SamRecordClipper.clip5PrimeEndOfAlignment(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Plus)
    SamRecordClipper.clip5PrimeEndOfAlignment(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "20S30M"

    val rec3 = r(10, "50M", strand=Minus)
    SamRecordClipper.clip5PrimeEndOfAlignment(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Minus)
    SamRecordClipper.clip5PrimeEndOfAlignment(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "30M20S"
  }

  "SamRecordClipper.clip3PrimeEndOfAlignment" should "add more clipping to the 3' end" in {
    val rec1 = r(10, "50M", strand=Minus)
    SamRecordClipper.clip3PrimeEndOfAlignment(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Minus)
    SamRecordClipper.clip3PrimeEndOfAlignment(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "20S30M"

    val rec3 = r(10, "50M", strand=Plus)
    SamRecordClipper.clip3PrimeEndOfAlignment(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Plus)
    SamRecordClipper.clip3PrimeEndOfAlignment(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "30M20S"
  }

  "SamRecordClipper.clip5PrimeEndOfRead" should "expand existing clipping at the 5' end" in {
    val rec1 = r(10, "50M", strand=Plus)
    SamRecordClipper.clip5PrimeEndOfRead(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Plus)
    SamRecordClipper.clip5PrimeEndOfRead(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "10S40M"

    val rec3 = r(10, "50M", strand=Minus)
    SamRecordClipper.clip5PrimeEndOfRead(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Minus)
    SamRecordClipper.clip5PrimeEndOfRead(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "40M10S"
  }

  "SamRecordClipper.clip3PrimeEndOfRead" should "expand existing clipping at the 3' end" in {
    val rec1 = r(10, "50M", strand=Minus)
    SamRecordClipper.clip3PrimeEndOfRead(rec1, 10, ClippingMode.Soft)
    rec1.getCigarString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Minus)
    SamRecordClipper.clip3PrimeEndOfRead(rec2, 10, ClippingMode.Soft)
    rec2.getCigarString shouldBe "10S40M"

    val rec3 = r(10, "50M", strand=Plus)
    SamRecordClipper.clip3PrimeEndOfRead(rec3, 10, ClippingMode.Soft)
    rec3.getCigarString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Plus)
    SamRecordClipper.clip3PrimeEndOfRead(rec4, 10, ClippingMode.Soft)
    rec4.getCigarString shouldBe "40M10S"
  }
}
