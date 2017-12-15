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

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.ClippingMode._
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus, Strand}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.TextCigarCodec

class SamRecordClipperTest extends UnitSpec {
  /** Returns a fragment SAM record with the start / cigar / strand requested. */
  def r(start: Int, cigar: String, strand: Strand = Plus, attrs:Map[String,Any] = Map.empty): SamRecord = {
    val cig = TextCigarCodec.decode(cigar)
    val len = cig.getReadLength
    new SamBuilder(readLength=len).addFrag(start=start, cigar = cigar, strand=strand, attrs=attrs).get
  }

  def clipper(mode: ClippingMode, autoClip: Boolean=false) = new SamRecordClipper(mode, autoClip)

  "SamRecordClipper.clipStartOfAlignment" should "soft-clip 10 matched bases" in {
    val rec = r(10, "50M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10S40M"
  }

  it should "soft-clip 10 matched/inserted bases" in {
    val rec = r(10, "4M2I44M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 18
    rec.cigar.toString shouldBe "10S40M"
  }

  it should "soft-clip 10 matched/deleted bases" in {
    val rec = r(10, "6M2D44M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 22
    rec.cigar.toString shouldBe "10S40M"
  }

  it should "soft-clip 10 additional bases" in {
    val rec = r(10, "10S40M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "20S30M"
  }

  it should "preserve hard-clipping while soft-clipping" in {
    val rec = r(10, "10H40M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10H10S30M"
  }

  it should "handle a cigar with a ton of elements in it while soft-clipping" in {
    val rec = r(10, "2H4S16M10I5M5I10M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "2H14S6M10I5M5I10M"
  }

  it should "also consume any extra bases if it ends in an insertion when soft-clipping" in {
    val rec = r(10, "8M4I38M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 12 // 8 + 4
    rec.start shouldBe 18
    rec.cigar.toString shouldBe "12S38M"
  }

  it should "preserve an insertion that starts immediately after the soft-clipping" in {
    val rec = r(10, "10M4I36M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10S4I36M"
  }

  it should "remove deletions that immediately follow the clipping" in {
    val rec = r(10, "10M4D40M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 24
    rec.cigar.toString shouldBe "10S40M"
  }

  it should "preseve deletions that are not close to the clipped region" in {
    val rec = r(10, "25M4D25M")
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10S15M4D25M"
  }

  it should "NOT throw an exception, but do nothing, with an unmapped read" in {
    val rec = r(10, "50M")
    rec.unmapped = true
    clipper(Soft).clipStartOfAlignment(rec, 10) shouldBe 0
    rec.start shouldBe 10
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask" in {
    val rec = r(10, "50M")
    clipper(SoftWithMask).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10S40M"
    rec.bases.take(10).forall(_ == 'N'.toByte) shouldBe true
    rec.quals.take(10).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask, including existing soft-clips" in {
    val rec = r(10, "10S40M")
    clipper(SoftWithMask).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "20S30M"
    rec.bases.take(20).forall(_ == 'N'.toByte) shouldBe true
    rec.quals.take(20).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "hard clip 10 bases" in {
    val rec = r(10, "50M")
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    clipper(Hard).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "10H40M"
    rec.bases shouldBe bases.drop(10)
    rec.quals shouldBe quals.drop(10)
  }

  it should "hard clip 10 extra bases and convert existing soft-clipping to hard clipping" in {
    val rec = r(10, "10S40M")
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    clipper(Hard).clipStartOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "20H30M"
    rec.bases shouldBe bases.drop(20)
    rec.quals shouldBe quals.drop(20)
  }

  it should "unmap a read when more clipping is requested than there is alignment to clip" in {
    val rec = r(10, "10S40M")
    clipper(Soft).clipStartOfAlignment(rec, 50) shouldBe 40 // 40M
    rec.unmapped shouldBe true
    rec.cigar shouldBe Cigar.empty
  }

  for (mode <- ClippingMode.values; auto <- Seq(true, false)) {
    it should s"correctly handle auto-trimming of attribute with auto=$auto and mode=$mode" in {
      val rec = r(10, "20M")
      rec("A1") =  "AB" * 10
      rec("A2") =  Range.inclusive(1, 20).toArray
      rec("B1") =  "A" * 10
      rec("B2") =  Range.inclusive(1, 10).toArray

      val cut = mode == ClippingMode.Hard && auto
      clipper(mode=mode, autoClip=auto).clipStartOfAlignment(rec, 5) shouldBe 5
      rec[String]("A1")     shouldBe { if (cut) "BABABABABABABAB" else "AB"*10 }
      rec[Array[Int]]("A2") shouldBe Range.inclusive(if (cut) 6 else 1, 20).toArray
      rec[String]("B1")     shouldBe  "A" * 10
      rec[Array[Int]]("B2") shouldBe Range.inclusive(1, 10).toArray
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  "SamRecordClipper.clipEndOfAlignment" should "soft-clip 10 matched bases" in {
    val rec = r(10, "50M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10S"
  }

  it should "soft-clip 10 matched/inserted bases" in {
    val rec = r(10, "44M2I4M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10S"
  }

  it should "soft-clip 10 matched/deleted bases" in {
    val rec = r(10, "44M2D6M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10S"
  }

  it should "soft-clip 10 additional bases" in {
    val rec = r(10, "40M10S")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "30M20S"
  }

  it should "preserve hard-clipping while soft-clipping" in {
    val rec = r(10, "40M10H")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "30M10S10H"
  }

  it should "handle a cigar with a ton of elements in it while soft-clipping" in {
    val rec = r(10, "10M5I5M10I16M4S2H")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "10M5I5M10I6M14S2H"
  }

  it should "also consume any extra bases if it ends in an insertion when soft-clipping" in {
    val rec = r(10, "38M4I8M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 12 // 4 + 8
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "38M12S"
  }

  it should "preserve an insertion that starts immediately after the soft-clipping" in {
    val rec = r(10, "36M4I10M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "36M4I10S"
  }

  it should "remove deletions that immediately precede the clipping" in {
    val rec = r(10, "40M4D10M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10S"
  }

  it should "preserve deletions that are not near the clipped region" in {
    val rec = r(10, "25M4D25M")
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "25M4D15M10S"
  }

  it should "NOT throw an exception, but do nothing, with an unmapped read" in {
    val rec = r(10, "50M")
    rec.unmapped = true
    clipper(Soft).clipEndOfAlignment(rec, 10) shouldBe 0
    rec.start shouldBe 10
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask" in {
    val rec = r(10, "50M")
    clipper(SoftWithMask).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10S"
    rec.bases.drop(40).forall(_ == 'N'.toByte) shouldBe true
    rec.quals.drop(40).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "mask bases and qualities when the clipping mode is SoftWithMask, including existing soft-clips" in {
    val rec = r(10, "40M10S")
    clipper(SoftWithMask).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "30M20S"
    rec.bases.drop(30).forall(_ == 'N'.toByte) shouldBe true
    rec.quals.drop(30).forall(_ == PhredScore.MinValue) shouldBe true
  }

  it should "hard clip 10 bases" in {
    val rec = r(10, "50M")
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    clipper(Hard).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "40M10H"
    rec.bases shouldBe bases.take(40)
    rec.quals shouldBe quals.take(40)
  }

  it should "hard clip 10 extra bases and convert existing soft-clipping to hard clipping" in {
    val rec = r(10, "40M10S")
    val bases = rec.bases.clone()
    val quals = rec.quals.clone()
    clipper(Hard).clipEndOfAlignment(rec, 10) shouldBe 10
    rec.start shouldBe 10
    rec.cigar.toString shouldBe "30M20H"
    rec.bases shouldBe bases.take(30)
    rec.quals shouldBe quals.take(30)
  }

  it should "unmap a read when more clipping is requested than there is alignment to clip" in {
    val rec = r(10, "40M10S")
    clipper(Soft).clipEndOfAlignment(rec, 50) shouldBe 40 // 40M
    rec.unmapped shouldBe true
    rec.cigar shouldBe Cigar.empty
  }

  for (mode <- ClippingMode.values; auto <- Seq(true, false)) {
    it should s"correctly handle auto-trimming of attribute with auto=$auto and mode=$mode" in {
      val rec = r(10, "20M")
      rec("A1") = "AB" * 10
      rec("A2") = Range.inclusive(1, 20).toArray
      rec("B1") = "A" * 10
      rec("B2") = Range.inclusive(1, 10).toArray

      val cut = mode == ClippingMode.Hard && auto
      clipper(mode=mode, autoClip=auto).clipEndOfAlignment(rec, 5) shouldBe 5
      rec[String]("A1")     shouldBe { if (cut) "ABABABABABABABA" else "AB"*10 }
      rec[Array[Int]]("A2") shouldBe Range.inclusive(1, if (cut) 15 else 20).toArray
      rec[String]("B1")     shouldBe  "A" * 10
      rec[Array[Int]]("B2") shouldBe Range.inclusive(1, 10).toArray
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // One or two test each for all the derived methods.
  //////////////////////////////////////////////////////////////////////////////

  "SamRecordClipper.clipStartOfRead" should "expand existing clipping" in {
    val rec1 = r(10, "50M")
    clipper(Soft).clipStartOfRead(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "10S40M"

    val rec2 = r(10, "5S45M")
    clipper(Soft).clipStartOfRead(rec2, 10) shouldBe 5
    rec2.cigar.toString shouldBe "10S40M"

    val rec3 = r(10, "20S30M")
    clipper(Soft).clipStartOfRead(rec3, 10) shouldBe 0
    rec3.cigar.toString shouldBe "20S30M"

    val rec4 = r(10, "20H30M")
    clipper(Soft).clipStartOfRead(rec4, 10) shouldBe 0
    rec4.cigar.toString shouldBe "20H30M"
  }

  it should "convert soft-clipping to hard clipping or masked clipping" in {
    val hard = r(10, "2H8S40M", attrs=Map("az" -> "345678901234567890123456789012345678901234567890"))
    val mask = r(10, "10S40M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    clipper(Hard, autoClip=true).clipStartOfRead(hard, 5)
    clipper(SoftWithMask, autoClip=true).clipStartOfRead(mask, 5)

    hard.cigar.toString shouldBe "5H5S40M"
    hard.bases.length   shouldBe 45
    hard[String]("az")  shouldBe "678901234567890123456789012345678901234567890"

    mask.cigar.toString shouldBe "10S40M"
    mask.basesString.take(5) shouldBe "NNNNN"
    mask[String]("az") shouldBe "12345678901234567890123456789012345678901234567890"
  }

  "SamRecordClipper.clipEndOfRead" should "expand existing clipping" in {
    val rec1 = r(10, "50M")
    clipper(Soft).clipEndOfRead(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "40M10S"

    val rec2 = r(10, "45M5S")
    clipper(Soft).clipEndOfRead(rec2, 10) shouldBe 5
    rec2.cigar.toString shouldBe "40M10S"

    val rec3 = r(10, "30M20S")
    clipper(Soft).clipEndOfRead(rec3, 10) shouldBe 0
    rec3.cigar.toString shouldBe "30M20S"

    val rec4 = r(10, "30M20H")
    clipper(Soft).clipEndOfRead(rec4, 10) shouldBe 0
    rec4.cigar.toString shouldBe "30M20H"
  }

  it should "convert soft-clipping to hard clipping or masked clipping" in {
    val hard = r(10, "40M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val mask = r(10, "40M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    clipper(Hard, autoClip=true).clipEndOfRead(hard, 5)
    clipper(SoftWithMask, autoClip=true).clipEndOfRead(mask, 5)

    hard.cigar.toString shouldBe "40M5S5H"
    hard.bases.length   shouldBe 45
    hard[String]("az")  shouldBe "123456789012345678901234567890123456789012345"

    mask.cigar.toString shouldBe "40M10S"
    mask.basesString.drop(45) shouldBe "NNNNN"
    mask[String]("az") shouldBe "12345678901234567890123456789012345678901234567890"
  }

  "SamRecordClipper.clip5PrimeEndOfAlignment" should "add more clipping to the 5' end" in {
    val rec1 = r(10, "50M", strand=Plus)
    clipper(Soft).clip5PrimeEndOfAlignment(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Plus)
    clipper(Soft).clip5PrimeEndOfAlignment(rec2, 10) shouldBe 10
    rec2.cigar.toString shouldBe "20S30M"

    val rec3 = r(10, "50M", strand=Minus)
    clipper(Soft).clip5PrimeEndOfAlignment(rec3, 10) shouldBe 10
    rec3.cigar.toString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Minus)
    clipper(Soft).clip5PrimeEndOfAlignment(rec4, 10) shouldBe 10
    rec4.cigar.toString shouldBe "30M20S"
  }

  "SamRecordClipper.clip3PrimeEndOfAlignment" should "add more clipping to the 3' end" in {
    val rec1 = r(10, "50M", strand=Minus)
    clipper(Soft).clip3PrimeEndOfAlignment(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Minus)
    clipper(Soft).clip3PrimeEndOfAlignment(rec2, 10) shouldBe 10
    rec2.cigar.toString shouldBe "20S30M"

    val rec3 = r(10, "50M", strand=Plus)
    clipper(Soft).clip3PrimeEndOfAlignment(rec3, 10) shouldBe 10
    rec3.cigar.toString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Plus)
    clipper(Soft).clip3PrimeEndOfAlignment(rec4, 10) shouldBe 10
    rec4.cigar.toString shouldBe "30M20S"
  }

  "SamRecordClipper.clip5PrimeEndOfRead" should "expand existing clipping at the 5' end" in {
    val rec1 = r(10, "50M", strand=Plus)
    clipper(Soft).clip5PrimeEndOfRead(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Plus)
    clipper(Soft).clip5PrimeEndOfRead(rec2, 10) shouldBe 0
    rec2.cigar.toString shouldBe "10S40M"

    val rec3 = r(10, "50M", strand=Minus)
    clipper(Soft).clip5PrimeEndOfRead(rec3, 10) shouldBe 10
    rec3.cigar.toString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Minus)
    clipper(Soft).clip5PrimeEndOfRead(rec4, 10) shouldBe 0
    rec4.cigar.toString shouldBe "40M10S"
  }

  it should "convert soft-clipping to hard clipping or masked clipping" in {
    val hard = r(10, "10S40M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val mask = r(10, "10S40M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    clipper(Hard, autoClip=true).clip5PrimeEndOfRead(hard, 5)
    clipper(SoftWithMask, autoClip=true).clip5PrimeEndOfRead(mask, 5)

    hard.cigar.toString shouldBe "5H5S40M"
    hard.bases.length   shouldBe 45
    hard[String]("az")  shouldBe "678901234567890123456789012345678901234567890"

    mask.cigar.toString shouldBe "10S40M"
    mask.basesString.take(5) shouldBe "NNNNN"
    mask[String]("az") shouldBe "12345678901234567890123456789012345678901234567890"
  }

  "SamRecordClipper.clip3PrimeEndOfRead" should "expand existing clipping at the 3' end" in {
    val rec1 = r(10, "50M", strand=Minus)
    clipper(Soft).clip3PrimeEndOfRead(rec1, 10) shouldBe 10
    rec1.cigar.toString shouldBe "10S40M"

    val rec2 = r(10, "10S40M", strand=Minus)
    clipper(Soft).clip3PrimeEndOfRead(rec2, 10) shouldBe 0
    rec2.cigar.toString shouldBe "10S40M"

    val rec3 = r(10, "50M", strand=Plus)
    clipper(Soft).clip3PrimeEndOfRead(rec3, 10) shouldBe 10
    rec3.cigar.toString shouldBe "40M10S"

    val rec4 = r(10, "40M10S", strand=Plus)
    clipper(Soft).clip3PrimeEndOfRead(rec4, 10) shouldBe 0
    rec4.cigar.toString shouldBe "40M10S"
  }

  it should "convert soft-clipping to hard clipping or masked clipping" in {
    val hard = r(10, "40M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val mask = r(10, "40M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    clipper(Hard, autoClip=true).clip3PrimeEndOfRead(hard, 5)
    clipper(SoftWithMask, autoClip=true).clip3PrimeEndOfRead(mask, 5)

    hard.cigar.toString shouldBe "40M5S5H"
    hard.bases.length   shouldBe 45
    hard[String]("az")  shouldBe "123456789012345678901234567890123456789012345"

    mask.cigar.toString shouldBe "40M10S"
    mask.basesString.drop(45) shouldBe "NNNNN"
    mask[String]("az") shouldBe "12345678901234567890123456789012345678901234567890"
  }

  "SamRecordClipper.upgradeAllClipping" should "convert leading and trailing soft clips" in {
    val noAuto   = r(10, "5S35M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val withAuto = r(10, "5S35M10S", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))

    clipper(Hard).upgradeAllClipping(noAuto) shouldBe (5, 10)
    noAuto.cigar.toString shouldBe "5H35M10H"
    noAuto.bases.length   shouldBe 35
    noAuto[String]("az")  shouldBe "12345678901234567890123456789012345678901234567890"

    clipper(Hard, autoClip=true).upgradeAllClipping(withAuto) shouldBe (5, 10)
    withAuto.cigar.toString shouldBe "5H35M10H"
    withAuto.bases.length   shouldBe 35
    withAuto[String]("az")  shouldBe "67890123456789012345678901234567890"
  }

  it should "convert soft clips that follow hard clips" in {
    val noAuto   = r(10, "5H5S35M10S5H", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val withAuto = r(10, "5H5S35M10S5H", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))

    clipper(Hard).upgradeAllClipping(noAuto) shouldBe (5, 10)
    noAuto.cigar.toString shouldBe "10H35M15H"
    noAuto.bases.length   shouldBe 35
    noAuto[String]("az")  shouldBe "12345678901234567890123456789012345678901234567890"

    clipper(Hard, autoClip=true).upgradeAllClipping(withAuto) shouldBe (5, 10)
    withAuto.cigar.toString shouldBe "10H35M15H"
    withAuto.bases.length   shouldBe 35
    withAuto[String]("az")  shouldBe "67890123456789012345678901234567890"
  }

  it should "not convert reads that have no soft-clipping" in {
    val noSoft = r(10, "55M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val hard   = r(10, "5H55M10H", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))

    Seq(noSoft, hard).foreach { rec =>
      val cigar = rec.cigar.toString
      clipper(Hard).upgradeAllClipping(rec) shouldBe (0, 0)
      rec.cigar.toString shouldBe cigar
      rec.bases.length   shouldBe 55
      rec[String]("az")  shouldBe "12345678901234567890123456789012345678901234567890"
    }
  }

  it should "not convert reads if they are unmapped or if the clipping mode is not Hard" in {
    val mapped   = r(10, "55M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    val unmapped = r(10, "55M", attrs=Map("az" -> "12345678901234567890123456789012345678901234567890"))
    unmapped.unmapped = true

    clipper(Soft).upgradeAllClipping(mapped) shouldBe (0, 0)
    clipper(SoftWithMask).upgradeAllClipping(mapped) shouldBe (0, 0)
    clipper(Hard).upgradeAllClipping(unmapped) shouldBe (0, 0)
  }
}
