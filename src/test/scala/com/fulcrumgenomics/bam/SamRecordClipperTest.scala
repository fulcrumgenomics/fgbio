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
import org.scalatest.OptionValues

class SamRecordClipperTest extends UnitSpec with OptionValues {
  /** Returns a fragment SAM record with the start / cigar / strand requested. */
  def r(start: Int, cigar: String, strand: Strand = Plus, attrs:Map[String,Any] = Map.empty): SamRecord = {
    new SamBuilder(readLength=Cigar(cigar).lengthOnQuery)
      .addFrag(start=start, cigar=cigar, strand=strand, attrs=attrs).value
  }

  /** Returns a pair of SAM records (read pair) with the respective starts / cigars / strands requested. */
  def pair(start1: Int, cigar1: String, strand1: Strand = Plus, start2: Int, cigar2: String, strand2: Strand = Minus): (SamRecord, SamRecord) = {
    val recs = new SamBuilder(readLength=Cigar(cigar1).lengthOnQuery)
      .addPair(start1=start1, cigar1=cigar1, strand1=strand1, start2=start2, cigar2=cigar2, strand2=strand2)
    recs match { case Seq(r1, r2) => (r1, r2) }
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

  "SamRecordClipper.clipOverlappingReads" should "not clip if the reads are not overlapping" in {
    val (rec, mate) = pair(1, "100M", Plus, 101, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 0)
    rec.start shouldBe 1
    rec.cigar.toString shouldBe "100M"
    mate.start shouldBe 101
    mate.cigar.toString shouldBe "100M"
  }

  it should "clip reads that overlap by one base" in {
    val (rec, mate) = pair(1, "100M", Plus, 100, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 1)
    rec.start shouldBe 1
    rec.cigar.toString shouldBe "100M"
    mate.start shouldBe 101
    mate.cigar.toString shouldBe "1S99M"
  }

  it should "clip reads that overlap by two bases" in {
    val (rec, mate) = pair(2, "100M", Plus, 100, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (1, 1)
    rec.start shouldBe 2
    rec.end   shouldBe 100
    rec.cigar.toString shouldBe "99M1S"
    mate.start shouldBe 101
    mate.end shouldBe 199
    mate.cigar.toString shouldBe "1S99M"
  }

  it should "clip reads that overlap with a deletion" in {
    val (rec, mate) = pair(2, "80M10D10M", Plus, 70, "10M10D80M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (10, 10)
    rec.start shouldBe 2
    rec.end shouldBe 81
    rec.cigar.toString shouldBe "80M10S"
    mate.start shouldBe 90
    mate.end shouldBe 169
    mate.cigar.toString shouldBe "10S80M"
  }

  it should "clip reads that overlap with soft-clipping on the forward read after the midpoint" in {
    val (rec, mate) = pair(1, "95M5S", Plus, 50, "100M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (20, 26)
    rec.start shouldBe 1
    rec.end shouldBe 75
    rec.cigar.toString shouldBe "75M25H"
    mate.start shouldBe 76
    mate.cigar.toString shouldBe "26H74M"
  }

  it should "clip reads that overlap with soft-clipping on the reverse read before the midpoint" in {
    val (rec, mate) = pair(1, "100M", Plus, 55, "5S95M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (25, 21)
    rec.start shouldBe 1
    rec.end shouldBe 75
    rec.cigar.toString shouldBe "75M25H"
    mate.start shouldBe 76
    mate.cigar.toString shouldBe "26H74M"
  }

  it should "clip reads that overlap 1 bp directly in the middle of the pair" in {
    val (rec, mate) = pair(1, "99M1S", Plus, 99, "1S99M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 1)
    rec.start shouldBe 1
    rec.end shouldBe 99
    rec.cigar.toString shouldBe "99M1H"
    mate.start shouldBe 100
    mate.cigar.toString shouldBe "2H98M"
  }

  // Template spans chr1:1-169, with mid-point=85
  // Second read doesn't overlap the mid-point so all clipping should be done to the first read
  it should "clip reads that overlap in the first half of the pair, not over the middle" in {
    val (rec, mate) = pair(1, "95M5S", Plus, 90, "20S80M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (6, 0)
    rec.start shouldBe 1
    rec.end shouldBe 89
    rec.cigar.toString shouldBe "89M11H"
    mate.start shouldBe 90
    mate.cigar.toString shouldBe "20H80M"
  }

  // Template spans chr1:1-164, with mid-point=82
  // First read doesn't overlap the mid-point so all clipping should be done to the first read
  it should "clip reads that overlap in the second half of the pair, not over the middle" in {
    val (rec, mate) = pair(1, "80M20S", Plus, 70, "5S95M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 11)
    rec.start shouldBe 1
    rec.end shouldBe 80
    rec.cigar.toString shouldBe "80M20H"
    mate.start shouldBe 81
    mate.cigar.toString shouldBe "16H84M"
  }

  it should "clip reads when the forward read is much longer" in {
    val (rec, mate) = pair(1, "100M", Plus, 30, "80S20M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (71, 0)
    rec.start shouldBe 1
    rec.end shouldBe 29
    rec.cigar.toString shouldBe "29M71H"
    mate.start shouldBe 30
    mate.cigar.toString shouldBe "80H20M"
  }

  it should "clip reads when the reverse read is much longer" in {
    val (rec, mate) = pair(50, "20M80S", Plus, 1, "100M", Minus)
    clipper(Hard).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 69)
    rec.start shouldBe 50
    rec.end shouldBe 69
    rec.cigar.toString shouldBe "20M80H"
    mate.start shouldBe 70
    mate.cigar.toString shouldBe "69H31M"
  }

  it should "clip reads that overlap with one end having a deletion" in {
    val (rec, mate) = pair(1, "60M10D40M", Plus, 50, "10M10D80M10D10M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (25, 26)
    rec.start shouldBe 1
    rec.end shouldBe 85
    rec.cigar.toString shouldBe "60M10D15M25S"
    mate.start shouldBe 86
    mate.cigar.toString shouldBe "26S64M10D10M"
  }

  it should "clip reads that overlap with one end having a deletion with mismatching cigars" in {
    val (rec, mate) = pair(1, "100M", Plus, 50, "10M10D80M10D10M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (15, 26)
    rec.start shouldBe 1
    rec.end shouldBe 85
    rec.cigar.toString shouldBe "85M15S"
    mate.start shouldBe 86
    mate.cigar.toString shouldBe "26S64M10D10M"
  }

  it should "clip reads that fully overlap with both ends having deletions" in {
    // NB: mid point occurs right in the middle of the deletion
    val (rec, mate) = pair(1, "50M10D50M", Plus, 1, "50M10D50M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (50, 50)
    rec.start shouldBe 1
    rec.end shouldBe 50
    rec.cigar.toString shouldBe "50M50S"
    mate.start shouldBe 61
    mate.cigar.toString shouldBe "50S50M"
  }

  it should "clip reads that overlap with both ends having deletions" in {
    // NB: mid point occurs right in the middle of the deletion.  Clipping loses the deletion!
    val (rec, mate) = pair(1, "50M10D50M", Plus, 3, "47M10D53M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (50, 47)
    rec.start shouldBe 1
    rec.end shouldBe 50
    rec.cigar.toString shouldBe "50M50S"
    mate.start shouldBe 60
    mate.cigar.toString shouldBe "47S53M"
  }

  it should "clip reads that fully overlap" in {
    val (rec, mate) = pair(1, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (50, 50)
    rec.start shouldBe 1
    rec.cigar.toString shouldBe "50M50S"
    mate.start shouldBe 51
    mate.cigar.toString shouldBe "50S50M"
  }

  it should "clip reads that extend past each other" in {
    // R1 is 50-149, while R2 is 1-100, so the reference midpoint is (50 + 100) / 2 = 75
    val (rec, mate) = pair(50, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (74, 75)
    rec.start shouldBe 50
    rec.end shouldBe 75
    rec.cigar.toString shouldBe "26M74S"
    mate.start shouldBe 76
    mate.end shouldBe 100
    mate.cigar.toString shouldBe "75S25M"
  }

  it should "clip reads that extend past each other with one read having deletions" in {
    // R1 is 50-149, while R2 is 1-120, so the reference midpoint is (50 + 120) / 2 = 85
    val (rec, mate) = pair(50, "100M", Plus, 1, "10M10D80M10D10M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (64, 75)
    rec.start shouldBe 50
    rec.end shouldBe 85
    rec.cigar.toString shouldBe "36M64S"
    mate.start shouldBe 86
    mate.end shouldBe 120
    mate.cigar.toString shouldBe "75S15M10D10M"
  }

  it should "clip reads that extend past each other with both read having deletions" in {
    // R1 is 50-149, while R2 is 1-120, so the reference midpoint is (50 + 120) / 2 = 85
    val (rec, mate) = pair(50, "50M10D50M", Plus, 1, "10M10D80M10D10M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (64, 75)
    rec.start shouldBe 50
    rec.end shouldBe 85
    rec.cigar.toString shouldBe "36M64S"
    mate.start shouldBe 86
    mate.end shouldBe 120
    mate.cigar.toString shouldBe "75S15M10D10M"
  }

  it should "not clip when the read pairs are mapped +/- with start(R1) > end(R2) but do not overlap" in {
    val (rec, mate) = pair(1000, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipOverlappingReads(rec=rec, mate=mate) shouldBe (0, 0)
    rec.start shouldBe 1000
    rec.cigar.toString shouldBe "100M"
    mate.start shouldBe 1
    mate.cigar.toString shouldBe "100M"
  }

  "SamRecordClipper.clipExtendingPastMateEnds" should "not clip reads that do not extend past each other" in {
    Seq((Plus, Minus), (Minus, Plus)).foreach { case (strand1, strand2) =>
      val (rec, mate) = pair(1, "100M", strand1, 1, "100M", strand2)
      clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0)
      rec.start shouldBe 1
      rec.cigar.toString shouldBe "100M"
      mate.start shouldBe 1
      mate.cigar.toString shouldBe "100M"
    }
  }

  it should "clip reads that extend one base past their mate's start" in {
    val (rec, mate) = pair(2, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (1, 1)
    rec.start shouldBe 2
    rec.cigar.toString shouldBe "99M1S"
    mate.start shouldBe 2
    mate.cigar.toString shouldBe "1S99M"
  }

  it should "clip reads that extend two bases past their mate's start" in {
    val (rec, mate) = pair(3, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (2, 2)
    rec.start shouldBe 3
    rec.cigar.toString shouldBe "98M2S"
    mate.start shouldBe 3
    mate.cigar.toString shouldBe "2S98M"
  }

  it should "clip reads that where both ends extends their mate's start" in {
    val (rec, mate) = pair(51, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (50, 50)
    rec.start shouldBe 51
    rec.cigar.toString shouldBe "50M50S"
    mate.start shouldBe 51
    mate.cigar.toString shouldBe "50S50M"
  }

  it should "clip reads that where only one end extends their mate's start" in {
    val (rec, mate) = pair(1, "100M", Plus, 1, "50S50M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (50, 0)
    rec.start shouldBe 1
    rec.cigar.toString shouldBe "50M50S" // added clipping
    mate.start shouldBe 1
    mate.cigar.toString shouldBe "50S50M" // clipping remains
  }

  it should "clip reads where only one end extends their mate's start that has insertions" in {
    val (rec, mate) = pair(1, "40M10I50M", Plus, 1, "50S50M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (40, 0)
    rec.start shouldBe 1
    rec.cigar.toString shouldBe "40M10I10M40S" // added clipping
    mate.start shouldBe 1
    mate.cigar.toString shouldBe "50S50M" // clipping remains
  }

  it should "clip the forward read when it ends before the mate's start but soft-clipped bases extend past" in {
    val (rec, mate) = pair(20, "30M20S", Plus, 20, "10S40M", Minus)
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0) // No aligned bases clipped
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "30M10S10H" // 10 of 20 bases extending past are hard-clipped, 10 remain soft-clipped
    mate.start shouldBe 20
    mate.cigar.toString shouldBe "10H40M" // hard-clip the 10S that extend past
  }

  it should "clip the forward read when it ends before the mate's start but soft-clipped bases extend past while there is a deletion" in {
    val (rec, mate) = pair(20, "15M1D15M20S", Plus, 20, "10S15M1D25M", Minus)
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0) // No aligned bases clipped
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "15M1D15M10S10H" // 10 of 20 bases extending past are hard-clipped, 10 remain soft-clipped
    mate.start shouldBe 20
    mate.cigar.toString shouldBe "10H15M1D25M" // hard-clip the 10S that extend past
  }

  it should "clip the forward read when it ends before the mate's start but soft-clipped bases extend past while there is an insertion" in {
    val (rec, mate) = pair(20, "15M1I15M20S", Plus, 20, "10S15M1I25M", Minus)
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0) // No aligned bases clipped
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "15M1I15M10S10H" // 10 of 20 bases extending past are hard-clipped, 10 remain soft-clipped
    mate.start shouldBe 20
    mate.cigar.toString shouldBe "10H15M1I25M" // hard-clip the 10S that extend past
  }

  it should "clip the reverse read when it ends before the mate's start but soft-clipped bases extend past" in {
    val (rec, mate) = pair(20, "40M10S", Plus, 30, "20S30M", Minus)
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0) // No aligned bases clipped
    rec.start shouldBe 20
    rec.cigar.toString shouldBe "40M10H" // hard-clip the 10S that extend past
    mate.start shouldBe 30
    mate.cigar.toString shouldBe "10H10S30M" // hard-clip the 10S that extend past
  }

  it should "not clip when the read pairs are mapped +/- with start(R1) > end(R2) but do not overlap" in {
    val (rec, mate) = pair(1000, "100M", Plus, 1, "100M", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0)
    rec.start shouldBe 1000
    rec.cigar.toString shouldBe "100M"
    mate.start shouldBe 1
    mate.cigar.toString shouldBe "100M"
  }

  it should "not clip when the reads do not extend past each other with insertions" in {
    val (rec, mate) = pair(start1=1, start2=1, strand1=Plus, strand2=Minus, cigar1="40M20I40M", cigar2="40M20I40M")
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (0, 0)
    rec.start shouldBe 1
    rec.cigar.toString() shouldBe "40M20I40M"
    mate.start shouldBe 1
    mate.cigar.toString shouldBe "40M20I40M"
  }

  // The last reference position shared by the two reads is 149.  The record has 50 query bases after it (aligned) and
  // the mate has 10 (soft-clipped), so 40 bases are removed from the record.  Symmetrically, the first shared
  // reference position is 100, before which the record has no query bases and the mate has 40, so 40 bases are removed
  // from the mate.  There are no indels here, so this is the control case that must be unaffected by the fix for
  // https://github.com/fulcrumgenomics/fgbio/issues/1090.
  it should "clip back to the mate's soft-clipped end when neither read contains an indel" in {
    val (rec, mate) = pair(100, "100M", Plus, 60, "90M10S", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (40, 40)
    rec.start shouldBe 100
    rec.cigar.toString shouldBe "60M40S"
    mate.start shouldBe 100
    mate.cigar.toString shouldBe "40S50M10S"
  }

  // Regression test for https://github.com/fulcrumgenomics/fgbio/issues/1090.  The mate's un-soft-clipped end
  // (223 + 2 = 225) is a reference position built by adding a query distance (the mate's trailing soft-clipping) to a
  // reference position, and it happens to land inside the record's deletion.  The reads share reference positions
  // through 223; after it the record has four query bases and the mate has two, so only two bases may be clipped.
  it should "clip using query distances when a deletion falls at the mate's un-soft-clipped end" in {
    val (rec, mate) = pair(101, "2S124M1D3M", Plus, 100, "3S124M2S", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (2, 0)
    rec.mapped shouldBe true
    rec.start shouldBe 101
    rec.cigar.toString shouldBe "2S124M1D1M2S"
    mate.mapped shouldBe true
    mate.start shouldBe 100
    mate.cigar.toString shouldBe "3S124M2S"
  }

  // A second regression test for https://github.com/fulcrumgenomics/fgbio/issues/1090 that pins the knock-on effect on
  // the mate: over-clipping the record un-maps it, and the mate is then compared against an un-mapped record and so is
  // left unclipped.  Both reads must be clipped here, and the amount clipped from each must not depend on the other
  // read having been clipped first.
  it should "clip the mate even when the record's deletion falls at the mate's un-soft-clipped end" in {
    val (rec, mate) = pair(101, "2S124M1D3M", Plus, 97, "115M14S", Minus)
    clipper(Soft).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (2, 2)
    rec.mapped shouldBe true
    rec.start shouldBe 101
    rec.cigar.toString shouldBe "2S124M1D1M2S"
    mate.mapped shouldBe true
    mate.start shouldBe 99
    mate.cigar.toString shouldBe "2S113M14S"
  }

  // The example from https://github.com/fulcrumgenomics/fgbio/issues/1090.  The reads share reference positions
  // through 169; after it the record has 80 query bases (10 inserted, 23 aligned and 47 soft-clipped) and the mate has
  // 30 (all soft-clipped), so 50 bases are removed from the record.  Comparing reference distances instead under-clips
  // here, since the insertion makes the record's alignment end 10 reference bases short of where its query bases reach.
  it should "clip using query distances when the record contains an insertion before the mate's end" in {
    val (rec, mate) = pair(100, "70M10I23M47S", Plus, 100, "50S70M30S", Minus)
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (3, 0) // 3 aligned bases clipped from the record
    rec.start shouldBe 100
    rec.cigar.toString shouldBe "70M10I20M50H"
    mate.start shouldBe 100
    mate.cigar.toString shouldBe "50H70M30S"
  }

  // Bases already hard-clipped from the 3' end are not present in the record, so they must not be counted towards the
  // clipping that is applied, otherwise the requested clipping is silently absorbed by the existing hard-clipping.
  it should "clip past the mate's end when the read has already been hard-clipped at its 3' end" in {
    val (rec, mate) = pair(100, "100M", Plus, 60, "90M10S", Minus)
    clipper(Hard).clip3PrimeEndOfRead(rec, 10) shouldBe 10
    clipper(Hard).clip3PrimeEndOfRead(mate, 10) shouldBe 10
    rec.cigar.toString shouldBe "90M10H"
    mate.cigar.toString shouldBe "10H80M10S"
    clipper(Hard).clipExtendingPastMateEnds(rec=rec, mate=mate) shouldBe (30, 30)
    rec.start shouldBe 100
    rec.cigar.toString shouldBe "60M40H"
    mate.start shouldBe 100
    mate.cigar.toString shouldBe "40H50M10S"
  }

  "SamRecordClipper.numBasesExtendingPastMate" should "return zero when reads do not extend past the end or are not FR pairs" in {
    val builder = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=100, start2=200, cigar1="20S80M", cigar2="10S90M")
    val Seq(r3, r4) = builder.addPair(start1=100, start2=200, cigar1="10S90M", cigar2="20S80M")
    val Seq(r5, r6) = builder.addPair(start1=100, start2=100, cigar1="10S90M", cigar2="20S80M", strand1=Plus, strand2=Plus)
    val Seq(r7, r8) = builder.addPair(start1=100, start2=100, cigar1="10S90M", cigar2="20S80M", strand1=Minus, strand2=Minus)

    Seq(r1, r2, r3, r4, r5, r6, r7, r8).foreach { r =>
      numBasesExtendingPastMate(clipper(Soft), r) shouldBe 0
    }
  }

  /** Convenience method for testing numBasesExtendingPastMate */
  def numBasesExtendingPastMate(clipper: SamRecordClipper, rec: SamRecord): Int = {
    val mateCigar = rec.mateCigar.getOrElse(throw new IllegalStateException(f"Mate cigar (MC SAM tag) needed for read: ${rec.name}"))
    clipper.numBasesExtendingPastMate(rec=rec, mateStart=rec.mateStart, mateCigar=mateCigar)
  }

  it should "return a return a positive value when reads extend past its mate" in {
    val builder = new SamBuilder()

    // both extend 50bp beyond the mates
    val Seq(r1, r2) = builder.addPair(start1=100, start2=50, cigar1="100M", cigar2="100M")
    numBasesExtendingPastMate(clipper(Soft), r1) shouldBe 50
    numBasesExtendingPastMate(clipper(Soft), r2) shouldBe 50

    // r3 ends at 149, which is the same as the end of r4 at 149
    // r4 starts at 51, which is the same as the start of r3 at 51
    val Seq(r3, r4) = builder.addPair(start1=100, start2=100, cigar1="50S50M", cigar2="50S50M")
    numBasesExtendingPastMate(clipper(Soft), r3) shouldBe 0
    numBasesExtendingPastMate(clipper(Soft), r4) shouldBe 0

    // r5 ends at 199, which is the same as the end of r6 at 199
    // r6 starts at 51, which is the same as the start of 45 at 51
    val Seq(r5, r6) = builder.addPair(start1=100, start2=100, cigar1="50M50S", cigar2="50M50S")
    numBasesExtendingPastMate(clipper(Soft), r5) shouldBe 0
    numBasesExtendingPastMate(clipper(Soft), r6) shouldBe 0

    // r7 ends at 149, which is before the end of r8 at 169
    // r8 starts at 71, which is after the start of r7 at 51
    val Seq(r7, r8) = builder.addPair(start1=100, start2=100, cigar1="50S50M", cigar2="30S70M")
    numBasesExtendingPastMate(clipper(Soft), r7) shouldBe 0
    numBasesExtendingPastMate(clipper(Soft), r8) shouldBe 0

    // r9 ends at 169, which is after the end of r10 at 149
    // r10 starts at 51, which is before the start of r9 at 71
    val Seq(r9, r10) = builder.addPair(start1=100, start2=100, cigar1="30S70M", cigar2="50S50M")
    numBasesExtendingPastMate(clipper(Soft), r9) shouldBe 20
    numBasesExtendingPastMate(clipper(Soft), r10) shouldBe 20
  }

  it should "count query bases, not reference bases, when an indel is present near the end of an alignment" in {
    val builder = new SamBuilder(readLength=129)

    // The reads share reference positions through 223, after which r1 has four query bases (the last base of its 124M
    // at 224 and the three bases of its trailing 3M at 226-228; the deleted position 225 contributes no query base)
    // and r2 has two (its trailing 2S).  The first shared reference position is 101, before which r1 has two query
    // bases (its leading 2S) and r2 has four (its leading 3S plus the base at 100).
    val Seq(r1, r2) = builder.addPair(start1=101, start2=100, cigar1="2S124M1D3M", cigar2="3S124M2S")
    numBasesExtendingPastMate(clipper(Soft), r1) shouldBe 2
    numBasesExtendingPastMate(clipper(Soft), r2) shouldBe 2

    // The reads share reference positions through 169, after which r3 has 80 query bases (10 inserted, 23 aligned and
    // 47 soft-clipped) and r4 has 30 (all soft-clipped).  The first shared reference position is 100, before which r3
    // has no query bases and r4 has 50.
    val builder2    = new SamBuilder(readLength=150)
    val Seq(r3, r4) = builder2.addPair(start1=100, start2=100, cigar1="70M10I23M47S", cigar2="50S70M30S")
    numBasesExtendingPastMate(clipper(Soft), r3) shouldBe 50
    numBasesExtendingPastMate(clipper(Soft), r4) shouldBe 50
  }

  // Both alignments lie within the insert, so a pair that reads through into the adapter can still leave them
  // disjoint: here each 100 base read aligns only its first 20 bases and soft-clips the rest, over an insert of 39
  // reference bases.  With the mate at 1019 the two alignments share exactly one reference position, 1019, past which
  // the record has all 80 of its soft-clipped bases and the mate has 19, so 61 bases of read-through are clipped from
  // each.  Moving the mate one base to the right removes that shared position, and the answer must not jump: the
  // insert is one base longer, so one fewer base of each read is read-through.  Comparing the record's soft-clipped
  // extent against the mate's un-soft-clipped end is what carries the count across that boundary, and no indel can lie
  // between two alignments that place no bases between them.
  it should "clip the read-through of a pair whose alignments share no reference position" in {
    val Seq(r1, r2) = new SamBuilder(readLength=100).addPair(start1=1000, start2=1019, cigar1="20M80S", cigar2="80S20M")
    numBasesExtendingPastMate(clipper(Soft), r1) shouldBe 61
    numBasesExtendingPastMate(clipper(Soft), r2) shouldBe 61

    val Seq(r3, r4) = new SamBuilder(readLength=100).addPair(start1=1000, start2=1020, cigar1="20M80S", cigar2="80S20M")
    numBasesExtendingPastMate(clipper(Soft), r3) shouldBe 60
    numBasesExtendingPastMate(clipper(Soft), r4) shouldBe 60

    // The count is bounded by the record's own soft clipping, so applying it takes no aligned base: hard-clipping
    // shows exactly which bases go, where soft-clipping them again would be invisible.
    clipper(Hard).clipExtendingPastMateEnds(rec=r3, mate=r4) shouldBe (0, 0) // No aligned bases clipped
    r3.start shouldBe 1000
    r3.cigar.toString shouldBe "20M20S60H"
    r4.start shouldBe 1020
    r4.cigar.toString shouldBe "60H20S20M"
  }

  // The mate's alignment is supplied by the caller — from the MC tag in practice — and is never cross-checked against
  // what the record's own flags say about the pair, so a record can be classified as one end of an FR pair while the
  // mate alignment in hand lies entirely to its left.  Its 3' end then points away from that mate, which describes an
  // outward-facing template rather than read-through, and nothing may be clipped: there is no estimate to make, and
  // the record's aligned bases are not the ones to spend on a geometry this function cannot verify.
  it should "not clip a record lying entirely past the mate alignment it is given" in {
    val (rec, _) = pair(301, "15M", Plus, 300, "10M5S", Minus)
    rec.isFrPair shouldBe true
    clipper(Soft).numBasesExtendingPastMate(rec=rec, mateStart=101, mateCigar=Cigar("10M5S")) shouldBe 0
  }
}
