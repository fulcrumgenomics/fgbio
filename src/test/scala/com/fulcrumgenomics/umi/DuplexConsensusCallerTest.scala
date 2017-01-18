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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.testing.SamRecordSetBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore

class DuplexConsensusCallerTest extends UnitSpec {
  // Function to create a caller without so much writing
  def caller(q: Int = 10, pre: Int = DuplexConsensusCaller.ErrorRatePreUmi, post: Int = DuplexConsensusCaller.ErrorRatePostUmi) =
    new DuplexConsensusCaller(readNamePrefix="test", minInputBaseQuality=q.toByte, errorRatePreUmi=pre.toByte, errorRatePostUmi=post.toByte)

  // A default caller that can be used by tests that don't care about the parameters
  private val c = caller()

  private val MI = ConsensusTags.MolecularId

  "DuplexConsensusCaller.sourceMoleculeId" should "strip of that last /suffix" in {
    val builder = new SamRecordSetBuilder()
    builder.addFrag(start=1, attrs=Map(MI -> "foo/A")).foreach { r => c.sourceMoleculeId(r) shouldBe "foo" }
    builder.addFrag(start=1, attrs=Map(MI -> "f/A/B")).foreach { r => c.sourceMoleculeId(r) shouldBe "f/A" }
    builder.addFrag(start=1, attrs=Map(MI -> "f/abc")).foreach { r => c.sourceMoleculeId(r) shouldBe "f"   }
    builder.addFrag(start=1, attrs=Map(MI -> "f///A")).foreach { r => c.sourceMoleculeId(r) shouldBe "f//" }
  }

  it should "throw an exception if the MI tag doesn't have a /suffix" in {
    an[Exception] shouldBe thrownBy {
      val builder = new SamRecordSetBuilder()
      builder.addFrag(start=1, attrs=Map(MI -> "foo")).foreach { r => c.sourceMoleculeId(r) }
    }
  }

  it should "throw an exception if the MI tag is null" in {
    an[Exception] shouldBe thrownBy {
      val builder = new SamRecordSetBuilder()
      builder.addFrag(start=1).foreach { r => c.sourceMoleculeId(r) }
    }
  }

  "DuplexConsensusCaller.consensusReadsFromSamRecords" should "not create records from fragments" in {
    val builder = new SamRecordSetBuilder(readLength=10)
    builder.addFrag(start=1, attrs=Map(MI -> "foo/A")).foreach(r => r.setReadString("AAAAAAAAAA"))
    builder.addFrag(start=1, attrs=Map(MI -> "foo/B")).foreach(r => r.setReadString("AAAAAAAAAA"))
    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs shouldBe empty
  }

  it should "create a simple double stranded consensus for a pair of A and a pair of B reads" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
    }

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("No first of pair."))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("No second of pair."))
    r1.getReadString shouldBe "AAAAAAAAAA"
    r2.getReadString shouldBe "GGGGGGGGGG" // after un-revcomping
    r1.getBaseQualities()(0) > 30 shouldBe true
    r2.getBaseQualities()(0) > 30 shouldBe true
  }

  it should "not create a a consensus with only A or B reads" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, attrs=Map(MI -> "foo/A")).foreach { r => r.setReadString("AAAAAAAAAA") }
    builder.addPair(name="q2", start1=100, start2=200, attrs=Map(MI -> "foo/A")).foreach { r => r.setReadString("AAAAAAAAAA") }
    builder.addPair(name="q3", start1=100, start2=200, attrs=Map(MI -> "foo/A")).foreach { r => r.setReadString("AAAAAAAAAA") }

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs shouldBe empty
  }

  it should "generate a sensible consensus from deep data" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    Range.inclusive(1, 50) foreach { i =>
      builder.addPair(name=s"a$i", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
        r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
      }
      builder.addPair(name=s"b$i", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
        r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
      }
    }

    val recs = caller(post=45.toByte).consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("No first of pair."))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("No second of pair."))
    r1.getReadString shouldBe "AAAAAAAAAA"
    r2.getReadString shouldBe "GGGGGGGGGG" // after un-revcomping

    Seq(r1, r2) foreach { r =>
      r.getBaseQualities.forall(_ == 90) shouldBe true
      r.getIntegerAttribute(ConsensusTags.PerRead.RawReadCount) shouldBe 100
      r.getIntegerAttribute(ConsensusTags.PerRead.AbRawReadCount) shouldBe 50
      r.getIntegerAttribute(ConsensusTags.PerRead.BaRawReadCount) shouldBe 50
      r.getFloatAttribute(ConsensusTags.PerRead.RawReadErrorRate) shouldBe 0.0f
    }
  }

  it should "not saturate the qualities with deep AB and light BA coverage" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    Range.inclusive(1, 50) foreach { i =>
      builder.addPair(name=s"a$i", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
        r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
      }
    }
    // Single BA read
    builder.addPair(name=s"b1", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
    }

    val recs = caller(post=45.toByte).consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("No first of pair."))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("No second of pair."))
    r1.getReadString shouldBe "AAAAAAAAAA"
    r2.getReadString shouldBe "GGGGGGGGGG" // after un-revcomping

    Seq(r1, r2) foreach { r =>
      r.getBaseQualities.forall(_ <= (45+20)) shouldBe true
    }
  }

  it should "generate a sensible consensus when AB and BAs are sequenced/clipped to different lengths" in {
    val ab = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    val ba = new SamRecordSetBuilder(readLength= 8, baseQuality=20)
    ab.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    ba.addPair(name="q2", start1=202, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCC" else "AAAAAAAA")
    }

    val recs = c.consensusReadsFromSamRecords(ab.toSeq ++ ba.toSeq)
    recs should have size 2
    recs.foreach { r => r.getReadLength shouldBe 8 }
  }

  it should "respect minimum quality when building consensuses" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
      if (r.getSecondOfPairFlag) r.getBaseQualities()(5) = 5
    }

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("Missing R1"))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("Missing R2"))
    r1.getReadString shouldBe "AAAAANAAAA"
    r1.getBaseQualities()(5) shouldBe PhredScore.MinValue
    r2.getReadString shouldBe "GGGGGGGGGG"
  }

  it should "emit a no-call when there are conflicting SS consensus bases with equal quality" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "TAAAAAAAAA")
    }

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("Missing R1"))
    r1.getReadString shouldBe "NAAAAAAAAA"
    r1.getBaseQualities()(0) shouldBe PhredScore.MinValue
  }

  it should "emit a low qual call when there are conflicting SS consensus bases with unequal quality" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=30)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "TAAAAAAAAA")
      if (r.getSecondOfPairFlag) r.getBaseQualities()(0) = 15
    }

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("Missing R1"))
    r1.getReadString shouldBe "AAAAAAAAAA"
    r1.getBaseQualities()(0).toInt shouldBe 15 +- 5
  }

  it should "generate correct summary and detail tags on consensus reads" in {
    val builder = new SamRecordSetBuilder(readLength=10, baseQuality=30)
    // 3 AB pairs:
    //   First pair is good
    //   Second pair has R1 with a no-call at base index 2, and a low quality base a index 3
    //   Third pair  has one error at index 4 in each read
    builder.addPair(name="ab1", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAAAAAAA" else "CCCCCCCCCC")
    }
    builder.addPair(name="ab2", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      if (r.getSecondOfPairFlag) r.setReadString("CCCCCCCCCC")
      else {
        r.setReadString("AANAAAAAAA")
        r.getBaseQualities()(2) = PhredScore.MinValue
        r.getBaseQualities()(3) = 10.toByte
      }
    }
    builder.addPair(name="ab3", start1=100, start2=200, strand1=Plus, strand2=Minus, attrs=Map(MI -> "foo/A")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "AAAAGAAAAA" else "CCCCCTCCCC")
    }

    // 2 BA pairs (recall that first of pair is R2 consensus and vice versa here):
    //    First pair is good
    //    Second pair is good
    builder.addPair(name="ba1", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
    }
    builder.addPair(name="ba2", start1=200, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "foo/B")).foreach { r =>
      r.setReadString(if (r.getFirstOfPairFlag) "CCCCCCCCCC" else "AAAAAAAAAA")
    }

    val recs = caller(q=20).consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.getFirstOfPairFlag).getOrElse(fail("R1 missing."))
    val r2 = recs.find(_.getSecondOfPairFlag).getOrElse(fail("R2 missing."))

    { // Check the per-read tags
      import ConsensusTags.PerRead._
      r1.getIntegerAttribute(RawReadCount)      shouldBe 5
      r1.getIntegerAttribute(AbRawReadCount)    shouldBe 3
      r1.getIntegerAttribute(BaRawReadCount)    shouldBe 2
      r1.getIntegerAttribute(MinRawReadCount)   shouldBe 4 // where the no-call or low-qual base are
      r1.getIntegerAttribute(AbMinRawReadCount) shouldBe 2 // where the no-call or low-qual base are
      r1.getIntegerAttribute(BaMinRawReadCount) shouldBe 2
      r1.getFloatAttribute(RawReadErrorRate)    shouldBe (1 / 48f) // 50 bases, minus 2 no/low-qual bases = 48
      r1.getFloatAttribute(AbRawReadErrorRate)  shouldBe (1 / 28f)
      r1.getFloatAttribute(BaRawReadErrorRate)  shouldBe 0f

      r2.getIntegerAttribute(RawReadCount)      shouldBe 5
      r2.getIntegerAttribute(AbRawReadCount)    shouldBe 3
      r2.getIntegerAttribute(BaRawReadCount)    shouldBe 2
      r2.getIntegerAttribute(MinRawReadCount)   shouldBe 5
      r2.getIntegerAttribute(AbMinRawReadCount) shouldBe 3
      r2.getIntegerAttribute(BaMinRawReadCount) shouldBe 2
      r2.getFloatAttribute(RawReadErrorRate)    shouldBe (1 / 50f)
      r2.getFloatAttribute(AbRawReadErrorRate)  shouldBe (1 / 30f)
      r2.getFloatAttribute(BaRawReadErrorRate)  shouldBe 0f
    }

    { // Check the per-base tags
      import ConsensusTags.PerBase._
      r1.getSignedShortArrayAttribute(AbRawReadCount)  shouldBe Array[Byte](3,3,2,2,3,3,3,3,3,3)
      r1.getSignedShortArrayAttribute(AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
      r1.getSignedShortArrayAttribute(BaRawReadCount)  shouldBe Array[Byte](2,2,2,2,2,2,2,2,2,2)
      r1.getSignedShortArrayAttribute(BaRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,0,0,0,0,0)
      r2.getSignedShortArrayAttribute(AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
      r2.getSignedShortArrayAttribute(AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
      r2.getSignedShortArrayAttribute(BaRawReadCount)  shouldBe Array[Byte](2,2,2,2,2,2,2,2,2,2)
      r2.getSignedShortArrayAttribute(BaRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,0,0,0,0,0)
    }
  }
}
