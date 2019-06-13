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

import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.ConsensusTags.PerRead.{AbRawReadErrorRate, BaRawReadErrorRate, RawReadErrorRate}
import com.fulcrumgenomics.util.NumericTypes.PhredScore

class DuplexConsensusCallerTest extends UnitSpec {
  // Function to create a caller without so much writing
  def caller(q: Int = 10, pre: Int = DuplexConsensusCaller.ErrorRatePreUmi, post: Int = DuplexConsensusCaller.ErrorRatePostUmi, minReads: Seq[Int] = Seq(1)) =
    new DuplexConsensusCaller(readNamePrefix="test", minInputBaseQuality=q.toByte, errorRatePreUmi=pre.toByte, errorRatePostUmi=post.toByte,
      minReads=minReads)

  // A default caller that can be used by tests that don't care about the parameters
  private val c = caller()

  private val MI = ConsensusTags.MolecularId
  private val RX = ConsensusTags.UmiBases


  "DuplexConsensusCaller.sourceMoleculeId" should "strip of that last /suffix" in {
    val builder = new SamBuilder()
    builder.addFrag(start=1, attrs=Map(MI -> "foo/A")).foreach { r => c.sourceMoleculeId(r) shouldBe "foo" }
    builder.addFrag(start=1, attrs=Map(MI -> "f/A/B")).foreach { r => c.sourceMoleculeId(r) shouldBe "f/A" }
    builder.addFrag(start=1, attrs=Map(MI -> "f/abc")).foreach { r => c.sourceMoleculeId(r) shouldBe "f"   }
    builder.addFrag(start=1, attrs=Map(MI -> "f///A")).foreach { r => c.sourceMoleculeId(r) shouldBe "f//" }
  }

  it should "throw an exception if the MI tag doesn't have a /suffix" in {
    an[Exception] shouldBe thrownBy {
      val builder = new SamBuilder()
      builder.addFrag(start=1, attrs=Map(MI -> "foo")).foreach { r => c.sourceMoleculeId(r) }
    }
  }

  it should "throw an exception if the MI tag is null" in {
    an[Exception] shouldBe thrownBy {
      val builder = new SamBuilder()
      builder.addFrag(start=1).foreach { r => c.sourceMoleculeId(r) }
    }
  }

  "DuplexConsensusCaller.consensusReadsFromSamRecords" should "not create records from fragments" in {
    val builder = new SamBuilder(readLength=10)
    builder.addFrag(start=1, attrs=Map(MI -> "foo/A"), bases="AAAAAAAAAA")
    builder.addFrag(start=1, attrs=Map(MI -> "foo/B"), bases="AAAAAAAAAA")
    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs shouldBe empty
  }

  it should "create a simple double stranded consensus for a pair of A and a pair of B reads" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("No first of pair."))
    val r2 = recs.find(_.secondOfPair).getOrElse(fail("No second of pair."))
    r1.basesString shouldBe "AAAAAAAAAA"
    r2.basesString shouldBe "GGGGGGGGGG" // after un-revcomping
    r1.quals(0) > 30 shouldBe true
    r2.quals(0) > 30 shouldBe true
  }

  Seq(true).foreach { allowSingleStrand =>
    val extraName: String = if (allowSingleStrand) "" else "not "
    Seq("A", "B").foreach { duplexStrand =>
      it should s"$extraName create a consensus with only $duplexStrand reads" in {
        val builder = new SamBuilder(readLength=10, baseQuality=20)
        builder.addPair(name="q1", start1=100, start2=200, attrs=Map(MI -> s"foo/$duplexStrand"), bases1="AAAAAAAAAA", bases2="CCCCCCCCCC")
        builder.addPair(name="q2", start1=100, start2=200, attrs=Map(MI -> s"foo/$duplexStrand"), bases1="AAAAAAAAAA", bases2="CCCCCCCCCC")
        builder.addPair(name="q3", start1=100, start2=200, attrs=Map(MI -> s"foo/$duplexStrand"), bases1="AAAAAAAAAA", bases2="CCCCCCCCCC")

        val minReads = if (allowSingleStrand) Seq(1, 1, 0) else Seq(1, 1, 1)
        val recs = caller(minReads=minReads).consensusReadsFromSamRecords(builder.toSeq)

        if (allowSingleStrand) {
          recs should have size 2
          val r1 = recs.find(_.firstOfPair).getOrElse(fail("No first of pair."))
          val r2 = recs.find(_.secondOfPair).getOrElse(fail("No second of pair."))
          r1.basesString shouldBe "AAAAAAAAAA"
          r2.basesString shouldBe "GGGGGGGGGG" // after un-revcomping
          r1.quals(0) > 30 shouldBe true
          r2.quals(0) > 30 shouldBe true

          { // Check a few per-read tags
            import ConsensusTags.PerRead._
            r1[Int](AbRawReadCount)        shouldBe 3
            r1[Int](BaRawReadCount)        shouldBe 0
            r2[Int](AbRawReadCount)        shouldBe 3
            r2[Int](BaRawReadCount)        shouldBe 0          }

          { // Check a few per-base tags
            import ConsensusTags.PerBase._
            r1[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
            r1.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
            r2[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
            r2.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
          }
        }
        else {
          recs shouldBe empty
        }
      }
    }
  }

  it should "throw validation exceptions if you provide --min-reads in increasing stringency order" in {
    an[Exception] shouldBe thrownBy { new DuplexConsensusCaller(readNamePrefix="pre", minReads=Seq(1,2,3)) }
    an[Exception] shouldBe thrownBy { new DuplexConsensusCaller(readNamePrefix="pre", minReads=Seq(9,4,6)) }
  }

  /** A builder where we have one duplex, with five total reads, three on one strand and two on the other. */
  private val builderForMinReads: SamBuilder = {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q3", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q4", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))
    builder.addPair(name="q5", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))
    builder
  }

  it should "support the --min-reads option when one value is provided" in {
    caller(minReads=Seq(3)).consensusReadsFromSamRecords(builderForMinReads.toSeq) shouldBe empty
    caller(minReads=Seq(2)).consensusReadsFromSamRecords(builderForMinReads.toSeq) should have size 2
  }

  it should "support the --min-reads option when two values are provided" in {
    caller(minReads=Seq(5, 3)).consensusReadsFromSamRecords(builderForMinReads.toSeq) shouldBe empty
    caller(minReads=Seq(5, 2)).consensusReadsFromSamRecords(builderForMinReads.toSeq) should have size 2
  }

  it should "support the --min-reads option when three values are provided" in {
    caller(minReads=Seq(5, 3, 3)).consensusReadsFromSamRecords(builderForMinReads.toSeq) shouldBe empty
    caller(minReads=Seq(5, 3, 2)).consensusReadsFromSamRecords(builderForMinReads.toSeq) should have size 2
  }

  it should "not saturate the qualities with deep AB and light BA coverage" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    Range.inclusive(1, 50) foreach { i =>
      builder.addPair(name=s"a$i", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    }
    // Single BA read
    builder.addPair(name=s"b1", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    val recs = caller(post=45.toByte).consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("No first of pair."))
    val r2 = recs.find(_.secondOfPair).getOrElse(fail("No second of pair."))
    r1.basesString shouldBe "AAAAAAAAAA"
    r2.basesString shouldBe "GGGGGGGGGG" // after un-revcomping

    Seq(r1, r2) foreach { r =>
      r.quals.forall(_ <= (45+20)) shouldBe true
    }
  }

  it should "generate a sensible consensus when AB and BAs are sequenced/clipped to different lengths" in {
    val ab = new SamBuilder(readLength=10, baseQuality=20)
    val ba = new SamBuilder(readLength= 8, baseQuality=20)
    ab.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    ba.addPair(name="q2", start1=202, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCC",   bases2="AAAAAAAA",   attrs=Map(MI -> "foo/B"))

    val recs = c.consensusReadsFromSamRecords(ab.toSeq ++ ba.toSeq)
    recs should have size 2
    recs.foreach { r => r.length shouldBe 8 }
  }

  it should "respect minimum quality when building consensuses" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    builder.iterator.find(r => r.name == "q2" && r.secondOfPair).foreach { r => r.quals(5) = 5}

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("Missing R1"))
    val r2 = recs.find(_.secondOfPair).getOrElse(fail("Missing R2"))
    r1.basesString shouldBe "AAAAANAAAA"
    r1.quals(5) shouldBe PhredScore.MinValue
    r2.basesString shouldBe "GGGGGGGGGG"
  }

  it should "emit a no-call when there are conflicting SS consensus bases with equal quality" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="TAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("Missing R1"))
    r1.basesString shouldBe "NAAAAAAAAA"
    r1.quals(0) shouldBe PhredScore.MinValue
  }

  it should "emit a low qual call when there are conflicting SS consensus bases with unequal quality" in {
    val builder = new SamBuilder(readLength=10, baseQuality=30)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="TAAAAAAAAA", attrs=Map(MI -> "foo/B"))
    builder.find(r => r.name == "q2" && r.secondOfPair).foreach(_.quals(0) = 15)

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("Missing R1"))
    r1.basesString shouldBe "AAAAAAAAAA"
    r1.quals(0).toInt shouldBe 15 +- 5
  }

  Seq.range(1, PhredScore.MinValue+1).foreach { qualDiff =>
    it should s"emit a the min qual call when there are conflicting SS consensus bases when they differ by $qualDiff" in {
      val builder = new SamBuilder(readLength=10, baseQuality=30)
      builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
      builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="TAAAAAAAAA", attrs=Map(MI -> "foo/B"))
      builder.find(r => r.name == "q2" && r.secondOfPair).foreach(_.quals(0) = (30-qualDiff).toByte)

      val recs = c.consensusReadsFromSamRecords(builder.toSeq)
      recs should have size 2
      val r1 = recs.find(_.firstOfPair).getOrElse(fail("Missing R1"))
      r1.basesString shouldBe "NAAAAAAAAA"
      r1.quals(0).toInt shouldBe PhredScore.MinValue
    }
  }

  it should "count errors against the actual consensus base, before it is masked to N" in {
    val builder = new SamBuilder(readLength=10, baseQuality=30)
    builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="ANAAAAAAAA", bases2="CCNCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    val recs = c.consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val Seq(r1, r2) = recs
    r1[Float](RawReadErrorRate)    shouldBe 0f
    r1[Float](AbRawReadErrorRate)  shouldBe 0f
    r1[Float](BaRawReadErrorRate)  shouldBe 0f
    r2[Float](RawReadErrorRate)    shouldBe 0f
    r2[Float](AbRawReadErrorRate)  shouldBe 0f
    r2[Float](BaRawReadErrorRate)  shouldBe 0f
  }

  it should "generate correct summary and detail tags on consensus reads" in {
    val builder = new SamBuilder(readLength=10, baseQuality=30)
    // 3 AB pairs:
    //   First pair is good
    //   Second pair has R1 with a no-call at base index 2, and a low quality base a index 3
    //   Third pair  has one error at index 4 in each read
    builder.addPair(name="ab1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
    builder.addPair(name="ab2", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AANAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A")).foreach { r =>
      if (r.firstOfPair) {
        r.quals(2) = PhredScore.MinValue
        r.quals(3) = 10.toByte
      }
    }

    builder.addPair(name="ab3", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAGAAAAA", bases2="CCCCCTCCCC", attrs=Map(MI -> "foo/A"))

    // 2 BA pairs (recall that first of pair is R2 consensus and vice versa here):
    //    First pair is good
    //    Second pair is good
    builder.addPair(name="ba1", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))
    builder.addPair(name="ba2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))

    // Add the original UMI bases to each read
    builder.foreach { rec =>
      val mi = rec[String](MI)
      // all "/A" have the same UMI bases, while the "/B" bases have them swapped
      if (mi.endsWith("/A")) rec(RX) = "AAT-CCG" else rec(RX) = "CCG-AAT"
    }

    val recs = caller(q=20).consensusReadsFromSamRecords(builder.toSeq)
    recs should have size 2
    val r1 = recs.find(_.firstOfPair).getOrElse(fail("R1 missing."))
    val r2 = recs.find(_.secondOfPair).getOrElse(fail("R2 missing."))

    { // Check the per-read tags
      import ConsensusTags.PerRead._
      r1[Int](RawReadCount)          shouldBe 5
      r1[Int](AbRawReadCount)        shouldBe 3
      r1[Int](BaRawReadCount)        shouldBe 2
      r1[Int](MinRawReadCount)       shouldBe 4 // where the no-call or low-qual base are
      r1[Int](AbMinRawReadCount)     shouldBe 2 // where the no-call or low-qual base are
      r1[Int](BaMinRawReadCount)     shouldBe 2
      r1[Float](RawReadErrorRate)    shouldBe (1 / 48f) // 50 bases, minus 2 no/low-qual bases = 48
      r1[Float](AbRawReadErrorRate)  shouldBe (1 / 28f)
      r1[Float](BaRawReadErrorRate)  shouldBe 0f
      r1[String](MI) shouldBe "foo"
      r1[String](RX) shouldBe "AAT-CCG"

      r2[Int](RawReadCount)          shouldBe 5
      r2[Int](AbRawReadCount)        shouldBe 3
      r2[Int](BaRawReadCount)        shouldBe 2
      r2[Int](MinRawReadCount)       shouldBe 5
      r2[Int](AbMinRawReadCount)     shouldBe 3
      r2[Int](BaMinRawReadCount)     shouldBe 2
      r2[Float](RawReadErrorRate)    shouldBe (1 / 50f)
      r2[Float](AbRawReadErrorRate)  shouldBe (1 / 30f)
      r2[Float](BaRawReadErrorRate)  shouldBe 0f
      r2[String](MI) shouldBe "foo"
      r2[String](RX) shouldBe "AAT-CCG"
    }

    { // Check the per-base tags
      import ConsensusTags.PerBase._
      r1[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,2,2,3,3,3,3,3,3)
      r1[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
      r1[Array[Short]](BaRawReadCount)  shouldBe Array[Byte](2,2,2,2,2,2,2,2,2,2)
      r1[Array[Short]](BaRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,0,0,0,0,0)
      r1[String](AbConsensusBases)   shouldBe "AAAAAAAAAA"
      r1[String](BaConsensusBases)   shouldBe "AAAAAAAAAA"
      r1[String](AbConsensusQuals)   shouldBe "NNMMBNNNNN"
      r1[String](BaConsensusQuals)   shouldBe "MMMMMMMMMM"
      r2[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
      r2[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
      r2[Array[Short]](BaRawReadCount)  shouldBe Array[Byte](2,2,2,2,2,2,2,2,2,2)
      r2[Array[Short]](BaRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,0,0,0,0,0)
      r2[String](AbConsensusBases)   shouldBe "GGGGGGGGGG"
      r2[String](BaConsensusBases)   shouldBe "GGGGGGGGGG"
      r2[String](AbConsensusQuals)   shouldBe "NNNNBNNNNN"
      r2[String](BaConsensusQuals)   shouldBe "MMMMMMMMMM"
    }
  }

  it should "generate the correct per-base and per-read tags when only one strand is present" in {

    // /A strand
    {
      val builder = new SamBuilder(readLength=10, baseQuality=20)
      builder.addPair(name="q1", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
      builder.addPair(name="q2", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAAAAAAA", bases2="CCCCCCCCCC", attrs=Map(MI -> "foo/A"))
      builder.addPair(name="q3", start1=100, start2=200, strand1=Plus, strand2=Minus, bases1="AAAAACAAAA", bases2="CCCCCGCCCC", attrs=Map(MI -> "foo/A"))

      val recs = caller(minReads=Seq(3, 3, 0)).consensusReadsFromSamRecords(builder.toSeq)

      recs should have size 2
      val r1 = recs.find(_.firstOfPair).getOrElse(fail("R1 missing."))
      val r2 = recs.find(_.secondOfPair).getOrElse(fail("R2 missing."))

      r1.basesString shouldBe "AAAAAAAAAA"
      r2.basesString shouldBe "GGGGGGGGGG"

      { // Check the per-read tags
        import ConsensusTags.PerRead._
        r1[Int](RawReadCount)          shouldBe 3
        r1[Int](AbRawReadCount)        shouldBe 3
        r1[Int](BaRawReadCount)        shouldBe 0
        r1[Int](MinRawReadCount)       shouldBe 3
        r1[Int](AbMinRawReadCount)     shouldBe 3
        r1[Int](BaMinRawReadCount)     shouldBe 0
        r1[Float](RawReadErrorRate)    shouldBe 1/30f
        r1[Float](AbRawReadErrorRate)  shouldBe 1/30f
        r1[Float](BaRawReadErrorRate)  shouldBe 0f

        r2[Int](RawReadCount)          shouldBe 3
        r2[Int](AbRawReadCount)        shouldBe 3
        r2[Int](BaRawReadCount)        shouldBe 0
        r2[Int](MinRawReadCount)       shouldBe 3
        r2[Int](AbMinRawReadCount)     shouldBe 3
        r2[Int](BaMinRawReadCount)     shouldBe 0f
        r2[Float](RawReadErrorRate)    shouldBe 1/30f
        r2[Float](AbRawReadErrorRate)  shouldBe 1/30f
        r2[Float](BaRawReadErrorRate)  shouldBe 0f
      }

      { // Check the per-base tags
        import ConsensusTags.PerBase._
        r1[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
        r1[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,1,0,0,0,0)
        r1[String](AbConsensusBases)      shouldBe r1.basesString
        r1[String](AbConsensusQuals)      shouldBe "MMMMM9MMMM"
        r2[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
        r2[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
        r2[String](AbConsensusBases)      shouldBe r2.basesString
        r2[String](AbConsensusQuals)      shouldBe "MMMM9MMMMM"
        r1.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
        r1.get[Array[Short]](BaRawReadErrors).isEmpty shouldBe true
        r1.get[String](BaConsensusBases).isEmpty shouldBe true
        r1.get[String](BaConsensusQuals).isEmpty shouldBe true
        r2.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
        r2.get[Array[Short]](BaRawReadErrors).isEmpty shouldBe true
        r2.get[String](BaConsensusBases).isEmpty shouldBe true
        r2.get[String](BaConsensusQuals).isEmpty shouldBe true
      }
    }

    // /B strand
    {
      val builder = new SamBuilder(readLength=10, baseQuality=20)
      builder.addPair(name="q1", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))
      builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCCCCCC", bases2="AAAAAAAAAA", attrs=Map(MI -> "foo/B"))
      builder.addPair(name="q2", start1=200, start2=100, strand1=Minus, strand2=Plus, bases1="CCCCCGCCCC", bases2="AAAAATAAAA", attrs=Map(MI -> "foo/B"))

      val recs = caller(minReads=Seq(3, 3, 0)).consensusReadsFromSamRecords(builder.toSeq)

      recs should have size 2
      val r1 = recs.find(_.firstOfPair).getOrElse(fail("R1 missing."))
      val r2 = recs.find(_.secondOfPair).getOrElse(fail("R2 missing."))

      r1.basesString shouldBe "GGGGGGGGGG"
      r2.basesString shouldBe "AAAAAAAAAA"

      { // Check the per-read tags
        import ConsensusTags.PerRead._
        r1[Int](RawReadCount)          shouldBe 3
        r1[Int](AbRawReadCount)        shouldBe 3
        r1[Int](BaRawReadCount)        shouldBe 0
        r1[Int](MinRawReadCount)       shouldBe 3
        r1[Int](AbMinRawReadCount)     shouldBe 3
        r1[Int](BaMinRawReadCount)     shouldBe 0
        r1[Float](RawReadErrorRate)    shouldBe 1/30f
        r1[Float](AbRawReadErrorRate)  shouldBe 1/30f
        r1[Float](BaRawReadErrorRate)  shouldBe 0f

        r2[Int](RawReadCount)          shouldBe 3
        r2[Int](AbRawReadCount)        shouldBe 3
        r2[Int](BaRawReadCount)        shouldBe 0
        r2[Int](MinRawReadCount)       shouldBe 3
        r2[Int](AbMinRawReadCount)     shouldBe 3
        r2[Int](BaMinRawReadCount)     shouldBe 0f
        r2[Float](RawReadErrorRate)    shouldBe 1/30f
        r2[Float](AbRawReadErrorRate)  shouldBe 1/30f
        r2[Float](BaRawReadErrorRate)  shouldBe 0f
      }

      { // Check the per-base tags
        import ConsensusTags.PerBase._
        r1[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
        r1[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,1,0,0,0,0,0)
        r1[String](AbConsensusBases)      shouldBe r1.basesString
        r1[String](AbConsensusQuals)      shouldBe "MMMM9MMMMM"
        r2[Array[Short]](AbRawReadCount)  shouldBe Array[Byte](3,3,3,3,3,3,3,3,3,3)
        r2[Array[Short]](AbRawReadErrors) shouldBe Array[Byte](0,0,0,0,0,1,0,0,0,0)
        r2[String](AbConsensusBases)      shouldBe r2.basesString
        r2[String](AbConsensusQuals)      shouldBe "MMMMM9MMMM"
        r1.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
        r1.get[Array[Short]](BaRawReadErrors).isEmpty shouldBe true
        r1.get[String](BaConsensusBases).isEmpty shouldBe true
        r1.get[String](BaConsensusQuals).isEmpty shouldBe true
        r2.get[Array[Short]](BaRawReadCount).isEmpty shouldBe true
        r2.get[Array[Short]](BaRawReadErrors).isEmpty shouldBe true
        r2.get[String](BaConsensusBases).isEmpty shouldBe true
        r2.get[String](BaConsensusQuals).isEmpty shouldBe true
      }
    }
  }
}
