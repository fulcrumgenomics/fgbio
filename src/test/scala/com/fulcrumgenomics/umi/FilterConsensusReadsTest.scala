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

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.SamPairUtil

class FilterConsensusReadsTest extends UnitSpec {
  private val temp = makeTempFile("meh.", ".bam")
  private val Q2 = PhredScore.MinValue

  /** Make a reference file that is 100 lines of 100 As. */
  private val ref = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 10 * 100) // 5000 bases
    builder.toTempFile()
  }

  /** Generates a FilterConsensusReads instance for testing non-end-to-end methods on vanilla reads. */
  private def fv(q: Int, d: Int, mq: Option[PhredScore], readErr: Double, baseErr: Double, nf: Double) : FilterConsensusReads =
    new FilterConsensusReads(input=temp, output=temp, ref=temp, reversePerBaseTags=false,
      minBaseQuality=q.toByte, minReads=Seq(d), minMeanBaseQuality=mq,
      maxReadErrorRate=Seq(readErr), maxBaseErrorRate=Seq(baseErr),
      maxNoCallFraction=nf)

  /** Generates a FilterConsensusReads instance for testing non-end-to-end methods on duplex reads. */
  private def fv(q: Int, d: Seq[Int], mq: Option[PhredScore], readErr: Seq[Double], baseErr: Seq[Double], nf: Double, ss: Boolean) : FilterConsensusReads =
    new FilterConsensusReads(input=temp, output=temp, ref=temp, reversePerBaseTags=false,
      minBaseQuality=q.toByte, minReads=d, minMeanBaseQuality=mq,
      maxReadErrorRate=readErr, maxBaseErrorRate=baseErr, maxNoCallFraction=nf, requireSingleStrandAgreement=ss)

  /** Tags up a SAMRecord for filtering. */
  private def tag(rec: SamRecord, minDepth: Int, depth: Int, readErr: Float, depths: Array[Short]=null, errors: Array[Short]=null): SamRecord = {
    rec.bases = "A" * rec.length
    rec(ConsensusTags.PerRead.RawReadCount)     = depth
    rec(ConsensusTags.PerRead.MinRawReadCount)  = minDepth
    rec(ConsensusTags.PerRead.RawReadErrorRate) = readErr
    rec(ConsensusTags.PerBase.RawReadCount)     = depths
    rec(ConsensusTags.PerBase.RawReadErrors)    = errors
    rec
  }

  /** Makes a SAMRecord with which to test filtering. */
  private def r(len: Int=10, q: Int=40, minDepth: Int, depth: Int, readErr: Float, depths: Array[Short]=null, errors: Array[Short]=null) : SamRecord = {
    val builder = new SamBuilder(readLength=len, baseQuality=q.toByte)
    tag(builder.addFrag(start=100).get, minDepth=minDepth, depth=depth, readErr=readErr, depths=depths, errors=errors)
  }

  /** Makes an array of Shorts of all the same value. */
  private def arr(value: Short, length: Int): Array[Short] = {
    val ss = new Array[Short](length)
    forloop (from=0, until=length) { i => ss(i) = value }
    ss
  }

  "FilterConsensusReads.filterRead" should "keep rec and not mask any bases" in {
    val filter = fv(q=10, d=2, mq=None, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=2, minDepth=2, readErr = 0.01f)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 0
  }

  it should "keep rec and not mask any bases when detail tags are present" in {
    val filter = fv(q=10, d=2, mq=None, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=2, minDepth=2, readErr=0.01f, depths=arr(2, 10), errors=arr(0, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 0
    rec.basesString shouldBe ("A" * 10)
  }

  it should "filter out a read with depth < minDepth" in {
    val filter = fv(q=10, d=5, mq=None, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=4, minDepth=4, readErr=0.00f, depths=arr(4, 10), errors=arr(0, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.basesString shouldBe ("A" * 10)
  }

  it should "filter out a read with too high an error rate" in {
    val filter = fv(q=10, d=3, mq=None, readErr=0.05, baseErr=0.5, nf=0.2)
    val rec    = r(depth=50, minDepth=50, readErr=0.20f, depths=arr(50, 10), errors=arr(10, 10))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.basesString shouldBe ("A" * 10)
  }

  it should "filter out bases that are below the required quality score" in {
    val filter = fv(q=40, d=3, mq=None, readErr=0.05, baseErr=0.1, nf=0.5)
    val rec    = r(depth=5, minDepth=5, readErr=0.00f, depths=arr(5, 10), errors=arr(0, 10))
    Seq(1,4,7).foreach(i => rec.quals(i) = 20)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 3
    rec.basesString shouldBe "ANAANAANAA"
  }

  it should "filter out bases that have too much disagreement (errors/depth > baseErr)" in {
    val filter = fv(q=40, d=10, mq=None, readErr=0.05, baseErr=0.1, nf=0.5)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f, depths=arr(20, 10), errors=Array[Short](0,1,2,3,4,3,2,1,0,0))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 3
    rec.basesString shouldBe "AAANNNAAAA"
  }

  it should "filter out the read if too many bases are masked in the read in the first place" in {
    val filter = fv(q=20, d=5, mq=None, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f)
    rec.bases = "AANNNNNAAA"
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.basesString shouldBe "AANNNNNAAA"
  }

  it should "filter out the read if too many bases are masked post-filtering" in {
    val filter = fv(q=20, d=5, mq=None, readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(depth=20, minDepth=20, readErr=0.00f, depths=arr(20, 10), errors=Array[Short](3,3,3,0,0,0,0,3,3,3))
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 6
    rec.basesString shouldBe "NNNAAAANNN"
  }

  it should "filter out the read if the mean base quality is too low" in {
    val filter = fv(q=20, d=5, mq=Some(31.toByte), readErr=0.05, baseErr=0.1, nf=0.2)
    val rec    = r(q=30, depth=20, minDepth=20, readErr=0.00f)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe false
    result.maskedBases shouldBe 0
    rec.basesString shouldBe ("A" * 10)
  }

  "FilterConsensusReads" should "not filter out any reads" in {
    val builder = new SamBuilder(readLength=10, baseQuality=45, sort=Some(SamOrder.Queryname))
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q2", start1=100, start2=200).foreach(r => tag(r, minDepth=5, depth=5, readErr=0f, depths=arr(5, 10), errors=arr(0,10)))
    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=Seq(3), maxReadErrorRate=Seq(0.025), maxBaseErrorRate=Seq(0.1), maxNoCallFraction=0.1).execute()

    val recs = SamSource(out).toSeq
    recs.size shouldBe 4
    recs.exists(_.basesString.contains("N")) shouldBe false
  }

  it should "filter out both reads of a pair whenever one fails" in {
    val builder = new SamBuilder(readLength=10, baseQuality=45, sort=Some(SamOrder.Queryname))
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q2", start1=100, start2=200).foreach(r => tag(r, minDepth=5, depth=5, readErr=0f, depths=arr(5, 10), errors=arr(0,10)))
    builder.filter(_.name == "q1").filter(_.firstOfPair).foreach(r => r(ConsensusTags.PerRead.RawReadErrorRate) =  0.2f)

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=Seq(3), maxReadErrorRate=Seq(0.025), maxBaseErrorRate=Seq(0.1), maxNoCallFraction=0.1).execute()

    val recs = SamSource(out).toSeq
    recs.size shouldBe 2
    recs.map(_.name).distinct shouldBe Seq("q2")
  }

  it should "filter out all reads with the same query name when any primary read fails" in {
    val builder = new SamBuilder(readLength=10, baseQuality=45, sort=Some(SamOrder.Queryname))
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q1", start1=100, start2=200).foreach { r =>
      tag(r, minDepth = 5, depth = 5, readErr = 0f, depths = arr(5, 10), errors = arr(0, 10))
      r.supplementary = true
    }

    builder.filterNot(_.supplementary).filter(_.firstOfPair).foreach(r => r(ConsensusTags.PerRead.RawReadErrorRate) = 0.2f)

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=Seq(3), maxReadErrorRate=Seq(0.025), maxBaseErrorRate=Seq(0.1), maxNoCallFraction=0.1).execute()

    val recs = SamSource(out).toSeq
    recs.size shouldBe 0
  }

  it should "not filter out primary alignments or other supplementals when a supplemental or secondary fails" in {
    val builder = new SamBuilder(readLength=10, baseQuality=45, sort=Some(SamOrder.Queryname))
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=4, depth=4, readErr=0f, depths=arr(4, 10), errors=arr(0,10)))
    builder.addPair(name="q1", start1=110, start2=190).foreach { r =>
      tag(r, minDepth = 5, depth = 5, readErr = 0f, depths = arr(5, 10), errors = arr(0, 10))
      r.supplementary = true
    }

    builder.filter(_.supplementary).filter(_.firstOfPair).foreach(r => r(ConsensusTags.PerRead.RawReadErrorRate) = 0.2f)

    val out = makeTempFile("filtered.", ".bam")
    new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=false,
      minBaseQuality=45.toByte, minReads=Seq(3), maxReadErrorRate=Seq(0.025), maxBaseErrorRate=Seq(0.1), maxNoCallFraction=0.1).execute()

    val recs = SamSource(out).toSeq
    recs.size shouldBe 3
    recs.count(r => !r.secondary && !r.supplementary) shouldBe 2
    recs.count(r => r.secondary || r.supplementary) shouldBe 1
  }

  it should "reverse the per-base tags only when requested" in {
    val builder = new SamBuilder(readLength=10, baseQuality=45, sort=Some(SamOrder.Queryname))
    builder.addPair(name="q1", start1=100, start2=200).foreach(r => tag(r, minDepth=7, depth=9, readErr=0f, depths=Array[Short](7,7,7,8,8,8,9,9,9,9), errors=arr(0,10)))

    Seq(true, false) foreach { reverseTags =>
      val out = makeTempFile("filtered.", ".bam")
      new FilterConsensusReads(input=builder.toTempFile(), output=out, ref=ref, reversePerBaseTags=reverseTags,
        minBaseQuality=45.toByte, minReads=Seq(3), maxReadErrorRate=Seq(0.025), maxBaseErrorRate=Seq(0.1), maxNoCallFraction=0.1).execute()

      val recs = SamSource(out).toSeq
      recs.size shouldBe 2
      recs.exists(_.negativeStrand) shouldBe true
      recs.foreach { rec =>
        val depths = rec[Array[Short]](ConsensusTags.PerBase.RawReadCount)
        if (reverseTags && rec.negativeStrand) depths shouldBe Array[Short](9,9,9,9,8,8,8,7,7,7) else depths shouldBe Array[Short](7,7,7,8,8,8,9,9,9,9)
      }
    }
  }

  it should "reverse the per-base duplex tags only when requested" in {
    val path = Files.createTempFile("SamRecordSet.", ".bam")

    {
      val good1 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.01f, err=0.01f,
      abBases="A"*5 + "C"*5, baBases="A"*5 + "C"*5, abQuals="I"*5 + "H"*5, baQuals="I"*5 + "H"*5)
      val good2 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.01f, err=0.01f,
        abBases="A"*5 + "C"*5, baBases="A"*5 + "C"*5, abQuals="I"*5 + "H"*5, baQuals="I"*5 + "H"*5)
      good1.paired = true
      good1.firstOfPair = true
      good1.secondOfPair = false
      good2.paired = true
      good2.firstOfPair = false
      good2.secondOfPair = true
      good2.negativeStrand = true
      good2.name = good1.name
      SamPairUtil.setMateInfo(good1.asSam, good2.asSam)
      path.toFile.deleteOnExit()
      val writer = SamWriter(path, DuplexBuilder.header, sort=DuplexBuilder.sort)
      writer += good1
      writer += good2
      writer.close()
    }

    Seq(true, false) foreach { reverseTags =>
      val out = makeTempFile("filtered.", ".bam")
      new FilterConsensusReads(input=path, output=out, ref=ref, reversePerBaseTags=reverseTags,
        minBaseQuality=2.toByte, minReads=Seq(0), maxReadErrorRate=Seq(1), maxBaseErrorRate=Seq(1), maxNoCallFraction=1).execute()

      val recs = SamSource(out).toSeq
      recs.size shouldBe 2
      recs.exists(_.negativeStrand) shouldBe true
      recs.foreach { rec =>
        if (reverseTags && rec.negativeStrand) {
          val abBases = rec[String](ConsensusTags.PerBase.AbConsensusBases)
          val abQuals = rec[String](ConsensusTags.PerBase.AbConsensusQuals)
          abBases shouldBe ("G"*5 + "T"*5)
          abQuals shouldBe ("H"*5 + "I"*5)
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Below this line are tests for filtering of duplex consensus reads.
  //////////////////////////////////////////////////////////////////////////////

  /** The type of all the detail tags. */
  type Arr = Array[Short]

  /** Builds a short[] that is of length n and filled value s. */
  def ss(s: Short, n: Int = 10): Arr = {
    val arr = new Array[Short](n)
    arr.indices.foreach(i => arr(i) = s)
    arr
  }

  /** An extension to SamBuilder to make building duplex pairs marginally less painful. */
  object DuplexBuilder extends SamBuilder(readLength=10, baseQuality=90) {
    /** Method to create/add a single read with all the duplex tags on it. */
    def add(name:String=nextName,
            dp: Int=10, minDp:Int=10, abDp:Int=5, baDp:Int=5, abMinDp:Int=5, baMinDp:Int=5,
            err: Float=0, abErr:Float=0, baErr:Float=0,
            abDps:Arr=ss(5), baDps:Arr=ss(5), abErrs:Arr=ss(0), baErrs:Arr=ss(0),
            abBases: String="A"*10, baBases: String="A"*10, abQuals:String="I"*10, baQuals:String="I"*10) = {

      val r1 = addFrag(name = name, start=100).get
      r1.bases = "AAAAAAAAAA"
      r1.paired = true
      r1.firstOfPair = true
      r1(ConsensusTags.PerRead.RawReadCount) = dp
      r1(ConsensusTags.PerRead.AbRawReadCount) = abDp
      r1(ConsensusTags.PerRead.BaRawReadCount) = baDp
      r1(ConsensusTags.PerRead.MinRawReadCount) = minDp
      r1(ConsensusTags.PerRead.AbMinRawReadCount) = abMinDp
      r1(ConsensusTags.PerRead.BaMinRawReadCount) = baMinDp
      r1(ConsensusTags.PerRead.RawReadErrorRate) = err
      r1(ConsensusTags.PerRead.AbRawReadErrorRate) = abErr
      r1(ConsensusTags.PerRead.BaRawReadErrorRate) = baErr
      r1(ConsensusTags.PerBase.AbRawReadCount) = abDps
      r1(ConsensusTags.PerBase.BaRawReadCount) = baDps
      r1(ConsensusTags.PerBase.AbRawReadErrors) = abErrs
      r1(ConsensusTags.PerBase.BaRawReadErrors) = baErrs
      r1(ConsensusTags.PerBase.AbConsensusBases) = abBases
      r1(ConsensusTags.PerBase.BaConsensusBases) = baBases
      r1(ConsensusTags.PerBase.AbConsensusQuals) = abQuals
      r1(ConsensusTags.PerBase.BaConsensusQuals) = baQuals
      r1
    }

    /** Simpler method that adds a duplex read where all bases have the same depth and no errors. */
    def addSimple(name:String=nextName, abDp:Int=5, baDp:Int=5): SamRecord = {
      add(name=name, dp=abDp+baDp, abDp=abDp, baDp=baDp, abDps=ss(abDp.toShort), baDps=ss(baDp.toShort))
    }
  }

  it should "throw validation exceptions if you provide options in increasing stringency order" in {
    val Seq(in, out) = Seq("in.", "out.").map(makeTempFile(_, ".bam"))
    an[Exception] shouldBe thrownBy { new FilterConsensusReads(input=in, output=out, ref=ref, minBaseQuality=Q2, minReads=Seq(1,2,3)) }
    an[Exception] shouldBe thrownBy { new FilterConsensusReads(input=in, output=out, ref=ref, minBaseQuality=Q2, minReads=Seq(9,4,6)) }
    an[Exception] shouldBe thrownBy { new FilterConsensusReads(input=in, output=out, ref=ref, minBaseQuality=Q2, minReads=Seq(1), maxReadErrorRate=Seq(0.1f, 0.2f, 0.01f)) }
    an[Exception] shouldBe thrownBy { new FilterConsensusReads(input=in, output=out, ref=ref, minBaseQuality=Q2, minReads=Seq(1), maxBaseErrorRate=Seq(0.1f, 0.2f, 0.01f)) }
  }

  it should "set filter values correctly when only one value is provided" in {
    val Seq(in, out) = Seq("in.", "out.").map(makeTempFile(_, ".bam"))
    val filter = new FilterConsensusReads(input=in, output=out, ref=ref,
      minReads=Seq(99), minBaseQuality=88.toByte, maxReadErrorRate=Seq(0.11f), maxBaseErrorRate=Seq(0.22f))

    Seq(filter.ccFilters, filter.abFilters, filter.baFilters).foreach { fs =>
      fs.minReads         shouldBe 99
      fs.maxReadErrorRate shouldBe 0.11f
      fs.maxBaseErrorRate shouldBe 0.22f
    }
  }

  it should "set filter values correctly when two values are provided" in {
    val Seq(in, out) = Seq("in.", "out.").map(makeTempFile(_, ".bam"))
    val filter = new FilterConsensusReads(input=in, output=out, ref=ref,
      minReads=Seq(99,98), minBaseQuality=88.toByte, maxReadErrorRate=Seq(0.11f, 0.10f), maxBaseErrorRate=Seq(0.22f, 0.2f))

    filter.ccFilters.minReads         shouldBe 99
    filter.ccFilters.maxReadErrorRate shouldBe 0.11f
    filter.ccFilters.maxBaseErrorRate shouldBe 0.22f

    Seq(filter.abFilters, filter.baFilters).foreach { fs =>
      fs.minReads         shouldBe 98
      fs.maxReadErrorRate shouldBe 0.10f
      fs.maxBaseErrorRate shouldBe 0.20f
    }
  }

  it should "set filter values correctly when three values are provided" in {
    val Seq(in, out) = Seq("in.", "out.").map(makeTempFile(_, ".bam"))
    val filter = new FilterConsensusReads(input=in, output=out, ref=ref,
      minReads=Seq(99,98,97), minBaseQuality=88.toByte, maxReadErrorRate=Seq(0.11f, 0.05f, 0.1f), maxBaseErrorRate=Seq(0.22f, 0.1f, 0.2f))

    filter.ccFilters.minReads         shouldBe 99
    filter.ccFilters.maxReadErrorRate shouldBe 0.11f
    filter.ccFilters.maxBaseErrorRate shouldBe 0.22f
    filter.abFilters.minReads         shouldBe 98
    filter.abFilters.maxReadErrorRate shouldBe 0.05f
    filter.abFilters.maxBaseErrorRate shouldBe 0.1f
    filter.baFilters.minReads         shouldBe 97
    filter.baFilters.maxReadErrorRate shouldBe 0.1f
    filter.baFilters.maxBaseErrorRate shouldBe 0.2f
  }

  "FilterConsensusReads.filterRead(duplex)" should "filter appropriately by total read depth" in {
    val good = DuplexBuilder.addSimple(abDp=3, baDp=2)
    val bad  = DuplexBuilder.addSimple(abDp=2, baDp=2)
    val filter = fv(q=1, d=Seq(5,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad ).keepRead shouldBe false
  }

  it should "filter appropriately by read error rate" in {
    val good = DuplexBuilder.add(err=0.024f)
    val bad  = DuplexBuilder.add(err=0.026f)
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(0.025,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad ).keepRead shouldBe false
  }

  it should "filter appropriately by AB/BA read depth when they are set to the same cutoff value" in {
    val good = DuplexBuilder.addSimple(abDp=3, baDp=3)
    val bad1 = DuplexBuilder.addSimple(abDp=3, baDp=2)
    val bad2 = DuplexBuilder.addSimple(abDp=2, baDp=3)
    val filter = fv(q=1, d=Seq(3,3,3), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
    filter.filterRecord(bad2).keepRead shouldBe false
  }

  it should "filter appropriately by AB/BA read depth when they are different cutoff values" in {
    val good1 = DuplexBuilder.addSimple(abDp=3, baDp=3)
    val good2 = DuplexBuilder.addSimple(abDp=3, baDp=2)
    val good3 = DuplexBuilder.addSimple(abDp=2, baDp=3)
    val bad1  = DuplexBuilder.addSimple(abDp=2, baDp=2)
    val filter = fv(q=1, d=Seq(3,3,2), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good1).keepRead shouldBe true
    filter.filterRecord(good2).keepRead shouldBe true
    filter.filterRecord(good3).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
  }

  it should "filter appropriately by AB/BA read error rate when they are set to the same cutoff value" in {
    val good = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.01f, err=0.01f)
    val bad1 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.03f, err=0.02f)
    val bad2 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.03f, baErr=0.01f, err=0.02f)
    val bad3 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.03f, baErr=0.03f, err=0.03f)
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,0.02,0.02), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
    filter.filterRecord(bad2).keepRead shouldBe false
    filter.filterRecord(bad3).keepRead shouldBe false
  }

  it should "filter appropriately by AB/BA read error rate when they are set to different cutoff values" in {
    val good1 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.01f, err=0.01f)
    val good2 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.03f, err=0.02f)
    val good3 = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.03f, err=0.02f)
    val bad1  = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.01f, baErr=0.04f, err=0.03f)
    val bad2  = DuplexBuilder.add(abDp=100, baDp=100, abErr=0.04f, baErr=0.01f, err=0.03f)
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,0.02,0.03), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good1).keepRead shouldBe true
    filter.filterRecord(good2).keepRead shouldBe true
    filter.filterRecord(good3).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
    filter.filterRecord(bad2).keepRead shouldBe false
  }

  it should "filter appropriately by the mean quality" in {
    val good1 = {
      val r = DuplexBuilder.add(dp=100, abDp=50, baDp=50, abDps=ss(50), baDps=ss(50))
      r.quals = ss(50, r.length).map(_.toByte)
      r
    }
    val good2 = {
      val r = DuplexBuilder.add(dp=100, abDp=50, baDp=50, abDps=ss(50), baDps=ss(50))
      r.quals = ss(5, r.length).map(_.toByte)
      r
    }
    val bad1 = {
      val r = DuplexBuilder.add(dp=100, abDp=50, baDp=50, abDps=ss(50), baDps=ss(50))
      r.quals = ss(4, r.length).map(_.toByte)
      r
    }
    val bad2 = {
      val r = DuplexBuilder.add(dp=100, abDp=50, baDp=50, abDps=ss(50), baDps=ss(50))
      r.quals = ss(5, r.length).map(_.toByte).drop(1) ++ Array[Byte](4.toByte)
      r
    }
    val filter = fv(q=1, d=Seq(1,1,1), mq=Some(5.toByte), readErr=Seq(1.0,0.02,0.03), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    filter.filterRecord(good1).keepRead shouldBe true
    filter.filterRecord(good2).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
    filter.filterRecord(bad2).keepRead shouldBe false
  }

  it should "mask bases when the total depth dips below the required minimum" in {
    val rec = DuplexBuilder.add(dp=15, abDp=8, baDp=7, abDps=Array[Short](8,7,8,8,8,8,8,8,8,8), baDps=Array[Short](7,7,7,7,2,7,7,7,1,7))
    val filter = fv(q=1, d=Seq(15,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 3
    rec.basesString shouldBe "AnAAnAAAnA".toUpperCase
  }

  it should "mask bases when the A/B depth dips below the symmetric required minimum" in {
    val rec = DuplexBuilder.add(dp=6, abDp=3, baDp=3, abDps=Array[Short](3,3,2,3,3,1,3,3,3,2), baDps=Array[Short](3,3,3,3,3,3,3,3,2,2))
    val filter = fv(q=1, d=Seq(3,3,3), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 4
    rec.basesString shouldBe "AAnAAnAAnn".toUpperCase
  }

  it should "mask bases when the A/B depth dips below the asymmetric required minimums" in {
    val rec = DuplexBuilder.add(dp=6, abDp=3, baDp=3, abDps=Array[Short](2,2,2,3,3,3,3,3,2,3), baDps=Array[Short](3,3,3,2,2,2,3,3,2,1))
    val filter = fv(q=1, d=Seq(3,3,2), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0,1.0,1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 2
    rec.basesString shouldBe "AAAAAAAAnn".toUpperCase
  }

  it should "mask bases when the total error rate rises above the required per-base maximum" in {
    val rec = DuplexBuilder.add(dp=100, abDp=50, baDp=50, abDps=ss(50), baDps=ss(50),
      abErrs=Array[Short](5,0,0,0,4,0,3,0,2,0), baErrs=Array[Short](0,0,5,0,1,0,2,0,3,0) )
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(0.04, 1.0, 1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 5
    rec.basesString shouldBe "nAnAnAnAnA".toUpperCase
  }

  it should "mask bases when the A/B error rates rise above the required symmetric maximum" in {
    val rec = DuplexBuilder.add(dp=200, abDp=100, baDp=100, abDps=ss(100), baDps=ss(100),
      abErrs=Array[Short](0,1,2,3,0,0,0,1,2,3), baErrs=Array[Short](0,3,2,1,0,0,0,3,2,3) )
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 0.025, 0.025), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 4
    rec.basesString shouldBe "AnAnAAAnAn".toUpperCase
  }

  it should "mask bases when the A/B error rates rise above the required asymmetric maximums" in {
    val rec = DuplexBuilder.add(dp=200, abDp=100, baDp=100, abDps=ss(100), baDps=ss(100),
      abErrs=Array[Short](0,1,2,3,0,3,3,3,2,1), baErrs=Array[Short](0,1,2,3,0,1,2,3,3,3) )
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 0.02, 0.03), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead shouldBe true
    result.maskedBases shouldBe 2
    rec.basesString shouldBe "AAAnAAAnAA".toUpperCase
  }
  
  it should "mask bases by base quality" in {
    val rec = DuplexBuilder.addSimple(abDp=3, baDp=3)
    rec.quals = Array(10, 20, 30, 40, 50, 60, 70, 80, 90, 99).map(_.toByte)
    val filter = fv(q=80, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead    shouldBe true
    result.maskedBases shouldBe 7
    rec.basesString  shouldBe "nnnnnnnAAA".toUpperCase
  }

  it should "mask bases if the AB and BA consensus reads disagree" in {
    val rec = DuplexBuilder.addSimple(abDp=3, baDp=3)
    rec(ConsensusTags.PerBase.AbConsensusBases) = rec.basesString
    rec(ConsensusTags.PerBase.BaConsensusBases) = "NT" + rec.basesString.drop(2)
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=1.0, ss=true)
    val result = filter.filterRecord(rec)
    result.keepRead    shouldBe true
    result.maskedBases shouldBe 2
    rec.basesString  shouldBe "nnAAAAAAAA".toUpperCase
  }

  it should "not mask bases if the AB and BA consensus reads disagree but single-strand-agreement is not specified" in {
    val rec = DuplexBuilder.addSimple(abDp=3, baDp=3)
    rec(ConsensusTags.PerBase.AbConsensusBases) = rec.basesString
    rec(ConsensusTags.PerBase.BaConsensusBases) = "NT" + rec.basesString.drop(2)
    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=1.0, ss=false)
    val result = filter.filterRecord(rec)
    result.keepRead    shouldBe true
    result.maskedBases shouldBe 0
    rec.basesString  shouldBe "AAAAAAAAAA".toUpperCase
  }

  it should "not mask bases if missing an AB or BA and requiring single strand agreement" in {
    val rec = DuplexBuilder.addSimple(abDp=3, baDp=0)
    rec(ConsensusTags.PerBase.AbConsensusBases) = rec.basesString
    val filter = fv(q=1, d=Seq(1,1,0), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=1.0, ss=true)
    val result = filter.filterRecord(rec)
    result.keepRead    shouldBe true
    result.maskedBases shouldBe 0
    rec.basesString  shouldBe "AAAAAAAAAA"
  }

  it should "reject reads that have too many Ns when the come in" in {
    val good = DuplexBuilder.addSimple(abDp=3, baDp=3)
    good.bases = "AAAAAAAANN"
    val bad  = DuplexBuilder.addSimple(abDp=3, baDp=3)
    bad.bases = "NNAAAAAANN"

    val filter = fv(q=1, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=0.2, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad ).keepRead shouldBe false
  }

  it should "reject reads that have too many Ns after masking" in {
    val good = DuplexBuilder.addSimple(abDp=3, baDp=3)
    val bad1 = DuplexBuilder.addSimple(abDp=3, baDp=3)
    val bad2 = DuplexBuilder.addSimple(abDp=3, baDp=3)

    good.bases = "AAAAAAAANN"
    bad1.bases = "NNAAAAAAAA"; bad1.quals = Array( 2, 2,20,30,90,90,90,90,90,90).map(_.toByte)
    bad2.bases = "AAAAAAAAAA"; bad2.quals = Array(90,90,40,35,48,90,90,90,90,90).map(_.toByte)

    val filter = fv(q=50, d=Seq(1,1,1), mq=None, readErr=Seq(1.0,1.0,1.0), baseErr=Seq(1.0, 1.0, 1.0), nf=0.2, ss=false)
    filter.filterRecord(good).keepRead shouldBe true
    filter.filterRecord(bad1).keepRead shouldBe false
    filter.filterRecord(bad2).keepRead shouldBe false
  }
}
