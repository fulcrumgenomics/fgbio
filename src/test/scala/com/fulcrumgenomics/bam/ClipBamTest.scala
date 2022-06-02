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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.bam.ClippingMode.{Hard, Soft, SoftWithMask}
import com.fulcrumgenomics.bam.api.SamOrder.Queryname
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.testing.SamBuilder._
import com.fulcrumgenomics.testing.{ErrorLogLevel, ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.SAMTag
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.SequenceUtil
import org.scalatest.OptionValues

import scala.util.Random

class ClipBamTest extends UnitSpec with ErrorLogLevel with OptionValues {
  private val ref = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 500) // 5000 bases
    builder.add("chr2").add("CCCCCCCCCC", 500) // 5000 bases
    builder.toTempFile()
  }

  private val dummyBam = makeTempFile("dummy.", ".bam")

  private object StartAndEnd {
    def apply(r: SamRecord): StartAndEnd = StartAndEnd(start=r.start, end=r.end)
  }

  private case class StartAndEnd(start: Int, end: Int) {
    // NB: assumes FR reads
    def checkClipping(r: SamRecord, fivePrime: Int, threePrime: Int): Unit = {
      if (r.negativeStrand) {
        r.start shouldBe (this.start + threePrime)
        r.end shouldBe (this.end - fivePrime)
      }
      else {
        r.start shouldBe (this.start + fivePrime)
        r.end shouldBe (this.end - threePrime)
      }
    }
  }

  private case class StartsAndEnds(r1: SamRecord, r2: SamRecord) {
    private val prior1 = StartAndEnd(r1)
    private val prior2 = StartAndEnd(r2)
    // NB: assumes FR reads
    def checkClipping(r1: SamRecord, r2: SamRecord,
                      readOneFivePrime: Int, readOneThreePrime: Int,
                      readTwoFivePrime: Int, readTwoThreePrime: Int): Unit = {
      prior1.checkClipping(r=r1, readOneFivePrime, readOneThreePrime)
      prior2.checkClipping(r=r2, readTwoFivePrime, readTwoThreePrime)
    }
  }

  "ClipBam.clip" should "should fail when no clipping options are given" in {
    an[Exception] should be thrownBy new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
  }

  it should "not clip reads where either read is unaligned" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, unmapped2=true)
    val expected = r1.cigar
    clipper.clipPair(r1, r2)
    r1.cigar shouldBe expected
  }

  it should "not clip reads that are on different chromosomes" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)
    r2.refIndex = 1
    r1.mateRefIndex = 1

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clipPair(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "not clip reads that are abutting but not overlapped" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=150)

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clipPair(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "not clip non-FR reads " in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, strand2=Plus)

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clipPair(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "clip reads that are fully overlapped" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)

    clipper.clipPair(r1, r2)
    r1.cigar.toString shouldBe "25M25H"
    r2.cigar.toString shouldBe "25H25M"
  }

  it should "clip reads that are overlapped by just one base" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=149)

    r1.end shouldBe r2.start
    clipper.clipPair(r1, r2)
    r1.end shouldBe (r2.start - 1)
  }

  it should "handle reads that contain insertions" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2I8M", cigar2="10M2I38M")

    r1.end >= r2.start shouldBe true
    clipper.clipPair(r1, r2)
    r1.end >= r2.start shouldBe false
  }

  it should "handle reads that contain deletions" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2D10M", cigar2="10M2D40M")

    r1.end >= r2.start shouldBe true
    clipper.clipPair(r1, r2)
    r1.end >= r2.start shouldBe false
  }

  it should "clip a fixed amount on the ends of the reads with reads that do not overlap" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
      readOneFivePrime=1, readOneThreePrime=2, readTwoFivePrime=3, readTwoThreePrime=4)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=150)

    val prior = StartsAndEnds(r1, r2)
    r1.end shouldBe (r2.start - 1)
    clipper.clipPair(r1, r2)
    prior.checkClipping(r1, r2, 1, 2, 3, 4)
  }

  it should "clip a fixed amount on the ends of the reads with reads with clipping present" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
      readOneFivePrime=5, readOneThreePrime=2, readTwoFivePrime=3, readTwoThreePrime=4)
    val Seq(r1, r2) = builder.addPair(start1=104, start2=150)

    r1.cigar = "4H46M"
    r2.cigar = "44M6H"

    val prior = StartsAndEnds(r1, r2)
    r1.end shouldBe (r2.start - 1)
    clipper.clipPair(r1, r2)
    // R1 5' end has one more base clipped
    // R1 3' end has two bases clipped
    // R2 5' end has has no more bases clipped
    // R2 e' end has four bases clipped
    prior.checkClipping(r1, r2, 1, 2, 0, 4)
  }

  it should "clip a fixed amount on the ends of the reads then clip overlapping reads" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
      readOneFivePrime=0, readOneThreePrime=1, readTwoFivePrime=0, readTwoThreePrime=1, clipOverlappingReads=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=146)

    val prior = StartsAndEnds(r1, r2)
    r1.end shouldBe (r2.start + 3) // four bases overlap!
    clipper.clipPair(r1, r2)
    r1.end shouldBe (r2.start - 1)
    prior.checkClipping(r1, r2, 0, 2, 0, 2) // clipping due to overlapping reads
  }

  Seq(Plus, Minus).foreach { strand1 =>
    Seq(Plus, Minus).foreach { strand2 =>
      it should s"clip a fixed amount on the ends of the reads in $strand1/$strand2" in {
        val builder = new SamBuilder(readLength=50)
        val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
          readOneFivePrime=1, readOneThreePrime=2, readTwoFivePrime=3, readTwoThreePrime=4)
        val Seq(r1, r2) = builder.addPair(start1=100, start2=150, strand1=strand1, strand2=strand2)

        val prior = StartsAndEnds(r1, r2)
        r1.end shouldBe (r2.start - 1)
        clipper.clipPair(r1, r2)
        prior.checkClipping(r1, r2, 1, 2, 3, 4)
      }
    }
  }

  "ClippingMetrics.add" should "add two metrics" in {
    val fragment = ClippingMetrics(ReadType.Fragment, 1, 2,  3,  4,  5,  6, 7,  8,   9, 10, 11, 12, 13, 14, 15)
    val readOne  = ClippingMetrics(ReadType.ReadTwo,  2, 3,  4,  5,  6,  7, 8,  9,  10, 11, 12, 13, 14, 15, 16)
    val readTwo  = ClippingMetrics(ReadType.ReadTwo,  3, 4,  5,  6,  7,  8, 9,  10, 11, 12, 13, 14, 15, 16, 17)
    val all      = ClippingMetrics(ReadType.All,      6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48)

    // Check that adding fragment, readOne, and readTwo is correct
    val added    = ClippingMetrics(ReadType.All)
    added.add(fragment, readOne, readTwo)
    added.productIterator.toSeq should contain theSameElementsInOrderAs all.productIterator.toSeq

    // Check that two different metrics have all different values
    added.productIterator.zip(fragment.productIterator).count { case (left, right) => left != right } shouldBe added.productIterator.length
  }

  "ClipBam" should "clip overlapping reads, update mate info, and reset NM, UQ & MD" in {
    val random = new Random(1)
    val builder = new SamBuilder(readLength=50, sort=Some(Queryname))
    builder.addPair(name="q1", start1=100, start2=140) // overlaps by ten
    builder.addPair(name="q2", start1=200, start2=242) // overlaps by eight
    builder.addPair(name="q3", start1=300, start2=344) // overlaps by six
    builder.addPair(name="q4", start1=400, start2=446) // overlaps by four
    builder.addPair(name="q5", start1=500, start2=548) // overlaps by two
    builder.addPair(name="q6", start1=600, start2=650) // does not overlap

    // Set appropriate bases and tags on each read
    val refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
    val chr1    = refFile.nextSequence().getBases
    builder.iterator.foreach { r =>
      val bases = ("A" * 50).getBytes
      val n     = random.nextInt(50)
      bases(n)  = 'C'
      r.bases = bases
      SequenceUtil.calculateMdAndNmTags(r.asSam, chr1, true, true)
      r(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(r.asSam, chr1, 0)
    }

    val out        = makeTempFile("out.", ".bam")
    val metricsOut = makeTempFile("out.", ".txt")
    new ClipBam(
      input=builder.toTempFile(), output=out, metrics=Some(metricsOut), ref=ref,
      readOneFivePrime=2, readOneThreePrime=0, readTwoFivePrime=3, readTwoThreePrime=0,
      clipOverlappingReads=true
    ).execute()
    val clipped = SamSource(out).toSeq

    // Check that some reads got clipped due to being overlapped
    clipped.filter(_.negativeStrand).count(_.cigar.head.operator.toString == "H") shouldBe 5
    clipped.filterNot(_.negativeStrand).count(_.cigar.last.operator.toString == "H") shouldBe 5

    // Check that all reads got clipped due to the fixed clipping
    clipped.filter(_.negativeStrand).count(_.cigar.last.toString == "3H") shouldBe 6
    clipped.filterNot(_.negativeStrand).count(_.cigar.head.toString == "2H") shouldBe 6

    // Check that everything is in coordinate order
    clipped.sliding(2).foreach { case Seq(lhs, rhs) => lhs.start <= rhs.start shouldBe true}

    // Validate that mate information didn't get mangled
    clipped.groupBy(_.name).values.foreach { case Seq(lhs, rhs) =>
      for ((a,b) <- Seq((lhs, rhs), (rhs, lhs))) {
        a.mateStart shouldBe b.start
        a.mateNegativeStrand shouldBe b.negativeStrand
        a[String]("MC") shouldBe b.cigar.toString()
      }
    }

    // Validate that the NM/UQ/MD tags got re-calculated
    clipped.foreach { r =>
      val Seq(nm, uq, md) = Seq("NM", "UQ", "MD").map(t => r[AnyRef](t))
      SequenceUtil.calculateMdAndNmTags(r.asSam, chr1, true, true)
      r(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(r.asSam, chr1, 0)
      val Seq(expNm, expUq, expMd) = Seq("NM", "UQ", "MD").map(t => r[AnyRef](t))

      nm shouldBe expNm
      uq shouldBe expUq
      md shouldBe expMd
    }

    val metrics = Metric.read[ClippingMetrics](metricsOut)
    val pair = ClippingMetrics(read_type=ReadType.Pair)
    pair.add(metrics.filter(m => m.read_type == ReadType.ReadOne || m.read_type == ReadType.ReadTwo):_*)
    val all = pair.copy(read_type=ReadType.All)
    metrics.foreach {
      case metric if metric.read_type == ReadType.Fragment =>
        metric shouldBe ClippingMetrics(read_type=ReadType.Fragment)
      case metric if metric.read_type == ReadType.ReadOne =>
        metric.reads shouldBe 6
        metric.reads_unmapped shouldBe 0
        metric.reads_clipped_pre shouldBe 0
        metric.reads_clipped_post shouldBe 6
        metric.reads_clipped_five_prime shouldBe 6
        metric.reads_clipped_three_prime shouldBe 0
        metric.reads_clipped_overlapping shouldBe 5 // all but q6
        metric.bases shouldBe 273 // 50*6 - (2+4+6+8+10)/2 - (2*6)
        metric.bases_clipped_pre shouldBe 0
        metric.bases_clipped_post shouldBe 27 // (2+4+6+8+10)/2 + (2*6)
        metric.bases_clipped_five_prime shouldBe 12 // 2*6
        metric.bases_clipped_three_prime shouldBe 0
        metric.bases_clipped_overlapping shouldBe 15 // (2+4+6+8+10)/2
      case metric if metric.read_type == ReadType.ReadTwo =>
        metric.reads shouldBe 6
        metric.reads_unmapped shouldBe 0
        metric.reads_clipped_pre shouldBe 0
        metric.reads_clipped_post shouldBe 6
        metric.reads_clipped_five_prime shouldBe 6
        metric.reads_clipped_three_prime shouldBe 0
        metric.reads_clipped_overlapping shouldBe 5 // all but q6
        metric.bases shouldBe 267 // 50*6 - (2+4+6+8+10)/2 - (3*6)
        metric.bases_clipped_pre shouldBe 0
        metric.bases_clipped_post shouldBe 33 // (2+4+6+8+10)/2 + (3*6)
        metric.bases_clipped_five_prime shouldBe 18 // 3*6
        metric.bases_clipped_three_prime shouldBe 0
        metric.bases_clipped_overlapping shouldBe 15 // (2+4+6+8+10)/2
      case metric if metric.read_type == ReadType.Pair =>
        metric.productIterator.toSeq should contain theSameElementsInOrderAs pair.productIterator.toSeq
      case metric if metric.read_type == ReadType.All =>
        metric.productIterator.toSeq should contain theSameElementsInOrderAs all.productIterator.toSeq
    }
  }

  it should "clip fragment reads, and reset NM, UQ & MD" in {
    val random = new Random(1)
    val builder = new SamBuilder(readLength=50, sort=Some(Queryname))
    builder.addFrag(start=100)
    builder.addFrag(start=200, strand=Minus)
    builder.addFrag(start=300, cigar="40S10M") // should be fully clipped

    // Set appropriate bases and tags on each read
    val refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
    val chr1    = refFile.nextSequence().getBases
    builder.iterator.foreach { r =>
      val bases = ("A" * 50).getBytes
      val n     = random.nextInt(50)
      bases(n)  = 'C'
      r.bases = bases
      SequenceUtil.calculateMdAndNmTags(r.asSam, chr1, true, true)
      r(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(r.asSam, chr1, 0)
    }

    val out     = makeTempFile("out.", ".bam")
    val metrics = makeTempFile("out.", ".txt")
    new ClipBam(
      input=builder.toTempFile(), output=out, metrics=Some(metrics), ref=ref,
      readOneFivePrime=2, readOneThreePrime=10, readTwoFivePrime=5, readTwoThreePrime=5, // read two info will be ignored
      clipOverlappingReads=true // overlapping will be ignored
    ).execute()
    val clipped = SamSource(out).toSeq

    // Check that all reads got clipped due to the fixed clipping
    clipped.filter(_.mapped).filter(_.negativeStrand).count(_.cigar.last.toString == "2H") shouldBe 1
    clipped.filter(_.mapped).filterNot(_.negativeStrand).count(_.cigar.head.toString == "2H") shouldBe 1

    // Check that everything is in coordinate order
    clipped.sliding(2).foreach { case Seq(lhs, rhs) =>
      if (lhs.mapped && rhs.mapped) lhs.start <= rhs.start shouldBe true
      else if (lhs.unmapped) rhs.unmapped shouldBe true
    }

    // Validate that the NM/UQ/MD tags got re-calculated
    clipped.filter(_.mapped).foreach { r =>
      val Seq(nm, uq, md) = Seq("NM", "UQ", "MD").map(t => r[AnyRef](t))
      SequenceUtil.calculateMdAndNmTags(r.asSam, chr1, true, true)
      r(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(r.asSam, chr1, 0)
      val Seq(expNm, expUq, expMd) = Seq("NM", "UQ", "MD").map(t => r[AnyRef](t))

      nm shouldBe expNm
      uq shouldBe expUq
      md shouldBe expMd
    }

    Metric.read[ClippingMetrics](metrics).foreach {
      case metric if metric.read_type == ReadType.Fragment | metric.read_type == ReadType.All =>
        metric.reads shouldBe 3
        metric.reads_unmapped shouldBe 1
        metric.reads_clipped_pre shouldBe 1
        metric.reads_clipped_post shouldBe 3
        metric.reads_clipped_five_prime shouldBe 2
        metric.reads_clipped_three_prime shouldBe 3
        metric.reads_clipped_overlapping shouldBe 0
        metric.bases shouldBe 76 // 50*2 - 12*2
        metric.bases_clipped_pre shouldBe 40 // from frag 3
        metric.bases_clipped_post shouldBe 74 // 12*2 + 50
        metric.bases_clipped_five_prime shouldBe 4
        metric.bases_clipped_three_prime shouldBe 30 // 10*3
        metric.bases_clipped_overlapping shouldBe 0
      case metric if metric.read_type == ReadType.ReadOne | metric.read_type == ReadType.ReadTwo | metric.read_type == ReadType.Pair =>
        metric shouldBe ClippingMetrics(read_type=metric.read_type)
    }
  }

  private def clipBases(rec: SamRecord, left: Int, right: Int): String = rec.basesString.slice(left, rec.length-right)
  private def clipQuals(rec: SamRecord, left: Int, right: Int): String = rec.qualsString.slice(left, rec.length-right)
  private def maskBases(rec: SamRecord, left: Int, right: Int): String = "N"*left + clipBases(rec, left, right) + "N"*right
  private def maskQuals(rec: SamRecord, left: Int, right: Int): String = "#"*left + clipQuals(rec, left, right) + "#"*right

  Seq((Soft, SoftWithMask), (Soft, Hard), (SoftWithMask, Hard)).foreach { case (prior, mode) =>
    it should s"upgrade existing clipping from $prior to $mode with --upgrade-clipping" in {
      val builder = new SamBuilder(readLength=50, sort=Some(Queryname))
      val frag1 = builder.addFrag(name="q1", start=100).value
      val frag2 = builder.addFrag(name="q2", start=200, strand=Minus).value

      // clip them
      val softClipper = new SamRecordClipper(mode=prior, autoClipAttributes=true)
      val priorClipper = new SamRecordClipper(mode=prior, autoClipAttributes=true)
      Seq(frag1, frag2).foreach { frag =>
        // always soft-clip bases of soft at the start and end
        softClipper.clip5PrimeEndOfRead(frag, 10) shouldBe 10
        softClipper.clip3PrimeEndOfRead(frag, 4) shouldBe 4
        // now do some more clipping based on the "prior" mode
        priorClipper.clip5PrimeEndOfRead(frag, 5) shouldBe 0
        priorClipper.clip3PrimeEndOfRead(frag, 2) shouldBe 0
      }

      val out     = makeTempFile("out.", ".bam")
      val metrics = makeTempFile("out.", ".txt")
      new ClipBam(input=builder.toTempFile(), output=out, metrics=Some(metrics), ref=ref, clippingMode=mode, upgradeClipping=true).execute()
      val clipped = SamSource(out).toSeq

      clipped.length shouldBe 2

      (prior, mode) match {
        case (Soft, SoftWithMask) =>
          clipped.head.cigar.toString shouldBe "10S36M4S"
          clipped.head.basesString shouldBe maskBases(frag1, 10, 4)
          clipped.head.qualsString shouldBe maskQuals(frag1, 10, 4)
          clipped.last.cigar.toString shouldBe "4S36M10S"
          clipped.last.basesString shouldBe maskBases(frag2, 4, 10)
          clipped.last.qualsString shouldBe maskQuals(frag2, 4, 10)
        case (_, Hard) =>
          clipped.head.cigar.toString shouldBe "10H36M4H"
          clipped.head.basesString shouldBe clipBases(frag1, 10, 4)
          clipped.head.qualsString shouldBe clipQuals(frag1, 10, 4)
          clipped.last.cigar.toString shouldBe "4H36M10H"
          clipped.last.basesString shouldBe clipBases(frag2, 4, 10)
          clipped.last.qualsString shouldBe clipQuals(frag2, 4, 10)
        case _ => unreachable(s"Bug: $prior $mode")
      }
    }
  }

  Seq((Soft, Soft), (SoftWithMask, Soft), (SoftWithMask, SoftWithMask), (Hard, Hard), (Hard, SoftWithMask), (Hard, Soft)).foreach { case (prior, mode) =>
    it should s"not upgrade existing clipping from $prior to $mode with --upgrade-clipping" in {
      val builder = new SamBuilder(readLength=50, sort=Some(Queryname))
      val frag1 = builder.addFrag(name="q1", start=100).value
      val frag2 = builder.addFrag(name="q2", start=200, strand=Minus).value

      // Clip some bases based on the "prior" mode
      val priorClipper = new SamRecordClipper(mode=prior, autoClipAttributes=true)
      Seq(frag1, frag2).foreach { frag =>
        priorClipper.clip5PrimeEndOfRead(frag, 10) shouldBe 10
        priorClipper.clip3PrimeEndOfRead(frag, 4) shouldBe 4
      }

      val out     = makeTempFile("out.", ".bam")
      val metrics = makeTempFile("out.", ".txt")
      new ClipBam(input=builder.toTempFile(), output=out, metrics=Some(metrics), ref=ref, clippingMode=mode, upgradeClipping=true).execute()
      val clipped = SamSource(out).toSeq

      clipped.length shouldBe 2

      prior match {
        case Soft =>
          clipped.head.cigar.toString shouldBe "10S36M4S"
          clipped.head.basesString shouldBe frag1.basesString
          clipped.head.qualsString shouldBe frag1.qualsString
          clipped.last.cigar.toString shouldBe "4S36M10S"
          clipped.last.basesString shouldBe frag2.basesString
          clipped.last.qualsString shouldBe frag2.qualsString
        case SoftWithMask =>
          clipped.head.cigar.toString shouldBe "10S36M4S"
          clipped.head.basesString shouldBe maskBases(frag1, 10, 4)
          clipped.head.qualsString shouldBe maskQuals(frag1, 10, 4)
          clipped.last.cigar.toString shouldBe "4S36M10S"
          clipped.last.basesString shouldBe maskBases(frag2, 4, 10)
          clipped.last.qualsString shouldBe maskQuals(frag2, 4, 10)
        case Hard =>
          clipped.head.cigar.toString shouldBe "10H36M4H"
          clipped.head.basesString shouldBe frag1.basesString
          clipped.head.qualsString shouldBe frag1.qualsString
          clipped.last.cigar.toString shouldBe "4H36M10H"
          clipped.last.basesString shouldBe frag2.basesString
          clipped.last.qualsString shouldBe frag2.qualsString
      }
    }
  }

  it should "clip FR reads that extend past the mate" in {
    val builder = new SamBuilder(readLength=50, sort=Some(Queryname))
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref, clipBasesPastMate=true)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=90)

    r1.end shouldBe 149
    r2.end shouldBe 139
    clipper.clipPair(r1, r2)
    r1.start shouldBe r2.start
    r1.end shouldBe r2.end
  }
}
