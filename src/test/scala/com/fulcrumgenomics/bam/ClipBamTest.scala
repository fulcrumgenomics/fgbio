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

import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.fasta.ReferenceSetBuilder
import com.fulcrumgenomics.testing.SamBuilder._
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.SequenceUtil
import htsjdk.samtools.SAMTag

import scala.util.Random

class ClipBamTest extends UnitSpec {
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

  "ClipOverlappingReads.clip" should "not clip reads where either read is unaligned" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, unmapped2=true)
    val expected = r1.cigar
    clipper.clip(r1, r2)
    r1.cigar shouldBe expected
  }

  it should "not clip reads that are on different chromosomes" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)
    r2.refIndex = 1

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clip(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "not clip reads that are abutting but not overlapped" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=150)

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clip(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "not clip non-FR reads " in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, strand2=Plus)

    val (r1Exp, r2Exp) = (r1.cigar, r2.cigar)
    clipper.clip(r1, r2)
    r1.cigar shouldBe r1Exp
    r2.cigar shouldBe r2Exp
  }

  it should "clip reads that are fully overlapped" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)

    clipper.clip(r1, r2)
    r1.cigar.toString shouldBe "25M25H"
    r2.cigar.toString shouldBe "25H25M"
  }

  it should "clip reads that are overlapped by just one base" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=149)

    r1.end shouldBe r2.start
    clipper.clip(r1, r2)
    r1.end shouldBe (r2.start - 1)
  }

  it should "handle reads that contain insertions" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2I8M", cigar2="10M2I38M")

    r1.end >= r2.start shouldBe true
    clipper.clip(r1, r2)
    r1.end >= r2.start shouldBe false
  }

  it should "handle reads that contain deletions" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2D10M", cigar2="10M2D40M")

    r1.end >= r2.start shouldBe true
    clipper.clip(r1, r2)
    r1.end >= r2.start shouldBe false
  }

  // This is a weird test that ensures that things terminate correctly when due to highly clipped
  // input reads, one of the reads becomes unmapped during clipping
  it should "handle reads that get unmapped because they are fully overlapped" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1a, r2a) = builder.addPair(name="q1", start1=100, start2=100, cigar1="3M47S", cigar2="48S2M")
    val Seq(r1b, r2b) = builder.addPair(name="q2", start1=100, start2=100, cigar1="2M48S", cigar2="47S3M")

    for ((r1, r2) <- Seq((r1a, r2a), (r1b, r2b))) {
      clipper.clip(r1, r2)
      val ok = r1.unmapped || r2.unmapped || r1.end < r2.start
      ok shouldBe true
    }
  }

  it should "clip a fixed amount on the ends of the reads with reads that do not overlap" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
      readOneFivePrime=1, readOneThreePrime=2, readTwoFivePrime=3, readTwoThreePrime=4)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=150)

    val prior = StartsAndEnds(r1, r2)
    r1.end shouldBe (r2.start - 1)
    clipper.clip(r1, r2)
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
    clipper.clip(r1, r2)
    // R1 5' end has one more base clipped
    // R1 3' end has two bases clipped
    // R2 5' end has has no more bases clipped
    // R2 e' end has four bases clipped
    prior.checkClipping(r1, r2, 1, 2, 0, 4)
  }

  it should "clip a fixed amount on the ends of the reads then clip overlapping reads" in {
    val builder = new SamBuilder(readLength=50)
    val clipper = new ClipBam(input=dummyBam, output=dummyBam, ref=ref,
      readOneFivePrime=0, readOneThreePrime=1, readTwoFivePrime=0, readTwoThreePrime=1)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=146)

    val prior = StartsAndEnds(r1, r2)
    r1.end shouldBe (r2.start + 3) // four bases overlap!
    clipper.clip(r1, r2)
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
        clipper.clip(r1, r2)
        prior.checkClipping(r1, r2, 1, 2, 3, 4)
      }
    }
  }

  "ClipBam" should "clip overlapping reads, update mate info, and reset NM, UQ & MD" in {
    val random = new Random(1)
    val builder = new SamBuilder(readLength=50)
    builder.addPair(name="q1", start1=100, start2=140)
    builder.addPair(name="q2", start1=200, start2=242)
    builder.addPair(name="q3", start1=300, start2=344)
    builder.addPair(name="q4", start1=400, start2=446)
    builder.addPair(name="q5", start1=500, start2=548)
    builder.addPair(name="q6", start1=600, start2=650)

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

    val out = makeTempFile("out.", ".bam")
    new ClipBam(
      input=builder.toTempFile(), output=out, ref=ref,
      readOneFivePrime=2, readOneThreePrime=0, readTwoFivePrime=2, readTwoThreePrime=0
    ).execute()
    val clipped = SamSource(out).toSeq

    // Check that some reads got clipped due to being overlapped
    clipped.filter(_.negativeStrand).count(_.cigar.head.operator.toString == "H") shouldBe 5
    clipped.filterNot(_.negativeStrand).count(_.cigar.last.operator.toString == "H") shouldBe 5

    // Check that all reads got clipped due to the fixed clipping
    clipped.filter(_.negativeStrand).count(_.cigar.last.toString == "2H") shouldBe 6
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
  }
}
