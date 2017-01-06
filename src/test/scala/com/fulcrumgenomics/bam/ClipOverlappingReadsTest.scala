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

import com.fulcrumgenomics.fasta.ReferenceSetBuilder
import com.fulcrumgenomics.testing.SamRecordSetBuilder._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import htsjdk.samtools.{SAMTag, SamReaderFactory}
import htsjdk.samtools.SAMTag.{MD, NM}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

import htsjdk.samtools.util.SequenceUtil

import scala.collection.JavaConversions.iterableAsScalaIterable
import scala.util.Random

class ClipOverlappingReadsTest extends UnitSpec {
  private val ref = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 500) // 5000 bases
    builder.add("chr2").add("CCCCCCCCCC", 500) // 5000 bases
    builder.toTempFile()
  }

  private val dummyBam = makeTempFile("dummy.", ".bam")

  "ClipOverlappingReads.clip" should "not clip reads where either read is unaligned" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, record2Unmapped=true)
    val expected = r1.getCigarString
    clipper.clip(r1, r2)
    r1.getCigarString shouldBe expected
  }

  it should "not clip reads that are on different chromosomes" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)
    r2.setReferenceIndex(1)

    val (r1Exp, r2Exp) = (r1.getCigarString, r2.getCigarString)
    clipper.clip(r1, r2)
    r1.getCigarString shouldBe r1Exp
    r2.getCigarString shouldBe r2Exp
  }

  it should "not clip reads that are abutting but not overlapped" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=150)

    val (r1Exp, r2Exp) = (r1.getCigarString, r2.getCigarString)
    clipper.clip(r1, r2)
    r1.getCigarString shouldBe r1Exp
    r2.getCigarString shouldBe r2Exp
  }

  it should "not clip non-FR reads " in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100, strand2=Plus)

    val (r1Exp, r2Exp) = (r1.getCigarString, r2.getCigarString)
    clipper.clip(r1, r2)
    r1.getCigarString shouldBe r1Exp
    r2.getCigarString shouldBe r2Exp
  }

  it should "clip reads that are fully overlapped" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=100)

    clipper.clip(r1, r2)
    r1.getCigarString shouldBe "25M25H"
    r2.getCigarString shouldBe "25H25M"
  }

  it should "clip reads that are overlapped by just one base" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=149)

    r1.getAlignmentEnd shouldBe r2.getAlignmentStart
    clipper.clip(r1, r2)
    r1.getAlignmentEnd shouldBe (r2.getAlignmentStart - 1)
  }

  it should "handle reads that contain insertions" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2I8M", cigar2="10M2I38M")

    r1.getAlignmentEnd >= r2.getAlignmentStart shouldBe true
    clipper.clip(r1, r2)
    r1.getAlignmentEnd >= r2.getAlignmentStart shouldBe false
  }

  it should "handle reads that contain deletions" in {
    val builder = new SamRecordSetBuilder(readLength=50)
    val clipper = new ClipOverlappingReads(input=dummyBam, output=dummyBam, ref=ref)
    val Seq(r1, r2) = builder.addPair(start1=100, start2=130, cigar1="40M2D10M", cigar2="10M2D10M")

    r1.getAlignmentEnd >= r2.getAlignmentStart shouldBe true
    clipper.clip(r1, r2)
    r1.getAlignmentEnd >= r2.getAlignmentStart shouldBe false
  }

  "ClipOverlappingReads" should "clip overlapping reads, update mate info, and reset NM, UQ & MD" in {
    val random = new Random(1)
    val builder = new SamRecordSetBuilder(readLength=50)
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
      r.setReadBases(bases)
      SequenceUtil.calculateMdAndNmTags(r, chr1, true, true)
      r.setAttribute(SAMTag.UQ.name, SequenceUtil.sumQualitiesOfMismatches(r, chr1, 0))
    }

    val out = makeTempFile("out.", ".bam")
    new ClipOverlappingReads(input=builder.toTempFile(), output=out, ref=ref).execute()
    val clipped = SamReaderFactory.make().open(out).toSeq

    // Check that some reads got clipped
    clipped.count(_.getCigarString.contains("H")) shouldBe 10

    // Check that everything is in coordinate order
    clipped.sliding(2).foreach { case Seq(lhs, rhs) => lhs.getAlignmentStart <= rhs.getAlignmentStart shouldBe true}

    // Validate that mate information didn't get mangled
    clipped.groupBy(_.getReadName).values.foreach { case Seq(lhs, rhs) =>
      for ((a,b) <- Seq((lhs, rhs), (rhs, lhs))) {
        a.getMateAlignmentStart shouldBe b.getAlignmentStart
        a.getMateNegativeStrandFlag shouldBe b.getReadNegativeStrandFlag
        a.getAttribute("MC") shouldBe b.getCigarString
      }
    }

    // Validate that the NM/UQ/MD tags got re-calculated
    clipped.foreach { r =>
      val Seq(nm, uq, md) = Seq("NM", "UQ", "MD").map(r.getAttribute)
      SequenceUtil.calculateMdAndNmTags(r, chr1, true, true)
      r.setAttribute(SAMTag.UQ.name, SequenceUtil.sumQualitiesOfMismatches(r, chr1, 0))
      val Seq(expNm, expUq, expMd) = Seq("NM", "UQ", "MD").map(r.getAttribute)

      nm shouldBe expNm
      uq shouldBe expUq
      md shouldBe expMd
    }
  }
}
