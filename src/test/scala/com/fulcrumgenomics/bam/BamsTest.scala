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

import java.util

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.reference.{ReferenceSequence, ReferenceSequenceFile, ReferenceSequenceFileWalker}

class BamsTest extends UnitSpec {
  /* Dummy implementation of a reference sequence file walker that always returns chr1 w/2000 As. */
  private object DummyRefWalker extends ReferenceSequenceFileWalker(null.asInstanceOf[ReferenceSequenceFile]) {
    private val bases = new Array[Byte](2000)
    util.Arrays.fill(bases, 'A'.toByte)
    private val ref = new ReferenceSequence("chr1", 0, bases)
    override def get(sequenceIndex: Int): ReferenceSequence = ref
  }

  "Bams.queryGroupedIterator" should "work on a reader that is already query sorted" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.queryname)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toReader).map(_.getReadName).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  it should "not need to sort a reader that is query grouped" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.unsorted)
    builder.header.setGroupOrder(GroupOrder.query)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toReader).map(_.getReadName).toSeq shouldBe Seq("q1", "p1", "p1", "q2")
  }

  it should "sort a coordinate sorted input" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toReader).map(_.getReadName).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  it should "sort an unsorted input" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.unsorted)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toReader).map(_.getReadName).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  "Bams.templateIterator" should "return template objects in order" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.queryname)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="f1", start=100)
    builder.addPair(name="p2", start1=500, start2=200)
    builder.addPair(name="p0", start1=700, start2=999)

    val templates = Bams.templateIterator(builder.toReader).toSeq
    templates should have size 4
    templates.map(_.name) shouldBe Seq("f1", "p0", "p1", "p2")

    templates.foreach {t =>
      (t.r1Supplementals ++ t.r1Secondaries ++ t.r2Supplementals ++ t.r2Secondaries) shouldBe 'empty
      if (t.name startsWith "f") {
        t.r1.exists(r => !r.getReadPairedFlag) shouldBe true
        t.r2 shouldBe 'empty
      }
      else {
        t.r1.exists(r => r.getFirstOfPairFlag) shouldBe true
        t.r2.exists(r => r.getSecondOfPairFlag) shouldBe true
      }
    }
  }

  it should "handle reads with secondary and supplementary records" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
    Range.inclusive(1, 5).foreach { i=>
      builder.addPair(s"p$i", start1=i*100,  start2=i*100  + 500)
      builder.addPair(s"p$i", start1=i*1000, start2=i*1000 + 500).foreach(_.setNotPrimaryAlignmentFlag(true))
      builder.addPair(s"p$i", start1=i*10,   start2=i*10   + 500).foreach(_.setSupplementaryAlignmentFlag(true))
    }

    val templates = Bams.templateIterator(builder.toReader).toSeq
    templates.map(_.name) shouldBe Seq("p1", "p2", "p3", "p4", "p5")
    templates.foreach { t =>
      t.r1.exists(r => r.getFirstOfPairFlag  && !r.isSecondaryOrSupplementary) shouldBe true
      t.r2.exists(r => r.getSecondOfPairFlag && !r.isSecondaryOrSupplementary) shouldBe true
      Seq(t.r1Secondaries, t.r2Secondaries, t.r1Supplementals, t.r2Supplementals).foreach(rs => rs.size shouldBe 1)

      t.r1Secondaries.foreach  (r => (r.getFirstOfPairFlag  && r.getNotPrimaryAlignmentFlag) shouldBe true)
      t.r2Secondaries.foreach  (r => (r.getSecondOfPairFlag && r.getNotPrimaryAlignmentFlag) shouldBe true)
      t.r1Supplementals.foreach(r => (r.getFirstOfPairFlag  && r.getSupplementaryAlignmentFlag) shouldBe true)
      t.r2Supplementals.foreach(r => (r.getSecondOfPairFlag && r.getSupplementaryAlignmentFlag) shouldBe true)
    }
  }

  "Bams.regenerateNmUqMdTags" should "null out fields on unmapped reads" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate)
    val rec = builder.addFrag(unmapped=true).get
    rec.setAttribute("NM", 7)
    rec.setAttribute("MD", "6A7C8T9G")
    rec.setAttribute("UQ", 237)
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.getAttribute("NM") shouldBe null
    rec.getAttribute("MD") shouldBe null
    rec.getAttribute("UQ") shouldBe null
  }

  it should "regenerate tags on mapped reads" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate, readLength=10, baseQuality=20)
    val rec = builder.addFrag(start=100).get
    rec.setAttribute("NM", 7)
    rec.setAttribute("MD", "6A7C8T9G")
    rec.setAttribute("UQ", 237)
    rec.setReadString("AAACAAAATA")
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.getAttribute("NM") shouldBe 2
    rec.getAttribute("MD") shouldBe "3A4A1"
    rec.getAttribute("UQ") shouldBe 40
  }

  it should "only regenerate tags that are already present on mapped reads" in {
    val builder = new SamRecordSetBuilder(sortOrder=SortOrder.coordinate, readLength=10, baseQuality=20)
    val rec = builder.addFrag(start=100).get
    rec.setReadString("AAACAAAATA")

    rec.setAttribute("NM", 7)
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.getAttribute("NM") shouldBe 2
    rec.getAttribute("MD") shouldBe null
    rec.getAttribute("UQ") shouldBe null

    rec.setAttribute("NM", null)
    rec.setAttribute("MD", "6A7C8T9G")
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.getAttribute("NM") shouldBe null
    rec.getAttribute("MD") shouldBe "3A4A1"
    rec.getAttribute("UQ") shouldBe null

    rec.setAttribute("MD", null)
    rec.setAttribute("UQ", 237)
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.getAttribute("NM") shouldBe null
    rec.getAttribute("MD") shouldBe null
    rec.getAttribute("UQ") shouldBe 40
  }
}
