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

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.GroupOrder
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
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toSource).map(_.name).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  it should "not need to sort a reader that is query grouped" in {
    val builder = new SamBuilder(sort=None)
    builder.header.setGroupOrder(GroupOrder.query)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toSource).map(_.name).toSeq shouldBe Seq("q1", "p1", "p1", "q2")
  }

  it should "sort a coordinate sorted input" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate))
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toSource).map(_.name).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  it should "sort an unsorted input" in {
    val builder = new SamBuilder(sort=None)
    builder.addFrag(name="q1", start=100)
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="q2", start=200)

    Bams.queryGroupedIterator(in=builder.toSource).map(_.name).toSeq shouldBe Seq("p1", "p1", "q1", "q2")
  }

  "Bams.templateIterator" should "return template objects in order" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addPair(name="p1", start1=100, start2=300)
    builder.addFrag(name="f1", start=100)
    builder.addPair(name="p2", start1=500, start2=200)
    builder.addPair(name="p0", start1=700, start2=999)

    val templates = Bams.templateIterator(builder.toSource).toSeq
    templates should have size 4
    templates.map(_.name) shouldBe Seq("f1", "p0", "p1", "p2")

    templates.foreach {t =>
      (t.r1Supplementals ++ t.r1Secondaries ++ t.r2Supplementals ++ t.r2Secondaries) shouldBe 'empty
      if (t.name startsWith "f") {
        t.r1.exists(r => !r.paired) shouldBe true
        t.r2 shouldBe 'empty
      }
      else {
        t.r1.exists(r => r.firstOfPair) shouldBe true
        t.r2.exists(r => r.secondOfPair) shouldBe true
      }
    }
  }

  it should "handle reads with secondary and supplementary records" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate))
    Range.inclusive(1, 5).foreach { i=>
      builder.addPair(s"p$i", start1=i*100,  start2=i*100  + 500)
      builder.addPair(s"p$i", start1=i*1000, start2=i*1000 + 500).foreach(_.secondary = true)
      builder.addPair(s"p$i", start1=i*10,   start2=i*10   + 500).foreach(_.supplementary = true)
    }

    val templates = Bams.templateIterator(builder.toSource).toSeq
    templates.map(_.name) shouldBe Seq("p1", "p2", "p3", "p4", "p5")
    templates.foreach { t =>
      t.r1.exists(r => r.firstOfPair  && !r.secondary & !r.supplementary) shouldBe true
      t.r2.exists(r => r.secondOfPair && !r.secondary & !r.supplementary) shouldBe true
      Seq(t.r1Secondaries, t.r2Secondaries, t.r1Supplementals, t.r2Supplementals).foreach(rs => rs.size shouldBe 1)

      t.r1Secondaries.foreach  (r => (r.firstOfPair  && r.secondary) shouldBe true)
      t.r2Secondaries.foreach  (r => (r.secondOfPair && r.secondary) shouldBe true)
      t.r1Supplementals.foreach(r => (r.firstOfPair  && r.supplementary) shouldBe true)
      t.r2Supplementals.foreach(r => (r.secondOfPair && r.supplementary) shouldBe true)
    }
  }

  "Bams.regenerateNmUqMdTags" should "null out fields on unmapped reads" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate))
    val rec = builder.addFrag(unmapped=true).get
    rec("NM") = 7
    rec("MD") = "6A7C8T9G"
    rec("UQ") = 237
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec.get[Int]("NM")    shouldBe None
    rec.get[String]("MD") shouldBe None
    rec.get[Int]("UQ")    shouldBe None
  }

  it should "regenerate tags on mapped reads" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    val rec = builder.addFrag(start=100).get
    rec("NM") = 7
    rec("MD") = "6A7C8T9G"
    rec("UQ") = 237
    rec.bases = "AAACAAAATA"
    Bams.regenerateNmUqMdTags(rec, DummyRefWalker)
    rec[Int]("NM")    shouldBe 2
    rec[String]("MD") shouldBe "3A4A1"
    rec[Int]("UQ")    shouldBe 40
  }

  "Bams.insertCoordinates" should "fail on fragments and inappropriate pairs" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    builder.addFrag(start=100).foreach(r => an[Exception] shouldBe thrownBy { Bams.insertCoordinates(r) })
    builder.addPair(start1=100, start2=100, unmapped2=true).foreach(r => an[Exception] shouldBe thrownBy { Bams.insertCoordinates(r) })
    builder.addPair(start1=100, start2=200).foreach { rec =>
      rec.mateRefIndex = rec.refIndex + 1
      an[Exception] shouldBe thrownBy { Bams.insertCoordinates(rec) }
    }
  }

  it should "calculate insert coordinates correctly" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10, baseQuality=20)
    val recs = builder.addPair(start1=100, start2=191).foreach { r => Bams.insertCoordinates(r) shouldBe (100, 200) }
  }


  "Bams.positionFromOtherEndOfTemplate" should "return None for anything that's not an FR mapped pair" in {
    val builder = new SamBuilder()
    builder.addFrag(start=100).foreach(r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
    builder.addPair(start1=100, start2=200, unmapped2=true).foreach(r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
    builder.addPair(start1=100, start2=200, strand1=Plus,  strand2=Plus).foreach (r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Minus).foreach(r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
    builder.addPair(start1=100, start2=200, strand1=Plus,  strand2=Plus).foreach (r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
    builder.addPair(start1=100, start2=200, strand1=Minus, strand2=Plus).foreach (r => Bams.positionFromOtherEndOfTemplate(r, 10) shouldBe None)
  }

  it should "correctly calculate the position from the other end of the template for FR pairs" in {
    val builder = new SamBuilder(readLength=50)
    val Seq(r1, r2) = builder.addPair(start1=101, start2=151) // Insert = (101..200), Length = 100

    Bams.positionFromOtherEndOfTemplate(r1, 101) shouldBe Some(100)
    Bams.positionFromOtherEndOfTemplate(r1, 111) shouldBe Some( 90)
    Bams.positionFromOtherEndOfTemplate(r1, 151) shouldBe Some( 50)
    Bams.positionFromOtherEndOfTemplate(r1, 200) shouldBe Some(  1)

    Bams.positionFromOtherEndOfTemplate(r2, 200) shouldBe Some(100)
    Bams.positionFromOtherEndOfTemplate(r2, 190) shouldBe Some( 90)
    Bams.positionFromOtherEndOfTemplate(r2, 150) shouldBe Some( 50)
    Bams.positionFromOtherEndOfTemplate(r2, 101) shouldBe Some(  1)
  }
}
