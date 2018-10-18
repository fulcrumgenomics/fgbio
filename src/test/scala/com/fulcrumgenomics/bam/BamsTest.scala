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

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.SAMFileHeader.GroupOrder
import htsjdk.samtools.SamPairUtil
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools.reference.{ReferenceSequence, ReferenceSequenceFile, ReferenceSequenceFileWalker}

import scala.collection.Iterator

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

  "Bams.sorterByTag" should "should sort by a tag" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addFrag(name="q1", start=100).foreach { r => r("ZZ") = 3}
    builder.addPair(name="p1", start1=100, start2=300).foreach { r => r("ZZ") = 2}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 1}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 4}

    Bams.sortByTag[Int](iterator=builder.iterator, header=builder.header, tag="ZZ")
      .map(r => r[Int]("ZZ")).toSeq should contain theSameElementsInOrderAs Seq(1, 2, 2, 3, 4)
  }

  it should "fail if the tag is required (no missing value)" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addFrag(name="q1", start=100).foreach { r => r("ZZ") = 3}
    builder.addPair(name="p1", start1=100, start2=300) // *** NO TAGS ***
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 1}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 4}

    an[IllegalStateException] should be thrownBy Bams.sortByTag[Int](iterator=builder.iterator, header=builder.header, tag="ZZ")
  }

  it should "handle missing tags when a missing value is given" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addFrag(name="q1", start=100).foreach { r => r("ZZ") = 3}
    builder.addPair(name="p1", start1=100, start2=300) // *** NO TAGS ***
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 1}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 4}

    Bams.sortByTag[Int](iterator=builder.iterator, header=builder.header, tag="ZZ", defaultValue=Some(2))
      .map(r => r.get[Int]("ZZ").getOrElse(2)).toSeq should contain theSameElementsInOrderAs Seq(1, 2, 2, 3, 4)
  }


  it should "should sort by a tag with a transform applied to the tag" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname))
    builder.addFrag(name="q1", start=100).foreach { r => r("ZZ") = 3}
    builder.addPair(name="p1", start1=100, start2=300).foreach { r => r("ZZ") = 2}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 1}
    builder.addFrag(name="q2", start=200).foreach { r => r("ZZ") = 4}

    val f: Int => Float = a => -(a/2.0f)

    Bams.sortByTransformedTag[Int,Float](iterator=builder.iterator, header=builder.header, tag="ZZ", transform=f)
      .map(r => r[Int]("ZZ")).toSeq should contain theSameElementsInOrderAs Seq(4, 3, 2, 2, 1)
  }

  "Template.primaryReads" should "return an iterator over just the primary reads" in {
    val builder = new SamBuilder()
    builder.addPair("q1", contig=0, start1=100, start2=200)
    builder.addPair("q1", contig=1, start1=2000, start2=3000).foreach(_.supplementary = true)
    builder.addPair("q1", contig=2, start1=2000, start2=3000).foreach(_.secondary = true)
    val t1 = Template(builder.iterator)
    t1.allSupplementaryAndSecondary should have size 4
    t1.primaryReads should have size 2
    t1.primaryReads.toSeq should contain theSameElementsInOrderAs builder.iterator.filter(_.refIndex == 0).toSeq
  }

  "Template.pairOrientation" should "return None if the template is not a pair" in {
    val recs = new SamBuilder().addFrag(contig=0, start=100)
    val template = Template(recs.iterator)
    template.pairOrientation shouldBe None
  }

  it should "return None if either read is unmapped or the reads are on different chromosomes" in {
    val recs = new SamBuilder().addPair(contig=1, start1=1000, start2=2000)
    val Seq(r1, r2) = recs
    val template = Template(recs.iterator)
    template.pairOrientation.isDefined shouldBe true

    r1.unmapped = true
    template.pairOrientation shouldBe None
    r2.unmapped = true
    template.pairOrientation shouldBe None
    r1.mapped = true
    r2.mapped = true
    template.pairOrientation.isDefined shouldBe true

    r2.refIndex = 2
    template.pairOrientation shouldBe None
  }

  it should "return the correct pair orientation for pairs mapped on the same chrom" in {
    val Seq(r1, r2) = new SamBuilder().addPair(contig=1, start1=1000, start2=2000, strand1=Plus, strand2=Minus)
    val template = Template(Iterator(r1, r2))
    template.pairOrientation shouldBe Some(PairOrientation.FR)

    r1.negativeStrand = true
    r2.positiveStrand = true
    SamPairUtil.setMateInfo(r1.asSam, r2.asSam)
    template.pairOrientation shouldBe Some(PairOrientation.RF)

    r1.positiveStrand = true
    r2.positiveStrand = true
    SamPairUtil.setMateInfo(r1.asSam, r2.asSam)
    template.pairOrientation shouldBe Some(PairOrientation.TANDEM)

    r1.negativeStrand = true
    r2.negativeStrand = true
    SamPairUtil.setMateInfo(r1.asSam, r2.asSam)
    template.pairOrientation shouldBe Some(PairOrientation.TANDEM)
  }

  "Template.unmapped" should "unmap all primary reads and reset pair information" in {
    val recs = new SamBuilder().addPair(contig=1, start1=1000, start2=2000, strand1=Plus, strand2=Minus)
    val template = Template(recs.iterator)

    val unmapped = template.unmapped
    val Seq(r1, r2) = unmapped.primaryReads.toSeq
    r1.unmapped     shouldBe true
    r1.mateUnmapped shouldBe true
    r2.unmapped     shouldBe true
    r2.mateUnmapped shouldBe true

    r1.refIndex     shouldBe SamRecord.UnmappedReferenceIndex
    r1.mateRefIndex shouldBe SamRecord.UnmappedReferenceIndex
    r2.refIndex     shouldBe SamRecord.UnmappedReferenceIndex
    r2.mateRefIndex shouldBe SamRecord.UnmappedReferenceIndex

    r1.start     shouldBe SamRecord.UnmappedStart
    r1.mateStart shouldBe SamRecord.UnmappedStart
    r2.start     shouldBe SamRecord.UnmappedStart
    r2.mateStart shouldBe SamRecord.UnmappedStart

    r1.cigar shouldBe Cigar.empty
    r2.cigar shouldBe Cigar.empty

    r1.mapq shouldBe SamRecord.ZeroMappingQuality
    r2.mapq shouldBe SamRecord.ZeroMappingQuality

    r2.basesString shouldBe Sequences.revcomp(template.r2.get.basesString)
    r2.quals.reverseIterator.toSeq should contain theSameElementsInOrderAs template.r2.get.quals.toSeq
  }

  it should "discard supplementary and secondary alignment records" in {
    val builder = new SamBuilder()
    builder.addPair("q1", contig=0, start1=100, start2=200)
    builder.addPair("q1", contig=1, start1=2000, start2=3000).foreach(_.supplementary = true)
    builder.addPair("q1", contig=2, start1=2000, start2=3000).foreach(_.secondary = true)
    val template = Template(builder.iterator)

    val unmapped = template.unmapped
    unmapped.allReads should have size 2
    unmapped.r1Secondaries should have size 0
    unmapped.r2Secondaries should have size 0
    unmapped.r1Supplementals should have size 0
    unmapped.r2Supplementals should have size 0

    unmapped.allReads.foreach { r =>
      r.secondary shouldBe false
      r.supplementary shouldBe false
    }

  }
}
