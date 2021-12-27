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

import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

class PileupTest extends UnitSpec {
  // A reference sequence for use below
  private val ref = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1").add("A", 1000)
    builder.add("chr2").add("C", 1000)
    val path = builder.toTempFile()
    ReferenceSequenceFileFactory.getReferenceSequenceFile(path)
  }

  private val dict = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    ref.getSequenceDictionary.fromSam
  }

  "PileupBuilder.build" should "build a simple pileup from fragment reads" in {
    val piler = new PileupBuilder(dict, mappedPairsOnly=false)
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    1 to 5 foreach { _ => builder.addFrag(start=100).foreach(_.bases = "A" * 50) }
    1 to 4 foreach { _ => builder.addFrag(start=100).foreach(_.bases = "C" * 50) }
    1 to 3 foreach { _ => builder.addFrag(start=100).foreach(_.bases = "G" * 50) }
    1 to 2 foreach { _ => builder.addFrag(start=100).foreach(_.bases = "T" * 50) }
    1 to 1 foreach { _ => builder.addFrag(start=100).foreach(_.bases = "N" * 50) }

    piler.build(builder.iterator, "chr1", 50).pile shouldBe empty
    val pile = piler.build(builder.iterator, "chr1", 100).withoutIndels
    pile.depth shouldBe 5 + 4 + 3 + 2 + 1
    pile.count(_.base == 'A') shouldBe 5
    pile.count(_.base == 'C') shouldBe 4
    pile.count(_.base == 'G') shouldBe 3
    pile.count(_.base == 'T') shouldBe 2
    pile.count(_.base == 'N') shouldBe 1
  }

  it should "build pileups from reads that contain indels" in {
    val piler = new PileupBuilder(dict, mappedPairsOnly=false)
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    builder.addFrag(name="q1", start=101, cigar="10M2D40M").foreach(_.bases = "A" * 50)
    builder.addFrag(name="q2", start=101, cigar="10M2I38M").foreach(_.bases = "C" * 50)
    builder.addFrag(name="q3", start=101, cigar="31M9I10M").foreach(_.bases = "G" * 50)
    builder.addFrag(name="q4", start=101, cigar="30M9D20M").foreach(_.bases = "T" * 50)
    builder.addFrag(name="q5", start=141, cigar="10I40M"  ).foreach(_.bases = "N" * 50)

    // Before any of the indels, we should have 4 _bases_
    piler.build(builder.iterator, "chr1", 105).baseIterator.size shouldBe 4

    // The locus before the first insertion and deletion (should include the insertion)
    val p1 = piler.build(builder.iterator, "chr1", 110)
    p1.depth shouldBe 4
    p1.iterator.size shouldBe 5
    p1.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACGT"
    p1.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q2"
    p1.baseIterator.toSeq should contain theSameElementsAs p1.withoutIndels.iterator.toSeq

    // The locus of the first deleted base
    val p2 = piler.build(builder.iterator, "chr1", 111)
    p2.depth shouldBe 4
    p2.iterator.size shouldBe 4
    p2.baseIterator.size shouldBe 3
    p2.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "CGT"
    p2.iterator.collect{ case x: DeletionEntry => x }.map(_.rec.name).next() shouldBe "q1"
    p2.baseIterator.toSeq should contain theSameElementsAs p2.withoutIndels.iterator.toSeq

    // The locus of the bigger indels
    val p3 = piler.build(builder.iterator, "chr1", 131)
    p3.depth shouldBe 4
    p3.iterator.size shouldBe 5
    p3.baseIterator.size shouldBe 3
    p3.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACG"
    p3.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q3"
    p3.iterator.collect{ case x: DeletionEntry  => x }.map(_.rec.name).next() shouldBe "q4"
    p3.baseIterator.toSeq should contain theSameElementsAs p3.withoutIndels.iterator.toSeq

    // Locus where there is a read with a leading indel
    val p4 = piler.build(builder.iterator, "chr1", 140)
    p4.depth shouldBe 4  // shouldn't count the leading indel
    p4.iterator.size shouldBe 5
    p4.baseIterator.size shouldBe 4
    p4.baseIterator.map(_.base.toChar).toSeq.sorted.mkString shouldBe "ACGT"
    p4.iterator.collect{ case x: InsertionEntry => x }.map(_.rec.name).next() shouldBe "q5"
    p4.baseIterator.toSeq should contain theSameElementsAs p4.withoutIndels.iterator.toSeq
  }

  it should "filter out reads below the min mapq" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    builder.addFrag(name="q1", start=101, mapq=9)
    builder.addFrag(name="q2", start=101, mapq=10)
    builder.addFrag(name="q3", start=101, mapq=11)
    val pile = new PileupBuilder(dict, mappedPairsOnly=false, minMapQ=10).build(builder.iterator, "chr1", 105)
    pile.depth shouldBe 2
    pile.exists(_.rec.name == "q1") shouldBe false
  }

  it should "filter out reads below the min baseq" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict), baseQuality=19)
    builder.addFrag(name="q1", start=101)
    builder ++= new SamBuilder(readLength=50, sd=Some(dict), baseQuality=20).addFrag(name="q2", start=101)
    builder ++= new SamBuilder(readLength=50, sd=Some(dict), baseQuality=21).addFrag(name="q3", start=101)
    val pile = new PileupBuilder(dict, mappedPairsOnly=false, minBaseQ=20).build(builder.iterator, "chr1", 105)
    pile.depth shouldBe 2
    pile.exists(_.rec.name == "q1") shouldBe false
  }

  it should "filter out reads that are not part of a mapped pair" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    builder.addFrag(name="q1", start=101)
    builder.addPair(name="q2", start1=101, start2=101, unmapped2=true)
    builder.addPair(name="q3", start1=101, start2=300)
    val pile = new PileupBuilder(dict, mappedPairsOnly=true).build(builder.iterator, "chr1", 105)
    pile.depth shouldBe 1
    pile.exists(r => r.rec.name == "q3" && r.rec.firstOfPair) shouldBe true
  }

  it should "filter out reads and bases with custom filters" in {
    val piler = new PileupBuilder(dict, mappedPairsOnly=false)
        .withReadFilter(r => !r.name.startsWith("x"))
        .withBaseFilter(_.offset > 5)
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    builder.addFrag(name="q1", start=100)
    builder.addFrag(name="x2", start=104)
    builder.addFrag(name="q3", start=108)
    builder.addFrag(name="q4", start=112)

    val pile = piler.build(builder.iterator, "chr1", 115)
    pile.depth shouldBe 2
    pile.map(_.rec.name) should contain theSameElementsAs Seq("q1", "q3")
  }

  it should "remove one half of each overlapping pair" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict))
    builder.addPair(name="q1", start1=100, start2=110)
    builder.addPair(name="q2", start1=110, start2=100)
    builder.addPair(name="q3", start1= 50, start2=100)
    val pile = new PileupBuilder(dict, mappedPairsOnly=true).build(builder.iterator, "chr1", 125)

    pile.depth shouldBe 5
    pile.withoutOverlaps.depth shouldBe 3
    pile.withoutOverlaps.groupBy(_.rec.name).values.foreach(xs => xs.size shouldBe 1)
  }

  "PilupRec" should "report correct offsets/positions/bases" in {
    val builder = new SamBuilder(readLength=50, sd=Some(dict), baseQuality=35)
    builder.addPair(name="q1", start1=101, start2=201, bases1 = "A"*50, bases2 = "C"*50)

    val pile1 = new PileupBuilder(dict, mappedPairsOnly=true).build(builder.iterator, "chr1", 105)
    pile1.depth shouldBe 1
    val p1 = pile1.pile.head.asInstanceOf[BaseEntry]
    p1.base shouldBe 'A'
    p1.baseInReadOrientation shouldBe 'A'
    p1.qual shouldBe  35
    p1.offset shouldBe 4
    p1.positionInRead shouldBe 5
    p1.positionInReadInReadOrder shouldBe 5

    val pile2 = new PileupBuilder(dict, mappedPairsOnly=true).build(builder.iterator, "chr1", 205)
    pile2.depth shouldBe 1
    val p2 = pile2.pile.head.asInstanceOf[BaseEntry]
    p2.base shouldBe 'C'
    p2.baseInReadOrientation shouldBe 'G'
    p2.qual shouldBe  35
    p2.offset shouldBe 4
    p2.positionInRead shouldBe 5
    p2.positionInReadInReadOrder shouldBe 46
  }
}
