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

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.Murmur3
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord, SAMRecordCoordinateComparator, SAMRecordQueryNameComparator}

import java.util.Random
import scala.collection.mutable.ListBuffer

object SamOrderTest {
  // Builder to for a set of records to be used in testing sorting
  val builder = new SamBuilder(sort=None)
  val random  = new Random(42)
  Range.inclusive(1, 1000).foreach { _ =>
    builder.addPair(name="q"+random.nextInt(), contig=random.nextInt(23), start1=random.nextInt(1e7.toInt)+1, start2=random.nextInt(1e7.toInt)+1)
  }
  Range.inclusive(1, 10).foreach { _ => builder.addPair(name="q"+random.nextInt(), unmapped1=true, unmapped2=true)}
}

class SamOrderTest extends UnitSpec {
  "SamOrder.apply(String)" should "find the appropriate SamOrder" in {
    SamOrder.values.foreach { order =>
      SamOrder(order.name)             shouldBe order
      SamOrder(order.name.toLowerCase) shouldBe order
      SamOrder(order.name.toUpperCase) shouldBe order
    }
  }

  "SamOrder.apply(SAMFileHeader)" should "find known orders from headers" in {
    val header = new SAMFileHeader()

    header.setSortOrder(SortOrder.coordinate)
    header.setAttribute("GO", null)
    header.setAttribute("SS", null)
    SamOrder(header) shouldBe Some(SamOrder.Coordinate)

    header.setSortOrder(SortOrder.queryname)
    header.setAttribute("GO", null)
    header.setAttribute("SS", null)
    SamOrder(header) shouldBe Some(SamOrder.Queryname)

    header.setSortOrder(SortOrder.unsorted)
    header.setAttribute("GO", null)
    header.setAttribute("SS", null)
    SamOrder(header) shouldBe Some(SamOrder.Unsorted)

    header.setSortOrder(SortOrder.unsorted)
    header.setAttribute("GO", null)
    header.setAttribute("SS", "unsorted:random")
    SamOrder(header) shouldBe Some(SamOrder.Random)

    header.setSortOrder(SortOrder.unsorted)
    header.setGroupOrder(GroupOrder.query)
    header.setAttribute("SS", "unsorted:random-query")
    SamOrder(header) shouldBe Some(SamOrder.RandomQuery)

    header.setSortOrder(SortOrder.unsorted)
    header.setGroupOrder(GroupOrder.query)
    header.setAttribute("GO", GroupOrder.query.name())
    header.setAttribute("SS", "unsorted:template-coordinate")
    SamOrder(header) shouldBe Some(SamOrder.TemplateCoordinate)
  }

  "SamOrder.Coordinate" should "sort reads into coordinate order" in {
    val f = SamOrder.Coordinate.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp = new SAMRecordCoordinateComparator()
    recs.sliding(2).foreach { case Seq(lhs,rhs) => comp.fileOrderCompare(lhs.asSam, rhs.asSam) <= 0 shouldBe true }
  }

  "SamOrder.Queryname" should "sort reads into queryname order" in {
    val f = SamOrder.Queryname.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp = new SAMRecordQueryNameComparator()
    recs.sliding(2).foreach { case Seq(lhs,rhs) => comp.fileOrderCompare(lhs.asSam, rhs.asSam) <= 0 shouldBe true }
  }

  "SamOrder.Random" should "randomize the order of the reads" in {
    val f = SamOrder.Random.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp1 = new SAMRecordCoordinateComparator()
    val comp2 = new SAMRecordQueryNameComparator()
    val counter1 = new SimpleCounter[Int]()
    val counter2 = new SimpleCounter[Int]()

    recs.sliding(2).foreach { case Seq(r1, r2) =>
      counter1.count(comp1.fileOrderCompare(r1.asSam, r2.asSam))
      counter2.count(comp2.fileOrderCompare(r1.asSam, r2.asSam))
    }

    counter1.countOf(0 ) / counter1.total.toDouble <= 0.1 shouldBe true
    counter1.countOf(-1) / counter1.total.toDouble <= 0.6 shouldBe true
    counter1.countOf(1 ) / counter1.total.toDouble <= 0.6 shouldBe true

    counter2.countOf(0 ) / counter2.total.toDouble <= 0.1 shouldBe true
    counter2.countOf(-1) / counter2.total.toDouble <= 0.6 shouldBe true
    counter2.countOf(1 ) / counter2.total.toDouble <= 0.6 shouldBe true
  }

  "SamOrder.RandomQuery" should "keep query names together" in {
    val f = SamOrder.RandomQuery.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val counts = new SimpleCounter[String]
    val iter = recs.iterator.bufferBetter
    while (iter.hasNext) {
      val name = iter.head.name
      val _ = iter.takeWhile(_.name == name).foreach { _ => () }
      counts.count(name)
    }

    counts.foreach { case (_, count) => count shouldBe 1}
  }

  it should "keep querynames together even when there are has collisions" in {
    val f      = SamOrder.RandomQuery.sortkey
    val builder = new SamBuilder()
    val name1 = "000000000-CHG2G:1:1102:14353:13008"
    val name2 = "000000000-CHG2G:1:2108:16511:13017"

    // Show that there is a hash collision
    val hasher = new Murmur3(SamOrder.RandomQuery.HashSeed)
    hasher.hashUnencodedChars(name1) shouldBe hasher.hashUnencodedChars(name2)

    builder.addFrag(name=name1, start=100)
    builder.addFrag(name=name1, start=150).foreach(_.supplementary = true)
    builder.addFrag(name=name2, start=100)
    builder.addFrag(name=name2, start=150).foreach(_.supplementary = true)

    val recs   = builder.toIndexedSeq.sortBy(f(_))
    val counts = new SimpleCounter[String]
    val iter   = recs.iterator.bufferBetter
    while (iter.hasNext) {
      val name = iter.head.name
      iter.takeWhile(_.name == name).foreach { _ => () }
      counts.count(name)
    }

    counts.foreach { case (_, count) => count shouldBe 1 }
  }

  "SamOrder.TemplateCoordinate" should "sort by molecular identifier then name" in {
    val addFuncs: Seq[SamBuilder => Unit] = Seq(
      b => { val _ = b.addPair(name="ab0", start1=200, start2=200, attrs=Map("MI" -> "0/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ab1", start1=100, start2=100, attrs=Map("MI" -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ab2", start1=100, start2=100, attrs=Map("MI" -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ab3", start1=100, start2=100, attrs=Map("MI" -> "2/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ba0", start1=200, start2=200, strand1=Minus, strand2=Plus, attrs=Map("MI" -> "0/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ba1", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map("MI" -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ba2", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map("MI" -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="ba3", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map("MI" -> "2/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
    )

    def seq(n: Int, str: String): Seq[String] = IndexedSeq.fill[String](n)(str)

    Range.inclusive(start=1, end=10).foreach { _ =>
      val builder = new SamBuilder(readLength=10)
      scala.util.Random.shuffle(addFuncs).foreach { func => func(builder) }
      val f = SamOrder.TemplateCoordinate.sortkey
      val recs = builder.iterator.toSeq.sortBy(f(_))
      recs should have length 16
      val moleclarIdentifiers = seq(4, "1/A") ++ seq(4, "1/B") ++ seq(2, "2/A") ++ seq(2, "2/B") ++ seq(2, "0/A") ++ seq(2, "0/B")
      recs.map(_.apply[String]("MI")) should contain theSameElementsInOrderAs moleclarIdentifiers
      val names = Seq("ab1", "ab2", "ba1", "ba2", "ab3", "ba3", "ab0", "ba0").flatMap { name => seq(2, name)}
      recs.map(_.name) should contain theSameElementsInOrderAs names
      recs.grouped(2).foreach { pair =>
        pair.count(_.firstOfPair) shouldBe 1
        pair.count(_.secondOfPair) shouldBe 1
        pair.foreach(_.name shouldBe pair.head.name)
      }
    }
  }

  it should "group molecules within cells when the CB tag is defined" in {
    // Two cells (AAAA, BBBB), each with two molecules sharing the same coordinates.
    // Without CB the sort degenerates to MI-only ordering, which interleaves molecules
    // across cells (AAAA/mol1 adjacent to BBBB/mol1, etc.).  With CB present, reads
    // must be grouped by cell first, then by molecule within the cell.
    val addFuncs: Seq[SamBuilder => Unit] = Seq(
      b => { val _ = b.addPair(name="aa1", start1=100, start2=200, attrs=Map("MI" -> "1/A", "CB" -> "AAAA"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="bb1", start1=100, start2=200, attrs=Map("MI" -> "1/A", "CB" -> "BBBB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="aa2", start1=100, start2=200, attrs=Map("MI" -> "2/A", "CB" -> "AAAA"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
      b => { val _ = b.addPair(name="bb2", start1=100, start2=200, attrs=Map("MI" -> "2/A", "CB" -> "BBBB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA") },
    )

    def seq(n: Int, str: String): Seq[String] = IndexedSeq.fill[String](n)(str)

    Range.inclusive(start=1, end=10).foreach { _ =>
      val builder = new SamBuilder(readLength=10)
      scala.util.Random.shuffle(addFuncs).foreach { func => func(builder) }
      val f    = SamOrder.TemplateCoordinate.sortkey
      val recs = builder.iterator.toSeq.sortBy(f(_))
      recs should have length 8

      // All reads from cell AAAA must precede all reads from cell BBBB.
      recs.take(4).map(_.apply[String]("CB")).foreach(_ shouldBe "AAAA")
      recs.drop(4).map(_.apply[String]("CB")).foreach(_ shouldBe "BBBB")

      // Within each cell, molecules are ordered by MI.
      recs.take(4).map(_.apply[String]("MI")) should contain theSameElementsInOrderAs (seq(2, "1/A") ++ seq(2, "2/A"))
      recs.drop(4).map(_.apply[String]("MI")) should contain theSameElementsInOrderAs (seq(2, "1/A") ++ seq(2, "2/A"))
    }
  }

  it should "sort by library, then cell barcode, before molecular identifier and name" in {
    val builder = new SamBuilder(readLength=10)
    builder.header.getReadGroups.get(0).setLibrary("libA")
    val rgB = new SAMReadGroupRecord("B")
    rgB.setLibrary("libB")
    rgB.setSample("Sample")
    builder.header.addReadGroup(rgB)

    builder.addPair(name="q1", start1=100, start2=200, attrs=Map("RG" -> "B", "CB" -> "AAAA", "MI" -> "1"))
    builder.addPair(name="q2", start1=100, start2=200, attrs=Map("CB" -> "BBBB", "MI" -> "2"))
    builder.addPair(name="q3", start1=100, start2=200, attrs=Map("RG" -> "B", "CB" -> "AAAA", "MI" -> "1"))
    builder.addPair(name="q4", start1=100, start2=200, attrs=Map("CB" -> "AAAA", "MI" -> "3"))

    val actual = builder.toSeq.sortBy(r => SamOrder.TemplateCoordinate.sortkey(r)).map(_.id)
    // libA (q2, q4) before libB (q1, q3); within libA, CB AAAA (q4) before CB BBBB (q2) despite q4's larger MI
    actual should contain theSameElementsInOrderAs Seq("q4/1", "q4/2", "q2/1", "q2/2", "q1/1", "q1/2", "q3/1", "q3/2")
  }

  it should "sort reverse-strand templates before forward-strand templates at the same position" in {
    val builder = new SamBuilder(readLength=10)
    // Both fragments have their 5' end at position 100; the forward one has the earlier name
    builder.addFrag(name="a", start=100, strand=Plus)
    builder.addFrag(name="b", start=91, strand=Minus)

    val actual = builder.toSeq.sortBy(r => SamOrder.TemplateCoordinate.sortkey(r)).map(_.name)
    actual should contain theSameElementsInOrderAs Seq("b", "a")
  }

  it should "sort templates whose upper end is on the reverse strand before those on the forward strand" in {
    val builder = new SamBuilder(readLength=10)
    // Both pairs have a forward lower end at 100 and an upper end whose 5' position is 200; the pair with the
    // forward upper end has the earlier name
    builder.addPair(name="a", start1=100, start2=200, strand1=Plus, strand2=Plus)
    builder.addPair(name="b", start1=100, start2=191, strand1=Plus, strand2=Minus)

    val actual = builder.toSeq.sortBy(r => SamOrder.TemplateCoordinate.sortkey(r)).map(_.id)
    actual should contain theSameElementsInOrderAs Seq("b/1", "b/2", "a/1", "a/2")
  }

  it should "ignore the strand flag of an unmapped end" in {
    val builder = new SamBuilder(readLength=10)
    // Both pairs have a mapped forward R1 at 100 and an unmapped R2 whose strand flag differs
    val Seq(a1, a2) = builder.addPair(name="a", start1=100, unmapped2=true, strand1=Plus, strand2=Plus)
    val Seq(b1, b2) = builder.addPair(name="b", start1=100, unmapped2=true, strand1=Plus, strand2=Minus)
    a2.negativeStrand shouldBe false
    b2.negativeStrand shouldBe true

    val key = SamOrder.TemplateCoordinate.sortkey
    // Apart from the name and which end is upper, every record has the same key
    Seq(a1, a2, b1, b2).map(r => key(r).copy(name="", isUpperOfPair=false)).distinct should have length 1
  }

  it should "sort records with no library before records with a library, as samtools does" in {
    val builder = new SamBuilder(readLength=10)
    // A library name that sorts before "Unknown", the default this order used to give records with no library
    builder.header.getReadGroups.get(0).setLibrary("A")
    val noLibrary = new SAMReadGroupRecord("B")
    noLibrary.setSample("Sample")
    builder.header.addReadGroup(noLibrary)
    builder.addPair(name="a", start1=100, start2=200)
    builder.addPair(name="b", start1=100, start2=200, attrs=Map("RG" -> "B"))

    val actual = builder.toSeq.sortBy(r => SamOrder.TemplateCoordinate.sortkey(r)).map(_.id)
    actual should contain theSameElementsInOrderAs Seq("b/1", "b/2", "a/1", "a/2")
  }

  it should "sort pairs by the 'lower' 5' position of the pair" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    val exp = ListBuffer[SamRecord]()
    // Records are added to the builder in the order that we expect them to be sorted, but the builder
    // will coordinate sort them for us, so we can re-sort them and test the results
    exp ++= builder.addPair("q1", contig=0, start1=100, start2=300)
    exp ++= builder.addPair("q2", contig=0, start1=106, start2=300, cigar1="5S95M") // effective=101
    exp ++= builder.addPair("q3", contig=0, start1=102, start2=299)
    exp ++= builder.addPair("q4", contig=0, start1=300, start2=110, strand1=Minus, strand2=Plus)
    exp ++= builder.addPair("q5", contig=0, start1=120, start2=320)
    exp ++= builder.addPair("q6", contig=1, start1=1,   start2=200)

    // Order they are added in except for q4 gets it's mate's flipped because of strand order
    val expected = List("q1/1", "q1/2", "q2/1", "q2/2", "q3/1", "q3/2", "q4/2", "q4/1", "q5/1", "q5/2", "q6/1", "q6/2")
    val actual   = builder.toList.sortBy(r => SamOrder.TemplateCoordinate.sortkey(r)).map(_.id)

    actual should contain theSameElementsInOrderAs expected
  }
}
