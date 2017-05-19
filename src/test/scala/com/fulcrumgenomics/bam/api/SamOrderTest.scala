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

import java.util.Random

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.{SAMFileHeader, SAMRecordCoordinateComparator, SAMRecordQueryNameComparator}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}

object SamOrderTest {
  // Builder to for a set of records to be used in testing sorting
  val builder = new SamBuilder(sort=None)
  val random  = new Random(42)
  Range.inclusive(1, 1000).foreach { i =>
    builder.addPair(name="q"+random.nextInt(), contig=random.nextInt(23), start1=random.nextInt(1e7.toInt)+1, start2=random.nextInt(1e7.toInt)+1)
  }
  Range.inclusive(1, 10).foreach { i => builder.addPair(name="q"+random.nextInt(), unmapped1=true, unmapped2=true)}
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
    header.setAttribute("GO", GroupOrder.query.name())
    header.setAttribute("SS", "unsorted:random-query")
    SamOrder(header) shouldBe Some(SamOrder.RandomQuery)

    header.setSortOrder(SortOrder.unsorted)
    header.setAttribute("GO", GroupOrder.query.name())
    header.setAttribute("SS", "unsorted:template-coordinate")
    SamOrder(header) shouldBe Some(SamOrder.TemplateCoordinate)
  }

  "SamOrder.Coordinate" should "sort reads into coordinate order" in {
    val f = SamOrder.Coordinate.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp = new SAMRecordCoordinateComparator()
    recs.grouped(2).foreach { case Seq(lhs,rhs) => comp.fileOrderCompare(lhs.asSam, rhs.asSam) <= 0 shouldBe true }
  }

  "SamOrder.Queryname" should "sort reads into queryname order" in {
    val f = SamOrder.Queryname.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp = new SAMRecordQueryNameComparator()
    recs.grouped(2).foreach { case Seq(lhs,rhs) => comp.fileOrderCompare(lhs.asSam, rhs.asSam) <= 0 shouldBe true }
  }

  "SamOrder.Random" should "randomize the order of the reads" in {
    val f = SamOrder.Random.sortkey
    val recs = SamOrderTest.builder.iterator.toSeq.sortBy(f(_))
    val comp1 = new SAMRecordCoordinateComparator()
    val comp2 = new SAMRecordQueryNameComparator()
    val counter1 = new SimpleCounter[Int]()
    val counter2 = new SimpleCounter[Int]()

    recs.grouped(2).foreach { case Seq(r1, r2) =>
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
      val remaining = iter.takeWhile(_.name == name).foreach { r => Unit }
      counts.count(name)
    }

    counts.foreach { case (name, count) => count shouldBe 1}
  }
}
