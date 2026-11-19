/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.{Interval, IntervalList}

/**
  * Tests for ByIntervalListVariantIterator.
  */
class ByIntervalListVariantIteratorTest extends UnitSpec {

  private val dict = VcfBuilder.DefaultHeader.dict.toSam

  /** A builder for a single-sample VCF. */
  private def newBuilder(): VcfBuilder = VcfBuilder(Seq("Sample"))

  private def emtpyIntervalList(): IntervalList = {
    val header = new SAMFileHeader
    header.setSequenceDictionary(this.dict)
    new IntervalList(header)
  }

  private def toIterator(source: VcfSource,
                         intervalList: IntervalList,
                         useIndex: Boolean): Iterator[Variant] = {
    if (useIndex) {
      ByIntervalListVariantIterator(source=source, intervalList=intervalList)
    }
    else {
      ByIntervalListVariantIterator(iterator=source.iterator, intervalList=intervalList, dict=source.header.dict)
    }
  }

  "ByIntervalListVariantIterator" should "return no variants if the interval list is empty" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=1, alleles=Seq("A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val iterator = toIterator(source=builder.toSource, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "return no variants if the variants are empty" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder()
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "return a variant if it overlaps an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.head
      actual.chrom shouldBe expected.chrom
      actual.pos shouldBe expected.pos
      actual.end shouldBe expected.end
      iterator.isEmpty shouldBe true
    }
  }

  it should "not return a variant if it doesn't overlap an interval (same chromosome)" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 750, 1000, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "not return a variant if it doesn't overlap an interval (different chromosome)" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(1).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "throw an exception when next() is call but hasNext() is false" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=1, alleles=Seq("A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val iterator = toIterator(source=builder.toSource, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.hasNext shouldBe false
      an[Exception] should be thrownBy iterator.next()
    }
  }

  it should "return a variant if it encloses an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.head
      actual.chrom shouldBe expected.chrom
      actual.pos shouldBe expected.pos
      actual.end shouldBe expected.end
      iterator.isEmpty shouldBe true
    }
  }

  it should "return a variant only once if it overlaps multiple intervals" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="Sample", gt="A")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.head
      actual.chrom shouldBe expected.chrom
      actual.pos shouldBe expected.pos
      actual.end shouldBe expected.end
      iterator.isEmpty shouldBe true
    }
  }

  it should "throw an exception when intervals are given out of order when using the VCF index" in {
    val builder: VcfBuilder = newBuilder()
      .add(chrom=dict.getSequence(0).getSequenceName, pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="Sample", gt="A")))
      .add(chrom=dict.getSequence(0).getSequenceName, pos=595, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="Sample", gt="A")))
    val intervalList = emtpyIntervalList()
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 494, 500, false, "foo"))
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
    val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=true)
    // OK, since we are overlapping the first interval
    iterator.isEmpty shouldBe false
    // NOK, since the intervals were overlapping when we pre-fetch the second variant
    an[Exception] should be thrownBy iterator.next()
  }

  it should "ignore a variant if does not overlap an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val builder = newBuilder().add(chrom=dict.getSequence(0).getSequenceName, pos=495, alleles=Seq("A", "C"), gts=Seq(Gt(sample="Sample", gt="C")))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(source=builder.toSource, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }
}
