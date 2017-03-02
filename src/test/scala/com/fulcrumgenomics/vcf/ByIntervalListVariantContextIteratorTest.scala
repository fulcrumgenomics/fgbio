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

import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import com.fulcrumgenomics.FgBioDef._

/**
  * Tests for ByIntervalListVariantContextIterator.
  */
class ByIntervalListVariantContextIteratorTest extends UnitSpec {

  private val dict = new VariantContextSetBuilder().header.getSequenceDictionary

  private def emtpyIntervalList(): IntervalList = {
    val header = new SAMFileHeader
    header.setSequenceDictionary(this.dict)
    new IntervalList(header)
  }

  private def toIterator(reader: VCFFileReader,
                         intervalList: IntervalList,
                         useIndex: Boolean): Iterator[VariantContext] = {
    if (useIndex) {
      ByIntervalListVariantContextIterator(reader=reader, intervalList=intervalList)
    }
    else {
      ByIntervalListVariantContextIterator(iterator=reader.iterator(), intervalList=intervalList, dict=reader.getFileHeader.getSequenceDictionary)
    }
  }

  "ByIntervalListVariantContextIterator" should "return no variants if the interval list is empty" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=1, variantAlleles=List("A"), genotypeAlleles=List("A"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator shouldBe 'empty
    }
  }

  it should "return no variants if the variants are empty" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder()
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator shouldBe 'empty
    }
  }

  it should "return a variant context if it overlaps an interval" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=500, variantAlleles=List("A"), genotypeAlleles=List("A"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'nonEmpty
      val actual = iterator.next()
      val expected = builder.head
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator shouldBe 'empty
    }
  }

  it should "not return a variant context if it doesn't overlap an interval (same chromosome)" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=500, variantAlleles=List("A"), genotypeAlleles=List("A"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 750, 1000, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'empty
    }
  }

  it should "not return a variant context if it doesn't overlap an interval (different chromosome)" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=500, variantAlleles=List("A"), genotypeAlleles=List("A"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(1).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'empty
    }
  }

  it should "throw an exception when next() is call but hasNext() is false" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=1, variantAlleles=List("A"), genotypeAlleles=List("A"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.hasNext shouldBe false
      an[Exception] should be thrownBy iterator.next()
    }
  }

  it should "return a variant context if it encloses an interval" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=495, variantAlleles=List("AAAAA", "A"), genotypeAlleles=List("A"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'nonEmpty
      val actual = iterator.next()
      val expected = builder.head
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator shouldBe 'empty
    }
  }

  it should "return a variant context only once if it overlaps multiple intervals" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=495, variantAlleles=List("AAAAA", "A"), genotypeAlleles=List("A"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'nonEmpty
      val actual = iterator.next()
      val expected = builder.head
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator shouldBe 'empty
    }
  }

  it should "throw an exception when intervals are given out of order when using the VCF index" in {
    val builder = new VariantContextSetBuilder()
      .addVariant(refIdx=0, start=495, variantAlleles=List("AAAAA", "A"), genotypeAlleles=List("A"))
      .addVariant(refIdx=0, start=595, variantAlleles=List("AAAAA", "A"), genotypeAlleles=List("A"))
    val intervalList = emtpyIntervalList()
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 494, 500, false, "foo"))
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
    val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=true)
    // OK, since we are overlapping the first interval
    iterator shouldBe 'nonEmpty
    // NOK, since the intervals were overlapping when we pre-fetch the second variant context
    an[Exception] should be thrownBy iterator.next()
  }

  it should "ignore a variant context if does not overlap an interval" in {
    Stream(true, false).foreach { useIndex =>
      val builder = new VariantContextSetBuilder().addVariant(refIdx=0, start=495, variantAlleles=List("A", "C"), genotypeAlleles=List("C"))
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(reader=builder.toVcfFileReader(), intervalList=intervalList, useIndex=useIndex)
      iterator shouldBe 'empty
    }
  }
}
