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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

/**
  * Tests for ByIntervalListVariantContextIterator.
  */
class ByIntervalListVariantContextIteratorTest extends UnitSpec {

  private val dict = VcfBuilder(samples=Seq("s1")).header.dict.toSam
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
      import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
      val dict = reader.getFileHeader.getSequenceDictionary.fromSam
      ByIntervalListVariantContextIterator(iterator=reader.iterator(), intervalList=intervalList, dict=dict)
    }
  }

  "ByIntervalListVariantContextIterator" should "return no variants if the interval list is empty" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=1, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="0/0")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val iterator   = toIterator(reader=builder, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "return no variants if the variants are empty" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1"))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "return a variant context if it overlaps an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="0/0")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.iterator().next()
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator.isEmpty shouldBe true
    }
  }

  it should "not return a variant context if it doesn't overlap an interval (same chromosome)" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="0/0")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 750, 1000, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "not return a variant context if it doesn't overlap an interval (different chromosome)" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=500, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="0/0")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(1).getSequenceName, 1, 1000, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }

  it should "throw an exception when next() is call but hasNext() is false" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=1, alleles=Seq("A"), gts=Seq(Gt(sample="s1", gt="0/0")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val iterator = toIterator(reader=builder, intervalList=emtpyIntervalList(), useIndex=useIndex)
      iterator.hasNext shouldBe false
      an[Exception] should be thrownBy iterator.next()
    }
  }

  it should "return a variant context if it encloses an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="s1", gt="1/1")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.iterator().next()
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator.isEmpty shouldBe true
    }
  }

  it should "return a variant context only once if it overlaps multiple intervals" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="s1", gt="1/1")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 496, 496, false, "foo"))
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe false
      val actual = iterator.next()
      val expected = builder.iterator().next()
      actual.getContig shouldBe expected.getContig
      actual.getStart shouldBe expected.getStart
      actual.getEnd shouldBe expected.getEnd
      iterator.isEmpty shouldBe true
    }
  }

  it should "throw an exception when intervals are given out of order when using the VCF index" in {
    val vcfBuilder = VcfBuilder(samples=Seq("s1"))
      .add(chrom="chr1", pos=495, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="s1", gt="1/1")))
      .add(chrom="chr1", pos=595, alleles=Seq("AAAAA", "A"), gts=Seq(Gt(sample="s1", gt="1/1")))
    val builder    = new VCFFileReader(vcfBuilder.toTempFile())
    val intervalList = emtpyIntervalList()
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 494, 500, false, "foo"))
    intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
    val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=true)
    // OK, since we are overlapping the first interval
    iterator.isEmpty shouldBe false
    // NOK, since the intervals were overlapping when we pre-fetch the second variant context
    an[Exception] should be thrownBy iterator.next()
  }

  it should "ignore a variant context if does not overlap an interval" in {
    Iterator(true, false).foreach { useIndex =>
      val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=495, alleles=Seq("A", "C"), gts=Seq(Gt(sample="s1", gt="1/1")))
      val builder    = new VCFFileReader(vcfBuilder.toTempFile())
      val intervalList = emtpyIntervalList()
      intervalList.add(new Interval(dict.getSequence(0).getSequenceName, 500, 500, false, "foo"))
      val iterator = toIterator(reader=builder, intervalList=intervalList, useIndex=useIndex)
      iterator.isEmpty shouldBe true
    }
  }
}
