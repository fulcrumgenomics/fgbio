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

import com.fulcrumgenomics.testing.{UnitSpec, VariantContextSetBuilder, VcfBuilder}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.jdk.CollectionConverters.IteratorHasAsScala

/**
  * Tests for JointVariantContextIterator.
  */
class JointVariantContextIteratorTest extends UnitSpec {
  private val dict = VcfBuilder(samples=Seq("s1")).header.dict

  private def compareVariantContexts(actual: VariantContext, expected: VariantContext): Unit = {
    actual.getContig shouldBe expected.getContig
    actual.getStart shouldBe expected.getStart
    actual.getEnd shouldBe expected.getEnd
  }

  "JointVariantContextIterator" should "iterate variant contexts given a single iterator" in {
    val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=1, alleles=Seq("A"))
    val builder    = new VCFFileReader(vcfBuilder.toTempFile())

    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator().asScala), dict=dict)
    compareVariantContexts(actual=iterator.next().head.get, expected=builder.iterator().next())
  }

  it should "not return a variant context if all the iterators are empty" in {
    val vcfBuilder = VcfBuilder(samples=Seq("s1"))
    val builder    = new VCFFileReader(vcfBuilder.toTempFile())

    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator().asScala, builder.iterator().asScala), dict=dict)
    iterator.hasNext shouldBe false
    an[NoSuchElementException] should be thrownBy iterator.next()
  }

  it should "return a pair of variant contexts at the same position" in {
    val vcfBuilder = VcfBuilder(samples=Seq("s1")).add(chrom="chr1", pos=1, alleles=Seq("A"))
    val builder    = new VCFFileReader(vcfBuilder.toTempFile())

    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator().asScala, builder.iterator().asScala), dict=dict)
    iterator.hasNext shouldBe true
    val Seq(left, right) = iterator.next().flatten
    compareVariantContexts(left, right)
  }

  it should "return a None for an iterator that doesn't have a variant context for a given covered site" in {
    val vcfBuilderLeft = VcfBuilder(samples=Seq("s1"))
        .add(chrom="chr1", pos=10, alleles=Seq("A"))
        .add(chrom="chr1", pos=30, alleles=Seq("A"))
    val builderLeft    = new VCFFileReader(vcfBuilderLeft.toTempFile())

    val vcfBuilderRight = VcfBuilder(samples=Seq("s1"))
      .add(chrom="chr1", pos=10, alleles=Seq("A"))
      .add(chrom="chr1", pos=20, alleles=Seq("A"))
    val builderRight    = new VCFFileReader(vcfBuilderRight.toTempFile())

    val iterator = JointVariantContextIterator(iters=Seq(builderLeft.iterator().asScala, builderRight.iterator().asScala), dict=dict)
    // pos: 10 status: both
    iterator.hasNext shouldBe true
    iterator.next().flatten match {
      case Seq(left, right) => compareVariantContexts(left, right)
    }
    // pos: 20 status: right
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(None, Some(right)) => compareVariantContexts(right, builderRight.iterator().asScala.toSeq.last)
    }
    // pos: 30 status: left
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(Some(left), None) => compareVariantContexts(left, builderLeft.iterator().asScala.toSeq.last)
    }
  }
}
