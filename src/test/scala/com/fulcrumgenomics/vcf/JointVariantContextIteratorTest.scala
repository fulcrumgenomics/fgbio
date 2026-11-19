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

import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}
import com.fulcrumgenomics.vcf.api.Variant

/**
  * Tests for JointVariantContextIterator.
  */
class JointVariantContextIteratorTest extends UnitSpec {
  private val dict = VcfBuilder.DefaultHeader.dict

  /** A builder for a single-sample VCF. */
  private def newBuilder(): VcfBuilder = VcfBuilder(Seq("Sample"))

  private def compareVariantContexts(actual: Variant, expected: Variant) = {
    actual.chrom shouldBe expected.chrom
    actual.pos shouldBe expected.pos
    actual.end shouldBe expected.end
  }

  "JointVariantContextIterator" should "iterate variant contexts given a single iterator" in {
    val builder = newBuilder().add(chrom=dict(0).name, pos=1, alleles=Seq("A"))
    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator), dict=dict)
    compareVariantContexts(actual=iterator.next().head.get, expected=builder.head)
  }

  it should "not return a variant context if all the iterators are empty" in {
    val builder = newBuilder()
    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator, builder.iterator), dict=dict)
    iterator.hasNext shouldBe false
    an[NoSuchElementException] should be thrownBy iterator.next()
  }

  it should "return a pair of variant contexts at the same position" in {
    val builder = newBuilder().add(chrom=dict(0).name, pos=1, alleles=Seq("A"))
    val iterator = JointVariantContextIterator(iters=Seq(builder.iterator, builder.iterator), dict=dict)
    iterator.hasNext shouldBe true
    val Seq(left, right) = iterator.next().flatten
    compareVariantContexts(left, right)
  }

  it should "return a None for an iterator that doesn't have a variant context for a given covered site" in {
    import scala.language.existentials
    val builderLeft = newBuilder()
      .add(chrom=dict(0).name, pos=10, alleles=Seq("A"))
      .add(chrom=dict(0).name, pos=30, alleles=Seq("A"))
    val builderRight = newBuilder()
      .add(chrom=dict(0).name, pos=10, alleles=Seq("A"))
      .add(chrom=dict(0).name, pos=20, alleles=Seq("A"))
    val iterator = JointVariantContextIterator(iters=Seq(builderLeft.iterator, builderRight.iterator), dict=dict)
    // pos: 10 status: both
    iterator.hasNext shouldBe true
    iterator.next().flatten match {
      case Seq(left, right) => compareVariantContexts(left, right)
    }
    // pos: 20 status: right
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(None, Some(right)) => compareVariantContexts(right, builderRight.last)
    }
    // pos: 30 status: left
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(Some(left), None) => compareVariantContexts(left, builderLeft.last)
    }
  }

  it should "order variants with the given ordering" in {
    val left: VcfBuilder  = newBuilder().add(chrom=dict(0).name, pos=10, alleles=Seq("A")).add(chrom=dict(1).name, pos=5, alleles=Seq("A"))
    val right: VcfBuilder = newBuilder().add(chrom=dict(0).name, pos=10, alleles=Seq("A")).add(chrom=dict(1).name, pos=5, alleles=Seq("A"))
    // order by the reverse of the contig index, so that variants on the second contig come first
    val ordering = Ordering.by((v: Variant) => (-dict(v.chrom).index, v.pos))
    val iters    = Seq(left, right).map(_.toSeq.sorted(ordering).iterator)
    val iterator = JointVariantContextIterator(iters=iters, ordering=ordering)
    iterator.map(_.flatten.map(v => (v.chrom, v.pos))).toSeq shouldBe Seq(
      Seq((dict(1).name, 5), (dict(1).name, 5)),
      Seq((dict(0).name, 10), (dict(0).name, 10))
    )
  }
}
