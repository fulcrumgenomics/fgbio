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
  * Tests for JointVariantIterator.
  */
class JointVariantIteratorTest extends UnitSpec {
  private val dict = VcfBuilder.DefaultHeader.dict

  /** A builder for a single-sample VCF. */
  private def newBuilder(): VcfBuilder = VcfBuilder(Seq("Sample"))

  private def compareVariants(actual: Variant, expected: Variant) = {
    actual.chrom shouldBe expected.chrom
    actual.pos shouldBe expected.pos
    actual.end shouldBe expected.end
  }

  "JointVariantIterator" should "iterate variants given a single iterator" in {
    val builder = newBuilder().add(chrom=dict(0).name, pos=1, alleles=Seq("A"))
    val iterator = JointVariantIterator(iters=Seq(builder.iterator), dict=dict)
    compareVariants(actual=iterator.next().head.get, expected=builder.head)
  }

  it should "not return a variant if all the iterators are empty" in {
    val builder = newBuilder()
    val iterator = JointVariantIterator(iters=Seq(builder.iterator, builder.iterator), dict=dict)
    iterator.hasNext shouldBe false
    an[NoSuchElementException] should be thrownBy iterator.next()
  }

  it should "return a pair of variants at the same position" in {
    val builder = newBuilder().add(chrom=dict(0).name, pos=1, alleles=Seq("A"))
    val iterator = JointVariantIterator(iters=Seq(builder.iterator, builder.iterator), dict=dict)
    iterator.hasNext shouldBe true
    val Seq(left, right) = iterator.next().flatten
    compareVariants(left, right)
  }

  it should "return a None for an iterator that doesn't have a variant for a given covered site" in {
    import scala.language.existentials
    val builderLeft = newBuilder()
      .add(chrom=dict(0).name, pos=10, alleles=Seq("A"))
      .add(chrom=dict(0).name, pos=30, alleles=Seq("A"))
    val builderRight = newBuilder()
      .add(chrom=dict(0).name, pos=10, alleles=Seq("A"))
      .add(chrom=dict(0).name, pos=20, alleles=Seq("A"))
    val iterator = JointVariantIterator(iters=Seq(builderLeft.iterator, builderRight.iterator), dict=dict)
    // pos: 10 status: both
    iterator.hasNext shouldBe true
    iterator.next().flatten match {
      case Seq(left, right) => compareVariants(left, right)
    }
    // pos: 20 status: right
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(None, Some(right)) => compareVariants(right, builderRight.last)
    }
    // pos: 30 status: left
    iterator.hasNext shouldBe true
    iterator.next() match {
      case Seq(Some(left), None) => compareVariants(left, builderLeft.last)
    }
  }

  it should "order variants by the index of their contig in the sequence dictionary, then by position" in {
    val left: VcfBuilder  = newBuilder().add(chrom=dict(0).name, pos=10, alleles=Seq("A")).add(chrom=dict(1).name, pos=5, alleles=Seq("A"))
    val right: VcfBuilder = newBuilder().add(chrom=dict(1).name, pos=5, alleles=Seq("A")).add(chrom=dict(1).name, pos=7, alleles=Seq("A"))
    val iterator = JointVariantIterator(iters=Seq(left.iterator, right.iterator), dict=dict)
    iterator.map(_.map(_.map(v => (v.chrom, v.pos)))).toSeq shouldBe Seq(
      Seq(Some((dict(0).name, 10)), None),
      Seq(Some((dict(1).name, 5)), Some((dict(1).name, 5))),
      Seq(None, Some((dict(1).name, 7)))
    )
  }
}
