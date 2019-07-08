/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import java.io.Closeable

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import htsjdk.samtools.util.CloseableIterator
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.IterableView

// class SamSource private(private val reader: SamReader) extends IterableView[SamRecord, SamSource] with HeaderHelper with Closeable {
class VariantSource private (private val reader: VCFFileReader) extends IterableView[Variant, VariantSource] with Closeable {
  val header: VcfHeader = VcfConversions.toScalaHeader(reader.getFileHeader)
  type VariantIterator = SelfClosingIterator[Variant]

  override def close(): Unit = this.reader.safelyClose()
  override protected def underlying: VariantSource = this
  private def wrap(it: CloseableIterator[VariantContext]): VariantIterator = {
    new SelfClosingIterator(
      iter   = it.map(vc => VcfConversions.toScalaVariant(vc, header)),
      closer = () => it.close())
  }

  override def iterator: VariantIterator = wrap(reader.iterator())
  def query(chrom: String, start: Int, end: Int): Iterator[Variant] = wrap(reader.query(chrom, start, end))
}

object VariantSource {
  def apply(path: PathToVcf): VariantSource = {
    val reader = new VCFFileReader(path, false)
    new VariantSource(reader)
  }
}

