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
import com.fulcrumgenomics.commons.io.PathUtil
import htsjdk.samtools.util.CloseableIterator
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.variant.bcf2.BCF2Codec
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.{VCFCodec, VCFFileReader, VCFHeader}

/**
  * Provides a reader over a source of VCF-like data that could be a VCF file or a BCF file. Has facilities
  * for iterating over the entire stream of variants as well as querying by genomic location if an
  * index is present.
  *
  * @param reader the underlying HTSJDK [[VCFFileReader]].
  * @param headerTransformer a method to transform the header after being read in.
  * @param allowKindMismatch allow a mismatch between the actual kind and the kind defined in the VCF header.
  * @param allowExtraFields allow fields not defined in the VCF header.
  */
class VcfSource private(private val reader: AbstractFeatureReader[VariantContext, _],
                        private val headerTransformer: VcfHeader => VcfHeader,
                        allowKindMismatch: Boolean = false,
                        allowExtraFields: Boolean = false) extends View[Variant] with Closeable {
  /** The header associated with the VCF being read. */
  val header: VcfHeader = headerTransformer(VcfConversions.toScalaHeader(reader.getHeader.asInstanceOf[VCFHeader]))

  /**
    * The type of iterator returned by both the [[iterator]] method as well as the [[query()]] method. Note that
    * [[SelfClosingIterator]] both self-closes when it hits the end of the iterator, _and_ extends
    * [[com.fulcrumgenomics.commons.collection.BetterBufferedIterator]].
    */
  type VariantIterator = SelfClosingIterator[Variant]

  /** Wraps an iterator provided by HTSJDK into a SelfClosingIterator that transforms VariantContexts into Variants. */
  private def wrap(it: CloseableIterator[VariantContext]): VariantIterator = {
    new SelfClosingIterator(
      iter   = it.map(vc => VcfConversions.toScalaVariant(vc, header, allowKindMismatch=allowKindMismatch, allowExtraFields=allowExtraFields)),
      closer = () => it.close())
  }

  /**
    * Returns an iterator over the entire stream of variants. The returned iterator may be be closed by invoking
    * `close()` on it, and will automatically close itself when exhausted.  Only a single iterator at a time
    * is supported per [[VcfSource]], including iterators returned from [[query()]].
    */
  override def iterator: VariantIterator = wrap(reader.iterator())

  /** True if the VCF is sorted and indexed such that queries can be executed, false otherwise. */
  def isQueryable: Boolean = this.reader.isQueryable

  /**
    * Returns an iterator over variants overlapping the specified genomic region.
    *
    * The returned iterator may be be closed by invoking `close()` on it, and will automatically close itself
    * when exhausted.  Only a single iterator at a time is supported per [[com.fulcrumgenomics.vcf.api.VcfSource]], including iterators
    * returned from [[iterator()]].
    */
  def query(chrom: String, start: Int, end: Int): Iterator[Variant] = wrap(reader.query(chrom, start, end))

  /** Closes the underlying reader. */
  override def close(): Unit = this.reader.safelyClose()

  /** Required for 2.12 compat. */
  protected def underlying: VcfSource = this
}

object VcfSource {
  /**
    * Manufactures a variant source for reading from the specified path.  The index, if one exists, will be
    * auto-discovered based on the path to the VCF.
    *
    * @param path   the path to a VCF, gzipped VCF or BCF file
    * @return a VariantSource for reading from the path given
    */
  def apply(path: PathToVcf): VcfSource = {
    this.apply(path=path, headerTransformer=identity)
  }

  /**
    * Manufactures a variant source for reading from the specified path.  The index, if one exists, will be
    * auto-discovered based on the path to the VCF.
    *
    * @param path the path to a VCF, gzipped VCF or BCF file
    * @param headerTransformer a method to transform the header after being read in
    * @param allowKindMismatch allow a mismatch between the actual kind and the kind defined in the VCF header.
    * @param allowExtraFields allow fields not defined in the VCF header.
    * @return a VariantSource for reading from the path given
    */
  def apply(path: PathToVcf,
            headerTransformer: VcfHeader => VcfHeader,
            allowKindMismatch: Boolean = false,
            allowExtraFields: Boolean = false): VcfSource = {
    val codec  = if (PathUtil.extensionOf(path).contains(".bcf")) {
      new BCF2Codec
    }
    else {
      val c = new VCFCodec
      c.disableOnTheFlyModifications()
      c
    }

    val reader = AbstractFeatureReader.getFeatureReader(path.toUri.toString, codec, false)

    new VcfSource(
      reader            = reader,
      headerTransformer = headerTransformer,
      allowKindMismatch = allowKindMismatch,
      allowExtraFields  = allowExtraFields
    )
  }
}

