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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import htsjdk.samtools.util.CloseableIterator
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.variant.bcf2.BCF2Codec
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.{VCFCodec, VCFHeader, VariantFileType}

import java.io.Closeable

/**
  * Provides a reader over a source of VCF-like data that could be a VCF file or a BCF file. Has facilities
  * for iterating over the entire stream of variants as well as querying by genomic location if an
  * index is present.
  *
  * Reading a variant with a value of a typed INFO or FORMAT field that is not a valid value of its declared type (e.g.
  * `1.5` or `3000000000` for an Integer field) throws an [[java.lang.IllegalArgumentException]] naming the type, the
  * value, the field, the variant's position and the VCF.  This does not cover the standard `GQ`, `DP`, `AD` and `PL`
  * FORMAT fields, which HTSJDK parses itself: an invalid value of one of these either fails with HTSJDK's error, or is
  * changed or dropped without an error.
  *
  * @param reader the underlying HTSJDK [[VCFFileReader]]
  * @param headerTransformer a method to transform the header after being read in.
  * @param allowMissingFieldsInHeader if true, INFO and FORMAT fields that are not described in the header are read
  *                                   as [[VcfFieldType.String]] fields of [[VcfCount.Unknown]] count (or as a flag,
  *                                   and the INFO END field as an Int), otherwise reading a variant with such a field
  *                                   throws an exception.
  * @param description a description of the VCF (e.g. its path) used in error messages
  */
class VcfSource private(private val reader: AbstractFeatureReader[VariantContext, _],
                        private val headerTransformer: VcfHeader => VcfHeader,
                        private val allowMissingFieldsInHeader: Boolean,
                        private val description: String) extends View[Variant] with Closeable {
  /** The header associated with the VCF being read. */
  val header: VcfHeader = headerTransformer(VcfConversions.toScalaHeader(reader.getHeader.asInstanceOf[VCFHeader]))

  /**
    * The type of iterator returned by both the [[iterator]] method as well as the [[query]] method. Note that
    * [[com.fulcrumgenomics.commons.collection.SelfClosingIterator]] both self-closes when it hits the end of the iterator, _and_ extends
    * [[com.fulcrumgenomics.commons.collection.BetterBufferedIterator]].
    */
  type VariantIterator = SelfClosingIterator[Variant]

  /** Wraps an iterator provided by HTSJDK into a [[com.fulcrumgenomics.commons.collection.SelfClosingIterator]] that transforms VariantContexts into Variants. */
  private def wrap(it: CloseableIterator[VariantContext]): VariantIterator = {
    new SelfClosingIterator(
      iter   = it.map { vc =>
        try VcfConversions.toScalaVariant(vc, header, allowMissingFieldsInHeader) catch {
          case ex: InvalidVcfValueException => throw new InvalidVcfValueException(s"${ex.getMessage} in $description", ex)
        }
      },
      closer = () => it.close())
  }

  /**
    * Returns an iterator over the entire stream of variants. The returned iterator may be closed by invoking
    * `close()` on it, and will automatically close itself when exhausted.  Only a single iterator at a time
    * is supported per [[VcfSource]], including iterators returned from [[query]].
    */
  override def iterator: VariantIterator = wrap(reader.iterator())

  /** True if the VCF is sorted and indexed such that queries can be executed, false otherwise. */
  def isQueryable: Boolean = this.reader.isQueryable

  /**
    * Returns an iterator over variants overlapping the specified genomic region.
    *
    * The returned iterator may be closed by invoking `close()` on it, and will automatically close itself
    * when exhausted.  Only a single iterator at a time is supported per [[com.fulcrumgenomics.vcf.api.VcfSource]], including iterators
    * returned from [[iterator]].
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
    * @param path the path to a VCF, gzipped VCF or BCF file
    * @return a VariantSource for reading from the path given
    */
  def apply(path: PathToVcf): VcfSource = {
    this.apply(path=path, headerTransformer=identity, allowMissingFieldsInHeader=false)
  }

  /**
    * Manufactures a variant source for reading from the specified path.  The index, if one exists, will be
    * auto-discovered based on the path to the VCF.
    *
    * @param path the path to a VCF, gzipped VCF or BCF file
    * @param headerTransformer a method to transform the header after being read in
    * @return a VariantSource for reading from the path given
    */
  def apply(path: PathToVcf, headerTransformer: VcfHeader => VcfHeader): VcfSource = {
    this.apply(path=path, headerTransformer=headerTransformer, allowMissingFieldsInHeader=false)
  }

  /**
    * Manufactures a variant source for reading from the specified path.  The index, if one exists, will be
    * auto-discovered based on the path to the VCF.
    *
    * @param path the path to a VCF, gzipped VCF or BCF file
    * @param headerTransformer a method to transform the header after being read in
    * @param allowMissingFieldsInHeader if true, INFO and FORMAT fields that are not described in the header are read
    *                                   as [[VcfFieldType.String]] fields of [[VcfCount.Unknown]] count (or as a
    *                                   flag, and the INFO END field as an Int), as HTSJDK's `VCFFileReader` does,
    *                                   otherwise reading a variant with such a field throws an exception
    * @return a VariantSource for reading from the path given
    */
  def apply(path: PathToVcf,
            headerTransformer: VcfHeader => VcfHeader = identity,
            allowMissingFieldsInHeader: Boolean = false): VcfSource = {
    // NB: recognize a BCF by its name as HTSJDK does, ignoring case
    val codec  = if (VariantFileType.fromFileName(path.getFileName.toString).orElse(null) == VariantFileType.BCF) {
      new BCF2Codec
    }
    else {
      val c = new VCFCodec
      c.disableOnTheFlyModifications()
      c
    }

    val reader = AbstractFeatureReader.getFeatureReader(path.toUri.toString, codec, false)
    new VcfSource(reader, headerTransformer=headerTransformer, allowMissingFieldsInHeader=allowMissingFieldsInHeader, description=path.toString)
  }

  /** Return the only sample in the VCF source otherwise raise an exception. */
  def onlySample(source: VcfSource): String = {
    source.header.samples.toList match {
      case head :: Nil => head
      case _ => throw new IllegalArgumentException(s"Source is not single-sample. Found samples: ${source.header.samples.mkString(", ")}")
    }
  }
}
