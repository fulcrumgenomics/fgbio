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
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.Defaults
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}

/**
  * Writes [[Variant]]s to a file or other storage mechanism.
  *
  * @param writer the underlying HTSJDK writer
  * @param header the header of the VCF
  */
class VcfWriter private(private val writer: VariantContextWriter, val header: VcfHeader) extends Writer[Variant] {
  override def write(variant: Variant): Unit = writer.add(VcfConversions.toJavaVariant(variant, header))
  override def close(): Unit = writer.close()
}


object VcfWriter {
  var DefaultUseAsyncIo: Boolean = Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE

  /**
    * Creates a [[VcfWriter]] that will write to the give path.  The path must end in either
    *   - `.vcf` to create an uncompressed VCF file
    *   - `.vcf.gz` to create a block-gzipped VCF file
    *   - `.bcf` to create a binary BCF file
    *
    * @param path the path to write to
    * @param header the header of the VCF
    * @return a VariantWriter to write to the given path
    */
  def apply(path: PathToVcf, header: VcfHeader, async: Boolean = DefaultUseAsyncIo): VcfWriter = {
    val javaHeader = VcfConversions.toJavaHeader(header)

    val builder = new VariantContextWriterBuilder()
      .setOutputPath(path)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setReferenceDictionary(header.dict)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .setBuffer(Io.bufferSize)

    if (async) builder.setOption(Options.USE_ASYNC_IO) else builder.unsetOption(Options.USE_ASYNC_IO)
    val writer = builder.build()
    writer.writeHeader(javaHeader)
    new VcfWriter(writer, header)
  }
}
