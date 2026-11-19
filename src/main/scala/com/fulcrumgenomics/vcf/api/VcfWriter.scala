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

import java.nio.file.Files
import java.nio.file.attribute.BasicFileAttributes

/**
  * Writes [[Variant]]s to a file or other storage mechanism.
  *
  * @param writer the underlying HTSJDK writer
  * @param header the header of the VCF
  * @param isBcf true if writing a BCF, false if writing a VCF
  */
class VcfWriter private(private val writer: VariantContextWriter, val header: VcfHeader, isBcf: Boolean) extends Writer[Variant] {
  override def write(variant: Variant): Unit = writer.add(VcfConversions.toJavaVariant(variant, header, isBcf=isBcf))
  override def close(): Unit = writer.close()
}


object VcfWriter {
  var DefaultUseAsyncIo: Boolean = Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE

  /**
    * Creates a [[VcfWriter]] that will write to the given path. If the path is meant to point to a regular file, then
    * the path must end in either:
    *
    *   - `.vcf`: to create an uncompressed VCF file
    *   - `.vcf.gz`: to create a block-gzipped VCF file
    *   - `.bcf`: to create a binary BCF file
    *
    * If `index` is true and the path is meant to point to a regular file, then the output is indexed on the fly, which
    * requires that variants be written in coordinate order. However, if the path already exists and the path is not a
    * file or symbolic link, then this function will assume the path is a named pipe or device (such as `/dev/null`)
    * and indexing will not occur.
    *
    * @param path the path to write to
    * @param header the header of the VCF
    * @param index if true, index the output while writing; must be false to write variants out of coordinate order
    * @param writeFullFormatField if true, write every FORMAT field of each sample, writing `.` for missing values, rather
    *                             than leaving out the trailing missing FORMAT fields of each sample
    * @return a VCF writer to write to the given path
    */
  def apply(path: PathToVcf,
            header: VcfHeader,
            async: Boolean = DefaultUseAsyncIo,
            index: Boolean = true,
            writeFullFormatField: Boolean = false): VcfWriter = {
    val javaHeader = VcfConversions.toJavaHeader(header)
    require(!Files.isDirectory(path), s"Path cannot be a directory! Found $path")

    // NB: a header with no contig lines gives an empty dictionary, which HTSJDK cannot index a block-gzipped VCF
    // against (it needs no dictionary at all), but which it needs to index a plain-text VCF.
    val outputType = VariantContextWriterBuilder.determineOutputTypeFromFile(path)
    val dict       = if (header.dict.isEmpty && outputType != VariantContextWriterBuilder.OutputType.VCF) null else header.dict.toSam

    val builder = new VariantContextWriterBuilder()
      .setOutputPath(path)
      .setReferenceDictionary(dict)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .setBuffer(Io.bufferSize)

    if (async) builder.setOption(Options.USE_ASYNC_IO) else builder.unsetOption(Options.USE_ASYNC_IO)
    if (writeFullFormatField) builder.setOption(Options.WRITE_FULL_FORMAT_FIELD)

    // Don't index if asked not to, or if the path exists and is not a file or symbolic link (assumed a named pipe).
    if (!index || (Files.exists(path) && Files.readAttributes(path, classOf[BasicFileAttributes]).isOther)) {
      builder.unsetOption(Options.INDEX_ON_THE_FLY)
      builder.setIndexCreator(null)
    } else {
      builder.setOption(Options.INDEX_ON_THE_FLY)
    }

    val writer = builder.build()
    writer.writeHeader(javaHeader)
    // NB: decide whether the output is a BCF as the writer builder does, which ignores case and follows symbolic links
    new VcfWriter(writer, header, isBcf=outputType == VariantContextWriterBuilder.OutputType.BCF)
  }
}
