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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.FgBioDef.{PathToSequenceDictionary, PathToVcf, SafelyClosable, javaIterableToIterator}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.vcf.{VCFFileReader, VCFHeader}

@clp(description =
  """
    |Updates then contig names in a VCF.
    |
    |The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
    |new name will be the primary (non-alias) name in the sequence dictionary.
  """,
  group = ClpGroups.VcfOrBcf)
class UpdateVcfContigNames
(@arg(flag='i', doc="Input VCF.") val input: PathToVcf,
 @arg(flag='d', doc="The path to the sequence dictionary with contig aliases.") val dict: PathToSequenceDictionary,
 @arg(flag='o', doc="Output VCF.") val output: PathToVcf,
 @arg(doc="Skip missing contigs.") val skipMissing: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, dict))
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val dict   = SequenceDictionary(this.dict)
    val reader = new VCFFileReader(this.input)
    val header = {
      import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
      val h: VCFHeader = new VCFHeader(reader.getFileHeader)
      h.setSequenceDictionary(dict.asSam)
      h
    }
    val writer = {
      new VariantContextWriterBuilder()
        .setOutputPath(this.output)
        .setOption(Options.INDEX_ON_THE_FLY)
        .build()
    }
    writer.writeHeader(header)

    // go through all the records
    val progress = ProgressLogger(logger, noun = "variants", verb = "written")
    reader.foreach { v =>
      dict.get(v.getContig) match {
        case None =>
          if (skipMissing) logger.warning(s"Did not find contig ${v.getContig} in the sequence dictionary.")
          else throw new IllegalStateException(s"Did not find contig ${v.getContig} in the sequence dictionary.")
        case Some(info) =>
          val newV = new VariantContextBuilder(v).chr(info.name).make()
          progress.record(newV.getContig, newV.getStart)
          writer.add(newV)
       }
    }
    progress.logLast()

    reader.safelyClose()
    writer.close()
  }
}
