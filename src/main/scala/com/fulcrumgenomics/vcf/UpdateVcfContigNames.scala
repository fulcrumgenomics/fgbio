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

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToVcf, javaIterableToIterator}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.vcf.{VCFFileReader, VCFHeader}

@clp(description =
  """
    |Updates then contig names in a VCF.
    |
    |The input sequence dictionary should have the output contig names in the SN filed, and the source contig names in
    |the AN field (comma-separated for multiple values).  For example, the contigs with names "NC_000001.11" and
    |"CM000663.2" will be updated to "chr1" if the following line is present in the sequence dictionary:
    |`@SQ SN:chr1 LN:248956422 AN:NC_000001.11,CM000663.2`.
  """,
  group = ClpGroups.VcfOrBcf)
class UpdateVcfContigNames
(@arg(flag='i', doc="Input VCF.") val input: PathToVcf,
 @arg(flag='d', doc="The sequence dictionary") val dict: FilePath,
 @arg(flag='o', doc="Output VCF.") val output: PathToVcf,
 @arg(doc="Skip missing contigs.") val skipMissing: Boolean = false
) extends FgBioTool with LazyLogging{

  override def execute(): Unit = {
    // read in the dictionary
    val dict              = SAMSequenceDictionaryExtractor.extractDictionary(this.dict)
    val altNamesToPrimary = dict.getSequences.flatMap { seq =>
      val primary = seq.getSequenceName
      val alts    = seq.getAttribute("AN") match {
        case null => throw new IllegalStateException(s"Sequence '$primary' does not have the AN attribute set.")
        case value => value.split(',')
      }
      alts.map { alt => (alt, primary) }
    }.toMap

    // build the VCF reader and writer
    val reader = new VCFFileReader(this.input)
    val header = {
      val h: VCFHeader = new VCFHeader(reader.getFileHeader)
      h.setSequenceDictionary(dict)
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
       altNamesToPrimary.get(v.getContig) match {
        case None =>
          if (skipMissing) logger.warning(s"Did not find contig ${v.getContig} in the sequence dictionary.")
          else throw new IllegalStateException(s"Did not find contig ${v.getContig}  in the sequence dictionary.")
        case Some(contig) =>
          val newV = new VariantContextBuilder(v)
              .chr(contig)
              .make()
          progress.record(newV.getContig, newV.getStart)
          writer.add(newV)
       }
    }
    progress.logLast()

    writer.close()
  }
}
