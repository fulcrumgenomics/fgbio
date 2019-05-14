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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{FilePath, _}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import htsjdk.samtools.reference.{FastaSequenceIndex, ReferenceSequenceFileFactory}

@clp(description =
  """
    |Updates the sequence names in a FASTA.
    |
    |The first column of the input is the source name and the second column is an target name.  If there
    |is more than one target names (ex. multiple alternates), each alternate name may be on a separate line or as
    |additional columns.  Only the first target name will be considered.  This is mainly to support the output of
    |`CollectAlternateContigNames`.
  """,
  group = ClpGroups.Fasta)
class UpdateFastaContigNames
(@arg(flag='i', doc="Input FASTA.") val input: PathToFasta,
 @arg(flag='m', doc="The path to the source to target contig names.") val mapping: FilePath,
 @arg(flag='o', doc="Output FASTA.")val output: PathToFasta,
 @arg(flag='l', doc="Line length or sequence lines.") val lineLength: Int = 100,
 @arg(doc="Skip missing source contigs.") val skipMissing: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, mapping))
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, noun="bases", verb="written", unit=10e7.toInt)
    val refFile  = ReferenceSequenceFileFactory.getReferenceSequenceFile(input, true, true)
    val seen     = new scala.collection.mutable.HashSet[String]()
    val out      = Io.toWriter(output)

    val srcContigs = refFile.getSequenceDictionary match {
      case null =>
        require(refFile.isIndexed,
          "Reference sequence file must have a sequence dictionary or be indexed.  Try 'picard CreateSequenceDictionary' or 'samtools faidx <ref.fasta>'.")
        new FastaSequenceIndex(ReferenceSequenceFileFactory.getFastaIndexFileName(this.input)) map(_.getContig)
      case dict => dict.getSequences.map(_.getSequenceName)
    }

    Io.readLines(mapping).flatMap { line =>
      val fields = line.split('\t')
      require(fields.length <= 2, s"Malformed line: expected at least two columns: $line")
      val srcName = fields(0)
      val targetName = fields(1)
      if (seen.contains(srcName)) None else {
        seen.add(srcName)
        Some((srcName, targetName))
      }
    }.foreach { case (srcName, targetName) =>
        val ref = refFile.getSequence(srcName)
        out.append('>').append(targetName).append('\n')
        ref.getBases.grouped(lineLength).zipWithIndex.foreach { case (bases, groupIndex) =>
          val start = groupIndex * lineLength
          bases.zipWithIndex.foreach { case (base, baseIdx) =>
            progress.record(targetName, start + baseIdx + 1)
            out.write(base)
          }
          out.newLine()
        }
        progress.logLast()
        seen.add(srcName)
    }
    out.close()

    // log those that weren't seen
    srcContigs.foreach { srcName =>
      if (!seen.contains(srcName)) {
        if (skipMissing) logger.warning(s"Did not find contig $srcName in the list of original names.")
        else throw new IllegalStateException(s"Did not find contig $srcName in the list of original names.")
      }
    }
  }
}