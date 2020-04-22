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
import com.fulcrumgenomics.util.ContigNameMapping

@clp(description =
  """
    |Updates the sequence names in a FASTA.
    |
    |The first column of the input is the source name and the second column is an target name.  If there
    |are more than one target names (ex. multiple alternates), each alternate name may be on a separate line or as
    |additional columns.  Only the first target name will be considered.  This is mainly to support the output of
    |`CollectAlternateContigNames`. The output fasta will need to be indexed.
  """,
  group = ClpGroups.Fasta)
class UpdateFastaContigNames
(@arg(flag='i', doc="Input FASTA.") val input: PathToFasta,
 @arg(flag='m', doc="The path to the source the contig name mappings.") val mapping: FilePath,
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
    val out      = Io.toWriter(output)

    val srcContigs = refFile.getSequenceDictionary match {
      case null =>
        require(refFile.isIndexed,
          "Reference sequence file must have a sequence dictionary or be indexed.  Try 'picard CreateSequenceDictionary' or 'samtools faidx <ref.fasta>'.")
        new FastaSequenceIndex(ReferenceSequenceFileFactory.getFastaIndexFileName(this.input)) map(_.getContig)
      case dict => dict.getSequences.map(_.getSequenceName)
    }
    val srcToTarget = ContigNameMapping.parse(mapping)

    srcContigs.foreach { srcName =>
      srcToTarget.get(srcName) match {
        case None if skipMissing => logger.warning(s"Did not find contig $srcName in the list of original names.")
        case None                => throw new IllegalStateException(s"Did not find contig $srcName in the list of original names.")
        case Some(targetNames)   => {
          val ref = refFile.getSequence(srcName)
          out.append('>').append(targetNames.head).append('\n')
          val bases = ref.getBases
          var baseCounter = 0
          forloop(from = 0, until = bases.length) { baseIdx =>
            out.write(bases(baseIdx))
            progress.record(targetNames.head, baseIdx + 1)
            baseCounter += 1
            if (baseCounter >= lineLength) {
              out.newLine()
              baseCounter = 0
            }
          }
          if (baseCounter > 0) out.newLine()
        }
      }
    }
    out.close()
  }

}