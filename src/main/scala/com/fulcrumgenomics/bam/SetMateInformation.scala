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
package com.fulcrumgenomics.bam

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.SamPairUtil.SetMateInfoIterator
import htsjdk.samtools.{SAMFileWriterFactory, SamReaderFactory}

@clp(group=ClpGroups.SamOrBam, description=
"""
  |Adds and/or fixes mate information on paired-end reads. Sets the MQ (mate mapping quality),
  |MC (mate cigar string), ensures all mate-related flag fields are set correctly, and that
  |the mate reference and mate start position are correct.
  |
  |Supplementary records are handled correctly (updated with their mate's non-supplemental
  |attributes).  Secondary alignments are passed through but are not updated.
  |
  |The input file must be query-name sorted or query-name grouped (i.e. all records from the same
  |query sequence must be adjacent in the file, though the ordering between queries is unspecified).
""")
class SetMateInformation
(
  @arg(flag='i', doc="Input SAM/BAM/CRAM file.")                      val input: PathToBam = Io.StdIn,
  @arg(flag='o', doc="Output SAM/BAM/CRAM file.")                     val output: PathToBam = Io.StdOut,
  @arg(flag='r', doc="Reference fasta, only needed if writing CRAM.") val ref: Option[PathToFasta] = None,
  @arg(flag='x', doc="If specified, do not fail when reads marked as paired are missing their mate pairs.")
                                                                      val allowMissingMates: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  private val in = SamReaderFactory.make().open(input.toFile)
  if (in.getFileHeader.getSortOrder != SortOrder.queryname && in.getFileHeader.getGroupOrder != GroupOrder.query) {
    in.safelyClose()
    throw new ValidationException("Input is not queryname sorted or grouped.")
  }

  override def execute(): Unit = {
    val out = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, ref.map(_.toFile).orNull)
    val progress = new ProgressLogger(logger)
    val iterator = new SetMateInfoIterator(in.iterator(), true, allowMissingMates).toIterator

    for (rec <- iterator) {
      out.addAlignment(rec)
      progress.record(rec)
    }

    out.close()
  }
}
