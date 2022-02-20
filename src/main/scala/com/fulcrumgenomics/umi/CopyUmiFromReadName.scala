/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef.{PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}

@clp(group=ClpGroups.Umi, description=
  """
    |Copies the UMI at the end of the BAM's read name to the RX tag.
    |
    |The read name is split by the given name delimiter, and the last field is assumed to be the UMI sequence.  The UMI
    |will be copied to the `RX` tag as per the SAM specification.
  """)
class CopyUmiFromReadName
( @arg(flag='i', doc="The input BAM file") input: PathToBam,
  @arg(flag='o', doc="The output BAM file") output: PathToBam,
  @arg(doc="Remove the UMI from the read name") removeUmi: Boolean = false,
  @arg(doc="Replaces any occurrences of this delimiter found in the UMI with a dash ('-') as per the SAM specification")
  umiDelimiter: Option[Char] = None
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val source   = SamSource(input)
    val writer   = SamWriter(output, source.header)
    val progress = new ProgressLogger(logger)
    source.foreach { rec =>
      progress.record(rec)
      writer += Umis.copyUmiFromReadName(rec=rec, removeUmi=removeUmi, umiDelimiter=umiDelimiter)
    }
    progress.logLast()
    source.safelyClose()
    writer.close()
  }
}
