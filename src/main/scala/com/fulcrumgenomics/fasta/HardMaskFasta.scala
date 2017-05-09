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
package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import htsjdk.samtools.reference.{ReferenceSequenceFileFactory, ReferenceSequenceFileWalker}

/**
  * Tool to take in a FASTA file and convert soft-masked (i.e. lower-case) sequence to hard-masked
  * sequence (i.e. Ns), for tools that cannot correctly interpret soft-masking.
  *
  * @author Tim Fennell
  */
@clp(description =
  """
     |Converts soft-masked sequence to hard-masked in a FASTA file. All lower case bases are
     |converted to Ns, all other bases are left unchanged.  Line lengths are also standardized
     |to allow easy indexing with "samtools faidx".
  """,
  group = ClpGroups.Fasta)
class HardMaskFasta
( @arg(flag='i', doc="Input FASTA file.")              val input: PathToFasta,
  @arg(flag='o', doc="Output FASTA file.")             val output: PathToFasta,
  @arg(flag='l', doc="Line length or sequence lines.") val lineLength: Int = 100
 ) extends FgBioTool with LazyLogging {

  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)

    val out      = Io.toWriter(output)
    val progress = new ProgressLogger(logger, noun="sequences", verb="masked", unit=5000)

    for (seq <- ReferenceSequenceIterator(input)) {
      out.append('>').append(seq.getName).append('\n')
      seq.getBases.grouped(lineLength).foreach(bs => {
        bs.indices.foreach(i => if (bs(i).toChar.isLower) bs(i) = 'N')
        bs.foreach(out.write(_))
        out.newLine()
      })

      progress.record()
    }

    out.close()
  }
}
