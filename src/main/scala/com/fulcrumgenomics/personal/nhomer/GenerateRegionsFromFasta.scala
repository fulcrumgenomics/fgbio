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

package com.fulcrumgenomics.personal.nhomer

import java.nio.file.Path

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.ClpGroups.Personal
import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.commons.CommonsDef.PathToFasta
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

@clp(
  description="Generates a list of freebayes/bamtools region specifiers.",
  group = Personal
)
class GenerateRegionsFromFasta
( @arg(doc = "The input FASTA.") val input: PathToFasta,
  @arg(doc = "The output.") val output: Path = Io.StdOut,
  @arg(doc = "The size of the regions to output.") val regionSize: Int = 100000
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(input.toFile)
    val sequenceDictionary = referenceSequenceFile.getSequenceDictionary
    sequenceDictionary.getSequences.foreach { sequenceRecord =>
      for (regionStart <- 0 until sequenceRecord.getSequenceLength by regionSize) {
        val regionEnd = Math.min(regionStart + regionSize, sequenceRecord.getSequenceLength)
       println(s"${sequenceRecord.getSequenceName}:$regionStart-$regionEnd")
      }
    }
  }
}
