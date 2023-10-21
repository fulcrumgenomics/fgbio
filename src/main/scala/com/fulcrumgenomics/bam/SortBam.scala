/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sorter}

@clp(group=ClpGroups.SamOrBam, description =
  """
    |Sorts a SAM or BAM file.  Several sort orders are available:
    |
    |1. **Coordinate**: sorts reads by their reference sequence and left-most aligned coordinate
    |2. **Queryname**: sort the reads by their query (i.e. read) name
    |3. **Random**: sorts the reads into a random order. The output is deterministic for any given input.
    |and several
    |4. **RandomQuery**: sorts the reads into a random order but keeps reads with the same
    |   queryname together. The ordering is deterministic for any given input.
    |
    |Uses a temporary directory to buffer sets of sorted reads to disk. The number of reads kept in memory
    |affects memory use and can be changed with the `--max-records-in-ram` option.  The temporary directory
    |to use can be set with the fgbio global option `--tmp-dir`.
    |
    |An example invocation might look like:
    |
    |```bash
    |java -Xmx4g -jar fgbio.jar --tmp-dir=/my/big/scratch/volume \
    |  SortBam --input=queryname.bam --sort-order=Coordinate --output coordinate.bam
    |```
  """)
class SortBam
( @arg(flag='i', doc="Input SAM or BAM.")   val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM.")  val output: PathToBam,
  @arg(flag='s', doc="Order into which to sort the records.") val sortOrder: SamOrder = SamOrder.Coordinate,
  @arg(flag='m', doc="Max records in RAM.") val maxRecordsInRam: Int = SamWriter.DefaultMaxRecordsInRam,
  @arg(flag='t', doc="Number of threads to use.") val threads: Int = 4
) extends FgBioTool with LazyLogging {
  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)

    val in     = SamSource(input)
    val header = in.header.clone()
    sortOrder.applyTo(header)
    val sorter = Bams.sorter(sortOrder, header, maxRecordsInRam=maxRecordsInRam, threads=threads)
    val progress = ProgressLogger(logger, verb="sorted", unit=2e6.toInt)
    in.foreach { r =>
      sorter += r
      progress.record(r)
    }

    val out = SamWriter(output, header=header)
    sorter.foreach { r => out += r }
    out.close()
  }
}
