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

package com.fulcrumgenomics.basecalling

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{ProgressLogger, Sorter}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.{BAMRecordCodec, SAMFileWriterFactory, SamReaderFactory}

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Sorts a SAM or BAM file.
  """)
class SortSam
( @arg(flag="i", doc="Input SAM or BAM file.")  val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag="s", doc="Output sort order.")      val sortOrder: SortOrder
) extends FgBioTool with LazyLogging {

  override def execute(): Unit = {
    val in = SamReaderFactory.make().open(input)
    val header = in.getFileHeader.clone()
    header.setSortOrder(sortOrder)
    val sorter = new Sorter(new BAMRecordCodec(header), sortOrder.getComparatorInstance, memory=Some(1e9.toInt))
    val progress1 = new ProgressLogger(logger, verb="sorted")
    in.foreach { r =>
      sorter.append(r)
      progress1.record(r)
    }
    in.safelyClose()

    val out = new SAMFileWriterFactory().makeWriter(header, true, output.toFile, null)
    val progress2 = new ProgressLogger(logger, verb="written")
    sorter.foreach { r =>
      out.addAlignment(r)
      progress2.record(r)
    }
    out.close()
  }
}
