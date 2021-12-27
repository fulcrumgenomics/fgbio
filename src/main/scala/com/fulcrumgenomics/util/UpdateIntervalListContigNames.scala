/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef.PathToSequenceDictionary
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.util.Interval

@clp(description =
  """
    |Updates the sequence names in an Interval List file.
    |
    |The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
    |new name will be the primary (non-alias) name in the sequence dictionary.
  """,
  group = ClpGroups.Fasta)
class UpdateIntervalListContigNames
(@arg(flag='i', doc="Input interval list.") val input: PathToIntervals,
 @arg(flag='d', doc="The path to the sequence dictionary with contig aliases.") val dict: PathToSequenceDictionary,
 @arg(flag='o', doc="Output interval list.") val output: PathToIntervals,
 @arg(doc="Skip missing source contigs.") val skipMissing: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, dict))
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, noun="bases", verb="written", unit=10e7.toInt)
    val dict     = SequenceDictionary(this.dict)
    val source   = IntervalListSource(this.input)
    val writer   = IntervalListWriter(path=this.output, dict=dict)
    source.foreach { interval =>
      dict.get(interval.getContig) match {
        case Some(info)          =>
          val newInterval = new Interval(
            info.name,
            interval.getStart,
            interval.getEnd,
            interval.isNegativeStrand,
            interval.getName
          )
          progress.record(newInterval)
          writer += newInterval
        case None =>
          val message = s"Did not find contig '${interval.getContig}' in source mappings."
          if (skipMissing) logger.warning(message)
          else throw new IllegalStateException(message)
      }


    }
    progress.logLast()
    writer.close()
    source.safelyClose()
  }
}