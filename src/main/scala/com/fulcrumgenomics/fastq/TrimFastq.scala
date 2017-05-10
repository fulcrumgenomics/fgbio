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
package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.cmdline.{FgBioTool, ClpGroups}
import com.fulcrumgenomics.util.ProgressLogger
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._

@clp(
  description =
    """
      |Trims reads in one or more line-matched fastq files to a specific read length. The
      |individual fastq files are expected to have the same set of reads, as would be the
      |case with an `r1.fastq` and `r2.fastq` file for the same sample.
      |
      |Optionally supports dropping of reads across all files when one or more reads
      |is already shorter than the desired trim length.
      |
      |Input and output fastq files may be gzipped.
    """,
  group=ClpGroups.Fastq
)
class TrimFastq
(   @arg(flag='i', doc="One or more input fastq files.") val input:  Seq[PathToFastq],
    @arg(flag='o', doc="A matching number of output fastq files.") val output: Seq[PathToFastq],
    @arg(flag='l', doc="Length to trim reads to.")             val length: Int,
    @arg(flag='x', doc="Exclude reads below the trim length.") val exclude: Boolean = false
) extends FgBioTool with LazyLogging {

  validate(input.size == output.size, "Number of input and output files must match.")

  override def execute(): Unit = {
    var discarded: Long = 0
    val progress = new ProgressLogger(this.logger, noun="records", verb="Wrote")

    val sources = input.map(FastqSource(_))
    val writers = output.map(FastqWriter(_))
    while (allHaveNext(sources)) {
      val recs = sources.map(_.next())
      if (exclude && recs.exists(_.length < length)) {
        discarded += 1
      }
      else {
        writers.iterator.zip(recs.iterator).foreach { case(w, r) =>
          w.write(r.trimmedTo(length))
          progress.record()
        }
      }
    }

    writers.foreach(_.close())
    logger.info(s"Wrote ${progress.getCount} records and discarded ${discarded}.")
  }

  /** Checks to see that all the input sources have another record, or that none do. */
  private def allHaveNext(xs: Seq[Iterator[_ <: Any]]): Boolean = {
    xs.map(_.hasNext).distinct match {
      case Seq(true)  => true
      case Seq(false) => false
      case _          => throw new IllegalStateException("Fastq files did not contain same number of records.")
    }
  }
}
