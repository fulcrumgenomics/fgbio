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

/**
  * Created by tfenne on 6/23/16.
  */

// import com.fulcrumgenomics.FgBioDef._


import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sorter}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.Murmur3

@clp(group=ClpGroups.SamOrBam, description=
"""
  |Randomizes the order of reads in a SAM or BAM file. Randomization is done by sorting
  |on a hash of the `queryname` (and bases and quals if not query-grouping). By default
  |reads with the same query name are grouped together in the output file; this can be
  |turned off by specifying --query-group=false.
""")
class RandomizeBam(
  @arg(flag='i', doc="The input SAM or BAM file.")           val input: PathToBam = Io.StdIn,
  @arg(flag='o', doc="The output SAM or BAM file.")          val output: PathToBam = Io.StdOut,
  @arg(flag='s', doc="Random seed.")                         val seed: Int = 42,
  @arg(flag='q', doc="Group together reads by queryname.")   val queryGroup: Boolean = true
) extends FgBioTool with LazyLogging {

  Io.assertCanWriteFile(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in = SamSource(input)
    val hasher = new Murmur3(seed)

    val sorter = {
      if (queryGroup) {
        val f = (r: SamRecord) => SamOrder.RandomQueryKey(hasher.hashUnencodedChars(r.name), r.name, r.asSam.getFlags)
        new Sorter(2e6.toInt, new SamRecordCodec(in.header), f)
      }
      else {
        val f = (r: SamRecord) => {
          val key = s"${r.name}:${r.basesString}:${r.qualsString}"
          SamOrder.RandomKey(hasher.hashUnencodedChars(key), r.asSam.getFlags)
        }
        new Sorter(2e6.toInt, new SamRecordCodec(in.header), f)
      }
    }

    logger.info("Sorting reads into random order.")
    val sortProgress = ProgressLogger(logger, verb="sorted", unit=5e6.toInt)
    in.foreach(rec => {
      sorter += rec
      sortProgress.record()
    })

    logger.info("Writing reads out to output.")
    val outHeader = in.header.clone()
    (if (queryGroup) SamOrder.RandomQuery else SamOrder.Random).applyTo(outHeader)
    outHeader.setSortOrder(SortOrder.unsorted)
    if (queryGroup) outHeader.setGroupOrder(GroupOrder.query)
    val out = SamWriter(output, outHeader)

    out ++= sorter
    out.close()
  }
}
