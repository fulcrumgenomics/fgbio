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


import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.FgBioDef._
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.{Murmur3, SortingCollection}

@clp(group=ClpGroups.SamOrBam, description=
"""
  |Randomizes the order of reads in a SAM or BAM file. Randomization is done by sorting
  |on a hash of the queryname (and bases and quals if not query-grouping). By default
  |reads with the same query name are grouped together in the output file; this can be
  |turned off by specifying --query-group=false.
""")
class RandomizeBam(
  @arg(flag="i", doc="The input SAM or BAM file.")           val input: PathToBam = Io.StdIn,
  @arg(flag="o", doc="The output SAM or BAM file.")          val output: PathToBam = Io.StdOut,
  @arg(flag="s", doc="Random seed.")                         val seed: Int = 42,
  @arg(flag="q", doc="Group together reads by queryname.")   val queryGroup: Boolean = true,
  @arg(flag="t", doc="Temporary directory for sorting.")     val tempDirectory: DirPath = Io.defaultTempDir()
) extends FgBioTool with LazyLogging {

  Io.assertCanWriteFile(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val extractor  = if (queryGroup) (rec: SAMRecord) => rec.getReadName
                     else            (rec: SAMRecord) => rec.getReadName + ":" + rec.getReadString + ":" + rec.getBaseQualityString
    val in = SamReaderFactory.make().open(input.toFile)
    val header = in.getFileHeader
    val comparator = new HashComparator(seed, extractor)
    val sorter = SortingCollection.newInstance(classOf[SAMRecord], new BAMRecordCodec(header), comparator, 2e6.toInt, tempDirectory.toFile)

    logger.info("Sorting reads into random order.")
    val sortProgress = new ProgressLogger(logger, verb="sorted", unit=5e6.toInt)
    in.foreach(rec => {
      sorter.add(rec)
      sortProgress.record()
    })

    logger.info("Writing reads out to output.")
    val writeProgress = new ProgressLogger(logger, verb="written", unit=5e6.toInt)
    val outHeader = header.clone()
    outHeader.setSortOrder(SortOrder.unsorted)
    if (queryGroup) outHeader.setGroupOrder(GroupOrder.query)
    val out = new SAMFileWriterFactory().makeWriter(outHeader, true, output.toFile, null)

    sorter.foreach(rec => {
      out.addAlignment(rec)
      writeProgress.record(rec)
    })

    out.close()
  }
}

/** Comparator that compares record based on a hash of extracted fields. */
private class HashComparator(val seed: Int = 42, val extractor: SAMRecord => String) extends SAMRecordQueryNameComparator {
  private val hasher: Murmur3 = new Murmur3(seed)
  private val HashAttributeKey = hasher.hashUnencodedChars("SecretMagicWord")

  override def fileOrderCompare(lhs: SAMRecord, rhs: SAMRecord): Int = hash(lhs).compareTo(hash(rhs))

  /** Extracts fields from a SAMRecord and hashes them, caching the result for re-use. */
  def hash(rec: SAMRecord): Int = {
    var result = rec.getTransientAttribute(HashAttributeKey).asInstanceOf[Integer]
    if (result == null) {
      result = this.hasher.hashUnencodedChars(extractor(rec))
      rec.setTransientAttribute(HashAttributeKey, result)
    }

    result
  }
}
