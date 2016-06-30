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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{IlluminaAdapters, Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.util.SequenceUtil
import htsjdk.samtools.{SAMFileWriterFactory, SamReaderFactory}

import scala.annotation.tailrec
import scala.collection.JavaConversions.asScalaIterator

@clp(group = ClpGroups.SamOrBam, description =
  """
     |Find reads that are from technical or synthetic sequences in a BAM file. Takes in
     |a BAM file, extracts the read pairs and fragment reads that are unmapped, and tests
     |them to see if they are likely generated from a technical sequence (e.g. adapter
     |dimer).
     |
     |The identification of reads is done by testing the first N bases (controlled by the
     |match-length parameter) of each read against all sub-sequences of length N from the
     |technical sequences.  Sub-sequences are generated from both the sequences and the
     |reverse complement of the sequences, ignoring any sub-sequences that include 'N's.
     |
     |The output BAM file will contain all reads that matched to a sub-sequence of the
     |technical sequences and, if the read is paired, the read's mate pair.
     |
     |The default set of sequences include a range of different Illumina adapter sequences
     |with the sample index/barcode region masked to Ns.
  """
  )
class FindTechnicalReads
( @arg(flag="i", doc="Input SAM or BAM file") val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file") val output: PathToBam,
  @arg(flag="m", doc="The number of bases at the start of the read to match against.") val matchLength: Int = 15,
  @arg(flag="e", doc="The maximum number of errors in the matched region.") val maxErrors: Int = 1,
  @arg(flag="s", doc="The set of technical sequences to look for.") val sequences: Seq[String] = IlluminaAdapters.all.flatMap(_.both)
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in       = SamReaderFactory.make().open(input.toFile)
    val iter     = if (in.hasIndex) in.queryUnmapped().buffered
                   else in.iterator().filter(r => r.getReadUnmappedFlag && (!r.getReadPairedFlag || r.getMateUnmappedFlag)).buffered
    val out      = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, null)
    val matcher  = new Matcher(sequences=sequences, matchLength=matchLength, maxErrors=maxErrors)
    val progress = new ProgressLogger(logger, unit=250000)
    var n: Long  = 0

    while (iter.hasNext) {
      val r1 = iter.next()
      val r2 = if (iter.hasNext && iter.head.getReadName == r1.getReadName) Some(iter.next()) else None
      val recs = Seq(Some(r1), r2).flatten
      if (recs.exists(r => matcher.matches(r.getReadBases))) {
        recs.foreach(out.addAlignment)
        n += recs.size
      }

      recs.foreach(progress.record)
    }

    logger.info(s"Wrote out ${n} reads that matched to technical sequences.")
    in.safelyClose()
    out.close()
  }
}

/** Little utility class to match reads to technical sequences. */
private[bam] class Matcher(val matchLength: Int, val maxErrors: Int, sequences: Seq[String]) {
  private val kmers: Array[Array[Byte]] =
    sequences.map(_.toUpperCase).flatMap(a => Seq(a, SequenceUtil.reverseComplement(a)))
      .flatMap(a => a.sliding(matchLength))
      .filterNot(_.contains('N'))
      .map(_.getBytes)
      .toSet.toArray

  def matches(read: Array[Byte]): Boolean = {
    read.length >= matchLength && kmers.exists(kmer => matches(read=read, kmer=kmer, len=matchLength))
  }

  @tailrec
  private def matches(read: Array[Byte], kmer: Array[Byte], index: Int = 0, len: Int, errors: Int = 0): Boolean = {
    if (errors > maxErrors) false
    else if (index >= len)  true
    else matches(read, kmer, index+1, len, errors + (if(read(index) == kmer(index)) 0 else 1))
  }
}
