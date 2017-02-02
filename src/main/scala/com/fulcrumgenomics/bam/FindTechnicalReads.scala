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
import htsjdk.samtools.{SAMFileWriterFactory, SAMRecord, SamReaderFactory}

import scala.annotation.tailrec
import scala.collection.JavaConversions.asScalaIterator

@clp(group = ClpGroups.SamOrBam, description =
  """
     |Find reads that are from technical or synthetic sequences in a BAM file. Takes in
     |a BAM file, extracts the read pairs and fragment reads that are unmapped, and tests
     |them to see if they are likely generated from a technical sequence (e.g. adapter
     |dimer).
     |
     |The technical sequences, supplied using the --sequences parameter, can be supplied
     |one of two ways:
     |  --sequences name1:bases1 name2:bases2 ...
     |  --sequences bases1 bases2
     |If the latter style is used the sequence itself is re-used as the name of the sequence.
     |
     |The identification of reads is done by testing the first N bases (controlled by the
     |match-length parameter) of each read against all sub-sequences of length N from the
     |technical sequences.  Sub-sequences are generated from both the sequences and the
     |reverse complement of the sequences, ignoring any sub-sequences that include 'N's.
     |
     |By default the output BAM file will contain all reads that matched to a sub-sequence of the
     |technical sequences and, if the read is paired, the read's mate pair.  An option is
     |available to apply a tag to matched reads (--tag/-t), and if specified each matching
     |read will be tagged with the names of the sequences to which it matched. In
     |combination with tagging it is possible to output all reads (-a/--all-reads) which will
     |re-create the input BAM with the addition of tags on matching reads.
     |
     |The default set of sequences include a range of different Illumina adapter sequences
     |with the sample index/barcode region masked to Ns.
  """
  )
class FindTechnicalReads
( @arg(flag="i", doc="Input SAM or BAM file") val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file") val output: PathToBam,
  @arg(flag="m", doc="The number of bases at the start of the read to match against.") val matchLength: Int = 15,
  @arg(flag="e", doc="The maximum number of errors allowed in the matched region.") val maxErrors: Int = 1,
  @arg(flag="s", doc="The set of technical sequences to look for.") val sequences: Seq[String] =
    IlluminaAdapters.all.flatMap(a => a.both.map(bases => a.name + ":" + bases)),
  @arg(flag="a", doc="Output all reads.") val allReads: Boolean = false,
  @arg(flag="t", doc="Tag to set to indicate a read is a technical sequence.") val tag: Option[String] = None
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  if (allReads && tag.isEmpty) invalid("Tag option must be used when outputting all reads.")

  override def execute(): Unit = {
    // Built the appropriate input iterator
    val in = SamReaderFactory.make().open(input.toFile)
    val iter = (allReads, in.hasIndex) match {
      case (true, _)      => in.iterator().buffered
      case (false, true)  => in.queryUnmapped().buffered
      case (false, false) => in.iterator().filter(isUnmappedTemplate).buffered
    }

    val out      = new SAMFileWriterFactory().makeWriter(in.getFileHeader, true, output.toFile, null)
    val matcher  = new Matcher(sequences=sequences.map(NamedSequence(_)), matchLength=matchLength, maxErrors=maxErrors)
    val progress = new ProgressLogger(logger, unit=250000)
    var n: Long  = 0

    while (iter.hasNext) {
      val r1 = iter.next()
      val r2 = if (iter.hasNext && iter.head.getReadName == r1.getReadName) Some(iter.next()) else None
      val recs = Seq(Some(r1), r2).flatten

      val isTechnical = isUnmappedTemplate(r1) && {
        recs.flatMap(r => matcher.findMatches(r.getReadBases)) match {
          case Nil  =>
            false
          case hits =>
            for (t <- tag; r <- recs) r.setAttribute(t, hits.map(_.name).distinct.mkString(","))
            n += recs.size
            true
        }
      }

      if (isTechnical || allReads) recs.foreach(out.addAlignment)
      recs.foreach(progress.record)
    }

    logger.info(s"Found ${n} reads that matched to technical sequences.")
    in.safelyClose()
    out.close()
  }

  /** Returns true if all reads from the template are unmapped, otherwise false. */
  private def isUnmappedTemplate(rec: SAMRecord) : Boolean = rec.getReadUnmappedFlag && (!rec.getReadPairedFlag || rec.getMateUnmappedFlag)
}

private[bam] case class NamedSequence(name: String, bases: String)
private[bam] object NamedSequence {
  def apply(s: String): NamedSequence = s.indexOf(':') match {
    case -1 => NamedSequence(name=s, bases=s)
    case  n => NamedSequence(name=s.substring(0, n), bases=s.substring(n+1))
  }
}

/** Little utility class to match reads to technical sequences. */
private[bam] class Matcher(val matchLength: Int, val maxErrors: Int, sequences: Seq[NamedSequence]) {
  // Import the mutable version of a few collections that we need
  import scala.collection.mutable.{MultiMap, Set, LinkedHashMap => LinkedMap}

  // Map of kmer to the set of sequences in which it is found
  private val kmerToSequence = new LinkedMap[String,Set[NamedSequence]]() with MultiMap[String,NamedSequence]
  sequences.foreach(s => {
    Seq(s.bases.toUpperCase).flatMap(a => Seq(a, SequenceUtil.reverseComplement(a)))
      .flatMap(a => a.sliding(matchLength))
      .filterNot(_.contains('N'))
      .foreach(k => kmerToSequence.addBinding(k, s))
  })

  private val kmers: Array[Array[Byte]] = kmerToSequence.keySet.map(_.getBytes).toArray

  /** Just tests to see if any sequence matches. */
  def matches(read: Array[Byte]): Boolean = {
    read.length >= matchLength && kmers.exists(kmer => matches(read=read, kmer=kmer, len=matchLength))
  }

  /** Identifies the first sequence with a matching kmer. */
  def findMatches(read: Array[Byte]): Seq[NamedSequence] = {
    if (read.length < matchLength) {
      Nil
    }
    else {
      kmers.filter(kmer => matches(read=read, kmer=kmer, len=matchLength))
        .flatMap(kmer => kmerToSequence(new String(kmer))).toSeq.distinct
    }
  }

  @tailrec
  private def matches(read: Array[Byte], kmer: Array[Byte], index: Int = 0, len: Int, errors: Int = 0): Boolean = {
    if (errors > maxErrors) false
    else if (index >= len)  true
    else matches(read, kmer, index+1, len, errors + (if(read(index) == kmer(index)) 0 else 1))
  }
}
