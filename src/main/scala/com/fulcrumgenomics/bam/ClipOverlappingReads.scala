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
import com.fulcrumgenomics.bam.SamRecordClipper.ClippingMode
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker

import scala.collection.JavaConversions.iterableAsScalaIterable

@clp(group = ClpGroups.SamOrBam, description=
  """
    |Clips reads from the same template to eliminate overlap between the reads. Ensures that downstream
    |processes, particularly variant calling, cannot double-count evidence from the same template when
    |both reads span a variant.
    |
    |Clipping is only performed on FR read pairs, and is implemented by clipping approximately half the
    |overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be
    |substituted using the --soft-clip parameter.
    |
    |Secondary alignments and supplemental alignments are not clipped, but are passed through into the
    |output.
    |
    |If the input BAM is neither queryname sorted nor query-grouped, it will be sorted into queryname
    |order so that clipping can be performed on both ends of a pair simultaneously and so that mate
    |pair information can be reset across all reads for the template.  Post-clipping the reads are
    |resorted into coordinate order, any existing NM, UQ and MD tags are repaired, and the output is
    |written in coordinate order.
  """)
class ClipOverlappingReads
( @arg(flag="i", doc="Input SAM or BAM file of aligned reads in coordinate order.") val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag="s", doc="Soft clip reads instead of hard clipping.") val softClip: Boolean = false,
  @arg(flag="r", doc="Reference sequence fasta file.") val ref: PathToFasta
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  private val clippingMode = if (softClip) ClippingMode.Soft else ClippingMode.Hard

  override def execute(): Unit = {
    val in       = SamReaderFactory.make().open(input)
    val progress = new ProgressLogger(logger)
    val sorter   = Bams.sortingCollection(SortOrder.coordinate, in.getFileHeader)

    // Go through and clip reads and fix their mate information
    Bams.templateIterator(in).foreach { template =>
      (template.r1, template.r2) match {
        case (Some(r1), Some(r2)) =>
          clip(r1, r2)
          SamPairUtil.setMateInfo(r1, r2, true)
          template.r1Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s, r2, true))
          template.r2Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s, r1, true))
        case _ => Unit
      }

      // TODO: add an allReads iterator to Template
      (template.r1.iterator ++ template.r2.iterator ++ template.allSupplementaryAndSecondary).foreach { r =>
        sorter.add(r)
        progress.record(r)
      }
    }

    // Then go through the coordinate sorted reads and fix up tags
    logger.info("Re-sorting into coordinate order and writing output.")
    val header = in.getFileHeader.clone()
    header.setSortOrder(SortOrder.coordinate)
    val walker = new ReferenceSequenceFileWalker(ref.toFile)
    val out    = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(header, true, output.toFile, ref.toFile)
    val writeProgress = new ProgressLogger(logger, verb="wrote")

    sorter.foreach { rec =>
      Bams.regenerateNmUqMdTags(rec, walker)
      out.addAlignment(rec)
      writeProgress.record(rec)
    }

    out.close()
  }

  /** Returns true if the read is from a paired end insert with both reads mapped to the same
    * chromosome in FR orientation.  Does not check that the reads actually overlap!
    */
  private[bam] def clip(r1: SAMRecord, r2: SAMRecord): Unit = {
    if (r1.getReadUnmappedFlag || r2.getReadUnmappedFlag || r1.getReferenceIndex != r2.getReferenceIndex ||
      SamPairUtil.getPairOrientation(r1) != PairOrientation.FR) {
      // Do nothing
    }
    else {
      val (f,r) = if (r1.getReadNegativeStrandFlag) (r2, r1) else (r1, r2)

      // What we really want is to trim by the number of _reference_ bases not read bases,
      // in order to eliminate overlap.  We could do something very complicated here, or
      // we could just trim read bases in a loop until the overlap is eliminated!
      while (f.getAlignmentEnd >= r.getAlignmentStart && !f.getReadUnmappedFlag && !r.getReadUnmappedFlag) {
        val lengthToClip = f.getAlignmentEnd - r.getAlignmentStart + 1
        val firstHalf    = lengthToClip / 2
        val secondHalf   = lengthToClip - firstHalf // safe guard against rounding on odd lengths
        SamRecordClipper.clip3PrimeEndOfAlignment(r1, firstHalf,  mode=this.clippingMode)
        SamRecordClipper.clip3PrimeEndOfAlignment(r2, secondHalf, mode=this.clippingMode)
      }
    }
  }
}
