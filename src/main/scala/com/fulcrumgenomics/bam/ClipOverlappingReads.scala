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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker

@clp(group = ClpGroups.SamOrBam, description=
  """
    |Clips reads from the same template to eliminate overlap between the reads. Ensures that downstream
    |processes, particularly variant calling, cannot double-count evidence from the same template when
    |both reads span a variant site in the same template.
    |
    |Clipping is only performed on `FR` read pairs, and is implemented by clipping approximately half the
    |overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be
    |substituted using the `--soft-clip` parameter.
    |
    |Secondary alignments and supplemental alignments are not clipped, but are passed through into the
    |output.
    |
    |If the input BAM is neither `queryname` sorted nor `query` grouped, it will be sorted into queryname
    |order so that clipping can be performed on both ends of a pair simultaneously and so that mate
    |pair information can be reset across all reads for the template.  Post-clipping the reads are
    |resorted into coordinate order, any existing `NM`, `UQ` and `MD` tags are repaired, and the output is
    |written in coordinate order.
  """)
class ClipOverlappingReads
( @arg(flag='i', doc="Input SAM or BAM file of aligned reads in coordinate order.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='s', doc="Soft clip reads instead of hard clipping.") val softClip: Boolean = false,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: PathToFasta,
  @arg(flag='a', doc="Automatically clip extended attributes that are the same length as bases.") val autoClipAttributes: Boolean = false
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  private val clipper = new SamRecordClipper(mode=if (softClip) ClippingMode.Soft else ClippingMode.Hard, autoClipAttributes=autoClipAttributes)

  override def execute(): Unit = {
    val in       = SamSource(input)
    val progress = ProgressLogger(logger)
    val sorter   = Bams.sorter(SamOrder.Coordinate, in.header)

    // Go through and clip reads and fix their mate information
    Bams.templateIterator(in).foreach { template =>
      (template.r1, template.r2) match {
        case (Some(r1), Some(r2)) =>
          clip(r1, r2)
          SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
          template.r1Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s.asSam, r2.asSam, true))
          template.r2Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s.asSam, r1.asSam, true))
        case _ => Unit
      }

      template.allReads.foreach { r =>
        sorter += r
        progress.record(r)
      }
    }

    // Then go through the coordinate sorted reads and fix up tags
    logger.info("Re-sorting into coordinate order and writing output.")
    val header = in.header.clone()
    SamOrder.Coordinate.applyTo(header)
    header.setSortOrder(SortOrder.coordinate)
    val walker = new ReferenceSequenceFileWalker(ref.toFile)
    val out    = SamWriter(output, header, ref=Some(ref))

    sorter.foreach { rec =>
      Bams.regenerateNmUqMdTags(rec, walker)
      out += rec
    }

    out.close()
  }

  /** Returns true if the read is from a paired end insert with both reads mapped to the same
    * chromosome in FR orientation.  Does not check that the reads actually overlap!
    */
  private[bam] def clip(r1: SamRecord, r2: SamRecord): Unit = {
    if (r1.unmapped || r2.unmapped || r1.refIndex != r2.refIndex || r1.pairOrientation != PairOrientation.FR) {
      // Do nothing
    }
    else {
      val (f,r) = if (r1.negativeStrand) (r2, r1) else (r1, r2)

      // What we really want is to trim by the number of _reference_ bases not read bases,
      // in order to eliminate overlap.  We could do something very complicated here, or
      // we could just trim read bases in a loop until the overlap is eliminated!
      while (f.end >= r.start && f.mapped && r.mapped) {
        val lengthToClip = f.end - r.start + 1
        val firstHalf    = lengthToClip / 2
        val secondHalf   = lengthToClip - firstHalf // safe guard against rounding on odd lengths
        this.clipper.clip3PrimeEndOfAlignment(r1, firstHalf)
        this.clipper.clip3PrimeEndOfAlignment(r2, secondHalf)
      }
    }
  }
}
