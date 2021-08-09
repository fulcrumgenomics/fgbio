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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Amplicon, AmpliconDetector, Io, ProgressLogger}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.reference.ReferenceSequenceFileWalker
import htsjdk.samtools.util._

import scala.collection.BufferedIterator

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Trims primers from reads post-alignment.  Takes in a BAM file of aligned reads
    |and a tab-delimited file with five columns (`chrom`, `left_start`, `left_end`,
    |`right_start`, and `right_end`) which provide the 1-based inclusive start and end
    |positions of the primers for each amplicon.  The primer file must include headers, e.g:
    |
    |```
    |chrom  left_start  left_end  right_start right_end
    |chr1   1010873     1010894   1011118     1011137
    |```
    |
    |Paired end reads that map to a given amplicon position are trimmed so that the
    |alignment no-longer includes the primer sequences. All other aligned reads have the
    |_maximum primer length trimmed_!
    |
    |Reads that are trimmed will have the `NM`, `UQ` and `MD` tags cleared as they are no longer
    |guaranteed to be accurate.  If a reference is provided the reads will be re-sorted
    |by coordinate after trimming and the `NM`, `UQ` and `MD` tags recalculated.
    |
    |If the input BAM is not `queryname` sorted it will be sorted internally so that mate
    |information between paired-end reads can be corrected before writing the output file.
    |
    |The `--first-of-pair` option will cause only the first of pair (R1) reads to be trimmed
    |based solely on the primer location of R1.  This is useful when there is a target
    |specific primer on the 5' end of R1 but no primer sequenced on R2 (eg. single gene-specific
    |primer target enrichment).  In this case, the location of each target specific primer should
    |be specified in an amplicons left or right primer exclusively.  The coordinates of the
    |non-specific-target primer should be `-1` for both start and end, e.g:
    |
    |```
    |chrom  left_start  left_end  right_start right_end
    |chr1   1010873     1010894   -1          -1
    |chr2   -1          -1        1011118     1011137
    |```
  """)
class TrimPrimers
( @arg(flag='i', doc="Input BAM file.")  val input: PathToBam,
  @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
  @arg(flag='p', doc="File with primer locations.") val primers: FilePath,
  @arg(flag='H', doc="If true, hard clip reads, else soft clip.") val hardClip: Boolean = false,
  @arg(flag='S', doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
  @arg(flag='s', doc="Sort order of output BAM file (defaults to input sort order).") val sortOrder: Option[SamOrder] = None,
  @arg(flag='r', doc="Optional reference fasta for recalculating NM, MD and UQ tags.") val ref: Option[PathToFasta] = None,
  @arg(flag='a', doc="Automatically trim extended attributes that are the same length as bases.") val autoTrimAttributes: Boolean = false,
  @arg(doc="Trim only first of pair reads (R1s), otherwise both ends of a pair") val firstOfPair: Boolean = false

)extends FgBioTool with LazyLogging {
  private val clipper = new SamRecordClipper(mode=if (hardClip) ClippingMode.Hard else ClippingMode.Soft, autoClipAttributes=autoTrimAttributes)

  Io.assertReadable(input)
  Io.assertReadable(primers)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in = SamSource(input)
    val outSortOrder = sortOrder.orElse(SamOrder(in.header)).getOrElse(SamOrder.Unknown)
    val outHeader = in.header.clone()
    outSortOrder.applyTo(outHeader)

    // Setup the outputs depending on whether or not we have a reference file or not
    // In order to minimize (ha!) the amount of sorting going on, things are a little complex.
    // The logic is more or less as follows:
    //   a) For trimming we need things in queryname order, so if input isn't queryname sorted, sort it
    //   b) Then, if we were given a reference, we need to coordinate sort to reset the MD, NM, UQ tags
    //   c) Then finally for output, the SAMFileWriter will see things as pre-sorted if:
    //         i) No reference was given (so no coordinate sort) and the output sort order is queryname, or
    //        ii) A reference was given (so the last sort is coordinate) and the output sort order is coordinate
    //      else, we'll have to sort for potentially a third time in the SAMFileWriter
    val (sorter, write: (SamRecord => Any), out: SamWriter) = ref match {
      case Some(_) =>
        val sorter = Bams.sorter(SamOrder.Coordinate, outHeader)
        val out = SamWriter(output, outHeader, sort= if (outSortOrder == SamOrder.Coordinate) None else Some(outSortOrder))
        (Some(sorter), sorter.write _, out)
      case None =>
        val out = SamWriter(output, outHeader, sort= if (outSortOrder == SamOrder.Queryname) None else Some(outSortOrder))
        (None, out += _, out)
    }

    val detector = new AmpliconDetector(
      detector             = Amplicon.overlapDetector(path=primers),
      slop                 = slop,
      unclippedCoordinates = true
    )

    // Validate that if we are trimming the first of pair, all amplicons have -1 coordinates for either the left or the
    // right primer.  Otherwise, coordinates must be > 0.
    detector.detector.getAll.iterator().foreach { amplicon =>
      if (firstOfPair) {
        require(amplicon.leftStart == -1 || amplicon.rightStart == -1,
          f"Either the left or right amplicon coordinates must be -1 when using --first-of-pair. Found ${amplicon.mkString("\t")}"
        )
      }
      else {
        require(amplicon.leftStart > 0 && amplicon.rightStart > 0,
          f"Both the left and right amplicon coordinates must be > 0. Did you forget to set '--first-of-pair'? Found ${amplicon.mkString("\t")}"
        )
      }
    }

    // Main processing loop
    val iterator = queryNameOrderIterator(in)
    val trimProgress = ProgressLogger(this.logger, "Trimmed")
    while (iterator.hasNext) {
      val reads = nextTemplate(iterator)
      trimReadsForTemplate(detector, reads)
      reads.foreach(write)
      reads.foreach(trimProgress.record)
    }

    // If we had a reference and re-sorted above, reset the NM/UQ/MD tags as we push to the final output
    (sorter, ref) match {
      case (Some(sorter), Some(path)) =>
        val walker = new ReferenceSequenceFileWalker(path.toFile)
        val progress = ProgressLogger(this.logger, "Written")

        sorter.foreach { rec =>
          recalculateTags(rec, walker)
          out += rec
          progress.record(rec)
        }

        sorter.safelyClose()
      case _ => ()
    }

    out.close()
  }

  /** Recalculates the MD, NM, and UQ tags on aligned records. */
  def recalculateTags(rec: SamRecord, walker: ReferenceSequenceFileWalker): Unit = {
    if (rec.mapped) {
      val refBases = walker.get(rec.refIndex).getBases
      SequenceUtil.calculateMdAndNmTags(rec.asSam, refBases, true, true)
      if (rec.quals != null && rec.quals.length != 0) {
        rec(SAMTag.UQ.name) = SequenceUtil.sumQualitiesOfMismatches(rec.asSam, refBases, 0)
      }
    }
  }

  /** Trims all the reads for a given template. */
  def trimReadsForTemplate(detector: AmpliconDetector, reads: Seq[SamRecord]): Unit = {
    val rec1 = reads.find(r => r.paired && r.firstOfPair  && !r.secondary && !r.supplementary)
    val rec2 = reads.find(r => r.paired && r.secondOfPair && !r.secondary && !r.supplementary)

    val readsToClip = if (firstOfPair) reads.filter(_.firstOfPair) else reads

    (rec1, rec2) match {
      case (Some(r1), Some(r2)) =>
        // FR mapped pairs get the full treatment
        if (r1.isFrPair) {
          (if (firstOfPair) detector.findPrimer(rec=r1) else detector.find(r1=r1, r2=r2)) match {
            case Some(amplicon) =>
              val leftClip  = amplicon.leftPrimerLength
              val rightClip = amplicon.rightPrimerLength
              readsToClip.foreach { rec =>
                // Note: r1.positiveStrand means that r1 is the "left" read, so we should clip on the left
                val toClip = if (rec.firstOfPair == r1.positiveStrand) leftClip else rightClip
                this.clipper.clip5PrimeEndOfRead(rec, toClip)
              }
            case None =>
              readsToClip.foreach(r => this.clipper.clip5PrimeEndOfRead(r, detector.maxPrimerLength))
          }

          clipFullyOverlappedFrReads(r1, r2)
        }
        // Pairs without both reads mapped in FR orientation are just maximally clipped
        else {
          readsToClip.foreach(r => this.clipper.clip5PrimeEndOfRead(r, detector.maxPrimerLength))
        }

        SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
        reads.filter(_.supplementary).foreach { rec =>
          val mate = if (rec.firstOfPair) r2 else r1
          SamPairUtil.setMateInformationOnSupplementalAlignment(rec.asSam, mate.asSam, true);
        }
      case _ =>
        // Just trim each read independently
        readsToClip.foreach(r => this.clipper.clip5PrimeEndOfRead(r, detector.maxPrimerLength))
    }
  }

  /** Gets an iterator in query name order over the records. */
  private def queryNameOrderIterator(in: SamSource): BufferedIterator[SamRecord] = {
    if (in.header.getSortOrder == SortOrder.queryname) {
      in.iterator.bufferBetter
    }
    else {
      logger.info("Sorting into queryname order.")
      val progress = ProgressLogger(this.logger, "Queryname sorted")
      val sorter = Bams.sorter(SamOrder.Queryname, in.header)
      in.foreach { rec =>
        sorter += rec
        progress.record(rec)
      }
      sorter.iterator
    }
  }

  /** Fetches the next group of records that all share the same readname/template from the iterator. */
  private def nextTemplate(iterator: BufferedIterator[SamRecord]): Seq[SamRecord] = {
    val first    = iterator.next()
    val template = first.name
    first :: iterator.takeWhile(_.name == template).toList
  }

  /**
    * Adapted from Picard's AbstractAlignmentMerger - takes in mapped FR pairs only
    * and clips from the 3' ends of the reads if the reads are fully overlapped
    * and extend past each other's starts.
    */
  private def clipFullyOverlappedFrReads(r1: SamRecord, r2: SamRecord): Unit = {
    val (plus, minus) = if (r1.negativeStrand) (r2,r1) else (r1, r2)

    if (plus.start < minus.end) {
      val plusTrim  = plus.end   - minus.end
      val minusTrim = plus.start - minus.start

      if (plusTrim  > 0) this.clipper.clip3PrimeEndOfAlignment(plus, plusTrim)
      if (minusTrim > 0) this.clipper.clip3PrimeEndOfAlignment(minus, minusTrim)
    }
  }
}
