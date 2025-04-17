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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric}
import enumeratum.EnumEntry
import htsjdk.samtools.SamPairUtil

@clp(group = ClpGroups.SamOrBam, description=
  """
    |Clips reads from the same template. Ensures that at least N bases are clipped from any end of the read (i.e.
    |R1 5' end, R1 3' end, R2 5' end, and R2 3' end).  Optionally clips reads from the same template to eliminate overlap
    |between the reads.  This ensures that downstream processes, particularly variant calling, cannot double-count
    |evidence from the same template when both reads span a variant site in the same template.
    |
    |Clipping overlapping reads is only performed on `FR` read pairs, and is implemented by clipping approximately half
    |the overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be substituted
    |using the `--soft-clip` parameter.
    |
    |Secondary alignments and supplemental alignments are not clipped, but are passed through into the
    |output.
    |
    |In order to correctly clip reads by template and update mate information, the input BAM must be either
    |`queryname` sorted or `query` grouped.  If your input BAM is not in an appropriate order the sort can be
    |done in streaming fashion with, for example:
    |
    |```
    |samtools sort -n -u in.bam | fgbio ClipBam -i /dev/stdin ...
    |```
    |
    |The output sort order may be specified with `--sort-order`.  If not given, then the output will be in the same
    |order as input.
    |
    |Any existing `NM`, `UQ` and `MD` tags are repaired, and mate-pair information updated.
    |
    |Three clipping modes are supported:
    |1. `Soft` - soft-clip the bases and qualities.
    |2. `SoftWithMask` - soft-clip and mask the bases and qualities (make bases Ns and qualities the minimum).
    |3. `Hard` - hard-clip the bases and qualities.
    |
    |The `--upgrade-clipping` parameter will convert all existing clipping in the input to the given more stringent mode:
    |from `Soft` to either `SoftWithMask` or `Hard`, and `SoftWithMask` to `Hard`. In all other cases, clipping remains
    |the same prior to applying any other clipping criteria.
  """)
class ClipBam
( @arg(flag='i', doc="Input SAM or BAM file of aligned reads in coordinate order.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='m', doc="Optional output of clipping metrics.") val metrics: Option[FilePath] = None,
  @arg(flag='r', doc="Reference sequence fasta file.") val ref: PathToFasta,
  @arg(flag='c', doc="The type of clipping to perform.") val clippingMode: ClippingMode = ClippingMode.Hard,
  @arg(flag='a', doc="Automatically clip extended attributes that are the same length as bases.") val autoClipAttributes: Boolean = false,
  @arg(flag='H', doc="Upgrade all existing clipping in the input to the given clipping mode prior to applying any other clipping criteria.") val upgradeClipping: Boolean = false,
  @arg(          doc="Require at least this number of bases to be clipped on the 5' end of R1") val readOneFivePrime: Int  = 0,
  @arg(          doc="Require at least this number of bases to be clipped on the 3' end of R1") val readOneThreePrime: Int = 0,
  @arg(          doc="Require at least this number of bases to be clipped on the 5' end of R2") val readTwoFivePrime: Int  = 0,
  @arg(          doc="Require at least this number of bases to be clipped on the 3' end of R2") val readTwoThreePrime: Int = 0,
  @arg(          doc="Clip overlapping reads.") val clipOverlappingReads: Boolean = false,
  @arg(          doc="Clip reads in FR pairs that sequence past the far end of their mate.") val clipBasesPastMate: Boolean = false,
  @arg(flag='S', doc="The sort order of the output. If not given, output will be in the same order as input if the input.")
  val sortOrder: Option[SamOrder] = None
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  validate(upgradeClipping || clipOverlappingReads || clipBasesPastMate || Seq(readOneFivePrime, readOneThreePrime, readTwoFivePrime, readTwoThreePrime).exists(_ != 0),
    "At least one clipping option is required")

  if (clipBasesPastMate && clipOverlappingReads) {
    logger.info("Clipping overlapping reads supersedes clipping past the far end of their mate.")
  }

  private val clipper = new SamRecordClipper(mode=clippingMode, autoClipAttributes=autoClipAttributes)

  override def execute(): Unit = {
    val in     = SamSource(input)
    val header = in.header
    val out    = Bams.nmUqMdTagRegeneratingWriter(writer=SamWriter(output, header.clone(), sort=sortOrder), ref=ref)

    // Require queryname sorted or query grouped
    Bams.requireQueryGrouped(header=in.header, toolName="ClipBam")

    val metricsMap: Map[ReadType, ClippingMetrics] = this.metrics.map { _ =>
      ReadType.values.map { readType => readType -> ClippingMetrics(read_type=readType) }.toMap
    }.getOrElse(Map.empty)

    // Go through and clip reads and fix their mate information
    Bams.templateIterator(in).foreach { template =>
      if (this.upgradeClipping) template.allReads.foreach { r => this.clipper.upgradeAllClipping(r) }

      (template.r1, template.r2) match {
        case (Some(r1), Some(r2)) =>
          clipPair(r1=r1, r2=r2, r1Metric=metricsMap.get(ReadType.ReadOne), r2Metric=metricsMap.get(ReadType.ReadTwo))
          SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
          template.r1Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s.asSam, r2.asSam, true))
          template.r2Supplementals.foreach(s => SamPairUtil.setMateInformationOnSupplementalAlignment(s.asSam, r1.asSam, true))
        case (Some(frag), None) =>
          clipFragment(frag=frag, metric=metricsMap.get(ReadType.Fragment))
        case _ => ()
      }

      out ++= template.allReads
    }
    out.close()

    this.metrics.foreach { path =>
      // Update the metrics for "All" and "Pair" read types
      import ReadType._
      metricsMap.foreach {
        case (Fragment, metric)          => metricsMap(All).add(metric)
        case (ReadOne | ReadTwo, metric) => Seq(Pair, All).foreach { r => metricsMap(r).add(metric) }
        case _                           => ()
      }
      // Write it!
      Metric.write(path, ReadType.values.map { readType => metricsMap(readType)})
    }

  }

  /** Clips a fixed amount from the reads and then clips overlapping reads.
    */
  private[bam] def clipFragment(frag: SamRecord, metric: Option[ClippingMetrics] = None): Unit = {
    val priorBasesClipped = frag.cigar.clippedBases

    // Clip the read!
    val numFivePrime  = this.clipper.clip5PrimeEndOfRead(frag, readOneFivePrime)
    val numThreePrime = this.clipper.clip3PrimeEndOfRead(frag, readOneThreePrime)

    // Update metrics
    metric.foreach { m =>
      m.update(
        rec                 = frag,
        priorBasesClipped   = priorBasesClipped,
        numFivePrime        = numFivePrime,
        numThreePrime       = numThreePrime,
        numOverlappingBases = 0,
        numExtendingBases   = 0
      )
    }
  }

  /** Clips a fixed amount from the reads and then clips overlapping reads.
    */
  private[bam] def clipPair(r1: SamRecord, r2: SamRecord, r1Metric: Option[ClippingMetrics] = None, r2Metric: Option[ClippingMetrics] = None): Unit = {
    val priorBasesClippedReadOne = r1.cigar.clippedBases
    val priorBasesClippedReadTwo = r2.cigar.clippedBases

    // Clip the read!
    val numReadOneFivePrime  = this.clipper.clip5PrimeEndOfRead(r1, readOneFivePrime)
    val numReadOneThreePrime = this.clipper.clip3PrimeEndOfRead(r1, readOneThreePrime)
    val numReadTwoFivePrime  = this.clipper.clip5PrimeEndOfRead(r2, readTwoFivePrime)
    val numReadTwoThreePrime = this.clipper.clip3PrimeEndOfRead(r2, readTwoThreePrime)

    val (numOverlappingBasesReadOne, numOverlappingBasesReadTwo) = {
      if (clipOverlappingReads && r1.isFrPair) this.clipper.clipOverlappingReads(r1, r2)
      else (0, 0)
    }

    val (numExtendingPastMateStartReadOne, numExtendingPastMateStartReadTwo) = {
      if (clipBasesPastMate && r1.isFrPair) {
        this.clipper.clipExtendingPastMateEnds(rec=r1, mate=r2)
      }
      else (0, 0)
    }

    r1Metric.foreach { m =>
      m.update(
        rec                 = r1,
        priorBasesClipped   = priorBasesClippedReadOne,
        numFivePrime        = numReadOneFivePrime,
        numThreePrime       = numReadOneThreePrime,
        numOverlappingBases = numOverlappingBasesReadOne,
        numExtendingBases   = numExtendingPastMateStartReadOne
      )
    }

    r2Metric.foreach { m =>
      m.update(
        rec                 = r2,
        priorBasesClipped   = priorBasesClippedReadTwo,
        numFivePrime        = numReadTwoFivePrime,
        numThreePrime       = numReadTwoThreePrime,
        numOverlappingBases = numOverlappingBasesReadTwo,
        numExtendingBases   = numExtendingPastMateStartReadTwo
      )
    }
  }
}

sealed trait ReadType extends EnumEntry
object ReadType extends FgBioEnum[ReadType] {
  def values: IndexedSeq[ReadType] = findValues
  case object Fragment extends ReadType
  case object ReadOne extends ReadType
  case object ReadTwo extends ReadType
  case object Pair extends ReadType
  case object All extends ReadType
}


/** Metrics produced by [[ClipBam]] that detail how many reads and bases are clipped respectively.
  *
  * @param read_type The type of read (i.e. Fragment, ReadOne, ReadTwo).
  * @param reads The number of reads examined.
  * @param reads_clipped_pre The number of reads with any type of clipping prior to clipping with [[ClipBam]].
  * @param reads_clipped_post The number of reads with any type of clipping after clipping with [[ClipBam]], including reads that became unmapped.
  * @param reads_clipped_five_prime The number of reads with the 5' end clipped.
  * @param reads_clipped_three_prime The number of reads with the 3' end clipped.
  * @param reads_clipped_overlapping The number of reads clipped due to overlapping reads.
  * @param reads_clipped_extending The number of reads clipped due to a read extending past its mate.
  * @param reads_unmapped The number of reads that became unmapped due to clipping.
  * @param bases The number of aligned bases after clipping.
  * @param bases_clipped_pre The number of bases clipped prior to clipping with [[ClipBam]].
  * @param bases_clipped_post The number of bases clipped after clipping with [[ClipBam]], including bases from reads that became unmapped.
  * @param bases_clipped_five_prime The number of bases clipped on the 5' end of the read.
  * @param bases_clipped_three_prime The number of bases clipped on the 3 end of the read.
  * @param bases_clipped_overlapping The number of bases clipped due to overlapping reads.
  * @param bases_clipped_extending The number of bases clipped due to a read extending past its mate.
  */
case class ClippingMetrics
(read_type: ReadType,
 var reads: Long = 0,
 var reads_unmapped: Long = 0,
 var reads_clipped_pre: Long = 0,
 var reads_clipped_post: Long = 0,
 var reads_clipped_five_prime: Long = 0,
 var reads_clipped_three_prime: Long = 0,
 var reads_clipped_overlapping: Long = 0,
 var reads_clipped_extending: Long = 0,
 var bases: Long = 0,
 var bases_clipped_pre: Long = 0,
 var bases_clipped_post: Long = 0,
 var bases_clipped_five_prime: Long = 0,
 var bases_clipped_three_prime: Long = 0,
 var bases_clipped_overlapping: Long = 0,
 var bases_clipped_extending: Long = 0,
) extends Metric {
  def update(rec: SamRecord, priorBasesClipped: Int, numFivePrime: Int, numThreePrime: Int, numOverlappingBases: Int, numExtendingBases: Int): Unit = {
    this.reads                      += 1
    this.bases                      += rec.cigar.alignedBases
    if (priorBasesClipped > 0) {
      this.reads_clipped_pre        += 1
      this.bases_clipped_pre        += priorBasesClipped
    }
    if (numFivePrime > 0) {
      this.reads_clipped_five_prime += 1
      this.bases_clipped_five_prime += numFivePrime
    }
    if (numThreePrime > 0) {
      this.reads_clipped_three_prime += 1
      this.bases_clipped_three_prime += numThreePrime
    }
    if (numOverlappingBases > 0) {
      this.reads_clipped_overlapping += 1
      this.bases_clipped_overlapping += numOverlappingBases
    }
    if (numExtendingBases > 0) {
      this.reads_clipped_extending += 1
      this.bases_clipped_extending += numExtendingBases
    }
    val additionalClippedBases = numFivePrime + numThreePrime + numOverlappingBases + numExtendingBases
    val totalClippedBases = additionalClippedBases + priorBasesClipped
    if (totalClippedBases > 0) {
      this.reads_clipped_post        += 1
      this.bases_clipped_post        += totalClippedBases
      if (rec.unmapped && additionalClippedBases > 0) this.reads_unmapped += 1
    }
  }
  
  def add(metric: ClippingMetrics*): Unit = {
    this.reads                     += metric.sumBy(_.reads)
    this.reads_unmapped            += metric.sumBy(_.reads_unmapped)
    this.reads_clipped_pre         += metric.sumBy(_.reads_clipped_pre)
    this.reads_clipped_post        += metric.sumBy(_.reads_clipped_post)
    this.reads_clipped_five_prime  += metric.sumBy(_.reads_clipped_five_prime)
    this.reads_clipped_three_prime += metric.sumBy(_.reads_clipped_three_prime)
    this.reads_clipped_overlapping += metric.sumBy(_.reads_clipped_overlapping)
    this.reads_clipped_extending   += metric.sumBy(_.reads_clipped_extending)
    this.bases                     += metric.sumBy(_.bases)
    this.bases_clipped_pre         += metric.sumBy(_.bases_clipped_pre)
    this.bases_clipped_post        += metric.sumBy(_.bases_clipped_post)
    this.bases_clipped_five_prime  += metric.sumBy(_.bases_clipped_five_prime)
    this.bases_clipped_three_prime += metric.sumBy(_.bases_clipped_three_prime)
    this.bases_clipped_overlapping += metric.sumBy(_.bases_clipped_overlapping)
    this.bases_clipped_extending   += metric.sumBy(_.bases_clipped_extending)
  }
}
