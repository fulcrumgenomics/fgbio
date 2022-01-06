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
 *
 */

package com.fulcrumgenomics.rnaseq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.PathToBam
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.GeneAnnotations._
import com.fulcrumgenomics.util._
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.filter._
import htsjdk.samtools.util.{CoordMath, Interval, OverlapDetector}

import scala.collection.JavaConverters._

@clp(description =
  """Computes the insert size for RNA-Seq experiments.
    |
    |Computes the insert size by counting the # of bases sequenced in transcript space.  The insert size is defined
    |as the distance between the first bases sequenced of each pair respectively (5' sequencing ends).
    |
    |This tool skips reads that overlap multiple genes, reads that aren't fully enclosed in a gene, and reads where the
    |insert size would disagree across transcripts from the same gene.  Also skips reads that are unpaired, failed QC,
    |secondary, supplementary, pairs without both ends mapped, duplicates, and pairs whose reads map to different
    |chromosomes. Finally, skips transcripts where too few mapped read bases overlap exonic sequence.
    |
    |This tool requires each mapped pair to have the mate cigar (`MC`) tag.  Use `SetMateInformation` to add the mate cigar.
    |
    |The output metric file will have the extension `.rna_seq_insert_size.txt` and the output histogram file will have
    |the extension `.rna_seq_insert_size_histogram.txt`.  The histogram file gives for each orientation (`FR`, `RF`, `tandem`),
    |the number of read pairs that had the given insert size.
  """,
  group = ClpGroups.RnaSeq)
class EstimateRnaSeqInsertSize
(@arg(flag='i', doc="Input BAM file.") val input: PathToBam,
 @arg(flag='r', doc="Input gene annotations in [RefFlat](http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat) form")
 val refFlat: FilePath,
 @arg(flag='p', doc="Output prefix file.  The file will have the extension `.rna_seq_insert_size.txt` if not given")
 val prefix: Option[PathPrefix] = None,
 @arg(flag='d', doc="Include duplicates") val includeDuplicates: Boolean = false,
 @arg(flag='D', doc=
   """Generate mean and standard deviation by filtering to `median + deviations*median_absolute_deviation`.
      |This is done because insert size data typically includes enough anomalous values from chimeras
      |and other artifacts to make the mean and sd grossly misleading regarding the real distribution.
  """"
 ) val deviations: Double = 10.0,
 @arg(flag='q', doc="Ignore reads with mapping quality less than this value.") val minimumMappingQuality: Int = 30,
 @arg(flag='m', doc="The minimum fraction of read bases that must overlap exonic sequence in a transcript") val minimumOverlap: Double = 0.95
) extends FgBioTool with LazyLogging {
  import EstimateRnaSeqInsertSize._

  Io.assertReadable(input)
  Io.assertReadable(refFlat)
  prefix.foreach(Io.assertCanWriteFile(_))

  private val geneOverlapDetector = new OverlapDetector[GeneLocus](0, 0)

  override def execute(): Unit = {
    val progress            = ProgressLogger(logger, verb = "read", unit = 5e6.toInt)
    val pairOrientations    = PairOrientation.values()
    val in                  = SamSource(input)
    val refFlatSource       = RefFlatSource(refFlat, Some(in.dict))
    val counters            = pairOrientations.map { pairOrientation => (pairOrientation, new NumericCounter[Long]()) }.toMap
    val filter              = new AggregateFilter(EstimateRnaSeqInsertSize.filters(minimumMappingQuality=minimumMappingQuality, includeDuplicates=includeDuplicates).asJava)
    var numReadPairs        = 0L
    val recordIterator      = in.iterator.filter { rec =>
      progress.record(rec)
      if (rec.paired && rec.firstOfPair) numReadPairs += 1
      !filter.filterOut(rec.asSam)
    }

    for (gene <- refFlatSource; locus <- gene.loci) {
      geneOverlapDetector.addLhs(locus, locus)
    }

    recordIterator.foreach { rec =>
      calculateInsertSize(rec=rec) match {
        case None             => () // ignore
        case Some(insertSize) => counters(rec.pairOrientation).count(insertSize)
      }
    }

    val examinedPairs = pairOrientations.map(counters(_).total).sum
    logger.info(f"Used $examinedPairs out of $numReadPairs (${examinedPairs/numReadPairs.toDouble * 100}%.2f%%) read pair(s) to estimate the insert size.")

    refFlatSource.safelyClose()
    in.safelyClose()

    // Collate the metrics
    val metrics = pairOrientations.map { pairOrientation =>
      val counter = counters(pairOrientation)

      // Create a metric *without* mean and standard deviation, as we will filter after
      val median  = counter.median()
      val mad     = counter.mad(m=median)
      val metric = InsertSizeMetric(
        pair_orientation           = pairOrientation,
        read_pairs                = counter.total,
        median                    = median,
        min                       = if (counter.nonEmpty) counter.map(_._2).min else 0,
        max                       = if (counter.nonEmpty) counter.map(_._2).max else 0,
        median_absolute_deviation = mad
      )

      // Remove outliers based on the mean absolute deviation and the # of allowable deviations
      val maximumInsertSize = median + (deviations * mad)
      val filteredCounter   = NumericCounter.from(counter.filter { case (value, _) => value <= maximumInsertSize })
      val filteredMean      = filteredCounter.mean()

      metric.copy(mean=filteredMean, standard_deviation=filteredCounter.stddev(m=filteredMean))
    }

    // Write the metrics
    val actualPrefix = prefix.getOrElse(PathUtil.removeExtension(input))
    val metricPath = PathUtil.pathTo(s"${actualPrefix}${RnaSeqInsertSizeMetricExtension}")
    Metric.write(metricPath, metrics)

    // Write the histogram
    val histogramPath   = PathUtil.pathTo(s"${actualPrefix}${RnaSeqInsertSizeMetricHistogramExtension}")
    val histogramKeys   = pairOrientations.flatMap(counters(_).iterator.map(_._1)).distinct.sorted.toList
    val histogramHeader = ("insert_size" +: pairOrientations.map(_.name().toLowerCase)).mkString("\t")
    val histogramLines  = histogramHeader +: histogramKeys.map { insertSize =>
      (insertSize +: pairOrientations.map(counters(_).countOf(insertSize))).mkString("\t")
    }
    Io.writeLines(path=histogramPath, lines=histogramLines)
  }

  /** Calculates the insert size in transcript space. */
  private def calculateInsertSize(rec: SamRecord): Option[Int] = {
    // Get the overlapping genes
    val mateCigar        = getAndRequireMateCigar(rec)
    val mateAlignmentEnd = mateAlignmentEndFrom(mateCigar, rec.mateStart)
    val recInterval      = intervalFrom(rec, mateAlignmentEnd=mateAlignmentEnd)
    val overlappingGenes = geneOverlapDetector.getOverlaps(recInterval)

    if (overlappingGenes.size() == 1) {
      insertSizeFromGene(
        rec              = rec,
        gene             = overlappingGenes.iterator().next(),
        minimumOverlap   = minimumOverlap,
        recInterval      = recInterval,
        recBlocks        = rec.asSam.getAlignmentBlocks.toList,
        mateBlocks       = mateAlignmentBlocksFrom(mateCigar, rec.mateStart),
        mateAlignmentEnd = mateAlignmentEnd
      )
    }
    else None
  }
}

object EstimateRnaSeqInsertSize {
  val RnaSeqInsertSizeMetricExtension: String = ".rnaseq_insert_size.txt"
  val RnaSeqInsertSizeMetricHistogramExtension: String = ".rnaseq_insert_size_histogram.txt"

  private trait SingleEndSamRecordFilter extends SamRecordFilter {
    override def filterOut(r1: SAMRecord, r2: SAMRecord): Boolean = filterOut(r1) || filterOut(r2)
  }

  def filter(f: SamRecord => Boolean): SamRecordFilter = new SingleEndSamRecordFilter {
    override def filterOut(r: SAMRecord): Boolean = f(r.asInstanceOf[SamRecord])
  }

  private def filters(minimumMappingQuality: Int, includeDuplicates: Boolean) = List(
    MateMappedFilter, // also filters out unpaired reads
    new FailsVendorReadQualityFilter,
    new AlignedFilter(true),
    new SecondaryOrSupplementaryFilter,
    new MappingQualityFilter(minimumMappingQuality),
    DuplicatesFilter(includeDuplicates=includeDuplicates),
    FirstOfPairOnlyFilter,
    DifferentReferenceIndexFilter
  )

  private def MateMappedFilter                             = filter(r => !r.paired || r.mateUnmapped)
  private def DuplicatesFilter(includeDuplicates: Boolean) = filter(r => !includeDuplicates && r.duplicate)
  private def FirstOfPairOnlyFilter                        = filter(_.firstOfPair)
  private def DifferentReferenceIndexFilter                = filter(r => r.refIndex != r.mateRefIndex)

  private[rnaseq] def getAndRequireMateCigar(rec: SamRecord): Cigar = {
    rec.mateCigar.getOrElse {
      throw new IllegalStateException(s"Mate CIGAR (Tag 'MC') not found for $rec, consider using SetMateInformation.")
    }
  }

  /** Calculates an interval representing the records span.  If mateAlignmentEnd is not given, computes from
    * the mate cigar.
    */
  private[rnaseq] def intervalFrom(rec: SamRecord, mateAlignmentEnd: Int): Interval = {
    val leftMostBase     = Math.min(rec.start, rec.mateStart)
    val rightMostBase    = Math.max(rec.end, mateAlignmentEnd)
    new Interval(rec.refName, leftMostBase, rightMostBase)
  }

  /** Gets the mate's alignment end. */
  private[rnaseq] def mateAlignmentEndFrom(mateCigar: Cigar, mateAlignmentStart: Int): Int = {
    CoordMath.getEnd(mateAlignmentStart, mateCigar.lengthOnTarget)
  }

  /** Gets the mate's alignment blocks. */
  private[rnaseq] def mateAlignmentBlocksFrom(mateCigar: Cigar, mateAlignmentStart: Int): List[AlignmentBlock] = {
    SAMUtils.getAlignmentBlocks(mateCigar.toHtsjdkCigar, mateAlignmentStart, "mate cigar").toList
  }

  /** Calculates the insert size from a gene.  Returns None if the record's span is not enclosed in the gene or if
    * the insert size disagree across transcripts.  Assumes the record and gene are mapped to the same chromosome. */
  private[rnaseq] def insertSizeFromGene(rec: SamRecord,
                                         gene: GeneLocus,
                                         minimumOverlap: Double,
                                         recInterval: Interval,
                                         recBlocks: List[AlignmentBlock],
                                         mateBlocks: List[AlignmentBlock],
                                         mateAlignmentEnd: Int): Option[Int] = {
    if (!CoordMath.encloses(gene.start, gene.end, recInterval.getStart, recInterval.getEnd)) return None

    // Check the insert size of each transcript, making sure they all have the same value
    val transcripts             = gene.iterator
    var insertSize: Option[Int] = None
    var sameInsertSize          = true
    val mappedBases             = recBlocks.map(_.getLength).sum + mateBlocks.map(_.getLength).sum
    while (transcripts.hasNext && sameInsertSize) {
      val transcript = transcripts.next()

      // Ignore transcripts with too little overlap
      val overlap    = Seq(recBlocks, mateBlocks).map(numReadBasesOverlappingTranscript(_, transcript)).sum
      if (overlap / mappedBases.toDouble >= minimumOverlap) {
        // calculate the distance in transcript space
        val transcriptInsertSize = insertSizeFromTranscript(rec=rec, transcript=transcript, mateAlignmentEnd=mateAlignmentEnd)

        // Make sure that the insert size is the same as any previous result
        insertSize match {
          case None        => insertSize = Some(transcriptInsertSize)
          case Some(isize) => sameInsertSize = isize == transcriptInsertSize
        }
      }
    }

    if (sameInsertSize) insertSize
    else None
  }

  /** Computes the insert size (5' to 5') from a record and transcript, assuming that the record overlaps the transcript. */
  private[rnaseq] def insertSizeFromTranscript(rec: SamRecord, transcript: Transcript, mateAlignmentEnd: Int): Int = {
    // get the 5' position of the record and its mate
    val rec5Prime   = if (rec.negativeStrand)     rec.end          else rec.start
    val mate5Prime  = if (rec.mateNegativeStrand) mateAlignmentEnd else rec.mateStart

    // calculate the distance in transcript space
    val lower5Prime      = Math.min(rec5Prime, mate5Prime)
    val upper5Prime      = Math.max(rec5Prime, mate5Prime)
    val overlappingExons = transcript.genomicOrder.filter { exon => CoordMath.overlaps(lower5Prime, upper5Prime, exon.start, exon.end) }
    if (overlappingExons.isEmpty) 0
    else overlappingExons.map { exon => CoordMath.getOverlap(lower5Prime, upper5Prime, exon.start, exon.end) }.sum
  }

  /** Gets the number of read bases that overlap the exons.  This assumes the read and transcript are on the same chromosome. */
  private[rnaseq] def numReadBasesOverlappingTranscript(alignmentBlocks: List[AlignmentBlock], transcript: Transcript): Int = {
    val blocks  = alignmentBlocks.iterator.bufferBetter
    val exons   = transcript.genomicOrder.bufferBetter
    var overlap = 0

    while (blocks.hasNext && exons.hasNext) {
      val block               = blocks.head
      val exon                = exons.head
      val blockReferenceStart = block.getReferenceStart
      val blockReferenceEnd   = blockReferenceStart + block.getLength - 1
      val exonReferenceStart  = exon.start
      val exonReferenceEnd    = exon.end

      // Calculate the number of read bases in alignment block overlap the exon
      if (CoordMath.overlaps(blockReferenceStart, blockReferenceEnd, exonReferenceStart, exonReferenceEnd)) {
        overlap += CoordMath.getOverlap(blockReferenceStart, blockReferenceEnd, exonReferenceStart, exonReferenceEnd)
      }

      // Consume either the next block or next exon, based on which has the earlier end position
      if (blockReferenceEnd < exonReferenceEnd) {
        blocks.next()
      }
      else {
        exons.next()
      }
    }

    overlap
  }
}

/**
  * Metrics produced by `EstimateRnaSeqInsertSize` to describe the distribution of insert sizes within an
  * RNA-seq experiment.  The insert sizes are computed in "transcript space", accounting for spliced
  * alignments, in order to get a true estimate of the size of the DNA fragment, not just it's span on
  * the genome.
  *
  * @param pair_orientation The orientation of the reads within a read-pair relative to each other.
  *                         Possible values are FR, RF and TANDEM.
  * @param read_pairs The number of read pairs observed with the `pair_orientation`.
  * @param mean The mean insert size of the read pairs.
  * @param standard_deviation The standard deviation of the insert size of the read pairs.
  * @param median The median insert size of the read pairs.
  * @param min The minimum observed insert size of the read pairs.
  * @param max The maximum observed insert size of the read pairs.
  * @param median_absolute_deviation The median absolution deviation of the read pairs.
  */
case class InsertSizeMetric(pair_orientation: PairOrientation,
                            read_pairs: Long = 0,
                            mean: Double = 0,
                            standard_deviation: Double = 0,
                            median: Double = 0,
                            min: Long = 0,
                            max: Long = 0,
                            median_absolute_deviation: Double = 0
) extends Metric
