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

package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.ReadStructure.SubRead
import com.fulcrumgenomics.util.miseq._
import com.fulcrumgenomics.util.{SampleBarcode => _, _}
import dagr.commons.CommonsDef.{DirPath, FilePath, PathToBam, PathToFastq}
import dagr.commons.io.PathUtil
import dagr.commons.util.{LazyLogging, Logger}
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.util.{ProgressLogger => _, _}
import htsjdk.samtools.{SAMRecord, _}

object DemuxFromInlineSampleBarcode {
  /** The name of the sample for unmatched reads. */
  val UnmatchedSampleId: String = "unmatched"

  /** Creates the sample output BAM for the given sample. */
  private[fastq] def sampleOutputBam(output: DirPath, sample: Sample): PathToBam = {
    val sampleBarcode: SampleBarcode = sample.sampleBarcode.getOrElse(throw new IllegalArgumentException(s"Sample barcode missing for sample: ${sample.sampleName}"))
    output.resolve(PathUtil.sanitizeFileName(s"${sample.sampleId}-${sample.sampleName}-${sampleBarcode}.bam"))
  }

  /** Creates a simple [[SAMFileHeader]] for the given sample. */
  private[fastq] def createSamFileHeader(sample: Sample): SAMFileHeader = {
    val header: SAMFileHeader = new SAMFileHeader
    val libraryId: String     = sample.libraryId.getOrElse(sample.sampleName)
    val platformUnit: String  = sample.sampleBarcode.map(sb => s"${sample.sampleName}.$sb").getOrElse(sample.sampleName)
    val readGroup             = new SAMReadGroupRecord(sample.sampleId)
    readGroup.setSample(sample.sampleName)
    readGroup.setLibrary(libraryId)
    readGroup.setPlatformUnit(platformUnit)
    readGroup.setPlatform("illumina")
    sample.description.foreach(readGroup.setDescription)
    header.addReadGroup(readGroup)
    header.setSortOrder(SortOrder.queryname)
    header
  }

  /** Creates SAM file writers for each sample in the Sample Sheet, returning them in the same order as specified
    * by each sample's ordinal.  An extra writer is appended for records that do not match a sample barcode. */
  private[fastq] def createSamFileWriters(sampleSheet: SampleSheet, output: DirPath, unmatched: String): List[SAMFileWriter] = {
    val writers: List[SAMFileWriter] = sampleSheet.map { sample =>
      val outputBam: FilePath   = sampleOutputBam(output, sample)
      val header: SAMFileHeader = createSamFileHeader(sample)
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputBam.toFile)
      writer
    }.toList
    val unmatchedWriter = {
      val outputFile: FilePath  = output.resolve(unmatched)
      val header: SAMFileHeader = createSamFileHeader(new Sample(writers.length, UnmatchedSampleId, UnmatchedSampleId, Some(UnmatchedSampleId), None, None))
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile.toFile)
      writer
    }
    writers ++ List(unmatchedWriter)
  }

  /** Gets the quality format of the FASTQs. */
  private def determineQualityFormat(fastqs: Seq[PathToFastq], qualityFormat: Option[FastqQualityFormat] = None, logger: Option[Logger] = None): FastqQualityFormat = {
    val readers = fastqs.map { fastq => new FastqReader(fastq.toFile) }
    val detector: QualityEncodingDetector = new QualityEncodingDetector
    detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers:_*)
    readers.foreach(_.close())
    val format = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, qualityFormat.orNull)
    if (detector.isDeterminationAmbiguous) {
      logger.foreach(_.warning("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities."))
    }
    logger.foreach(_.info(String.format("Auto-detected quality format as: %s.", format)))
    format
  }
}

@clp(
  description =
    """
      |Demultiplexes samples from FASTQ given an inline sample barcode.
      |
      |The sample barcode for each sample in the sample sheet will be compared against the sample barcode bases in
      |the FASTQs, to assign each read to a sample.  Reads that to not match any sample within the given error tolerance
      |will be placed in the 'unmatched' file.
      |
      |The output directory will contain one BAM file per sample in the sample sheet, plus a BAM for reads that could
      |not be assigned to a sample given the criteria.  The output files will be the concatenation of sample id, sample
      |name, and sample barcode bases (expected not observed), delimited by "-".  A metrics file will also be output
      |providing analagous information to the metric desribed here:
      |https://broadinstitute.github.io/picard/picard-metric-definitions.html#SampleBarcodeMetric
      |
      |A single read structure for all bases in read one and read two should be given.  The read structures
      |may contain sample barcode bases ('B'), molecular identifier bases ('M'), template bases ('T'), and bases to
      |skip ('S'). Both reads must have teplate bases.  Any molecular identifiers will be concatenated using the '-'
      |delimiter and placed in the given SAM record tag ("RX" by default).  Similarly, the sample barcode bases from the
      |given read will be placed in the "BC" tag.
      |
      |The read group's sample id, sample name, and library id all correspond to the similarly named values in the
      |sample sheet.  Library id will be the sample id if not found, and the platform unit will be the sample name
      |concatenated with the sample barcode bases delimited by a ".".
      |
      |The sample section of the sample sheet should contain information related to each sample with the following
      |columns
      |  - Sample Identifier:  Sample_ID
      |  - Sample Name:        Sample_Name
      |  - Library Identifier: Library_ID
      |  - Sample Project:     Sample_Project
      |  - Description:        Description
      |  - Read One Inline Sample Barcode Bases: Index or R1_Barcode_Bases*
      |  - Read Two Inline Sample Barcode Bases: Index2 or R2_Barcode_Bases*
      |If the --use-index-columns is true, the sample sheet's 'index' and 'index2' columns will be use for the given
      |sample's barcodes for read one and read two respectively.  Otherwise, the columns with names 'R1_Barcode_Bases'
      |and 'R2_Barcode_Bases' are used.
      |
      |The read structures will be used to extract the observed sample barcode and molecular identifiers from each
      |read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in the sample sheet
      |and associated read structures.  Please see the following link for the read structure:
      |  https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/parser/ReadStructure.html
    """,
  group=ClpGroups.Fastq
)
class DemuxFromInlineSampleBarcode
(@arg(flag="1", doc="One or more input fastq files for read one.")       val inputReadOne:  Seq[PathToFastq],
 @arg(flag="2", doc="One or more input fastq files for read two (omit for fragment reads).", minElements=0) val inputReadTwo: Seq[PathToFastq]=Seq.empty,
 @arg(flag="o", doc="The output directory in which to place sample BAMs.")   val output: DirPath,
 @arg(flag="s", doc="The sample sheet.") val sampleSheet: FilePath,
 @arg(flag="r", doc="The read structure for the first read.") val readStructureReadOne: ReadStructure,
 @arg(flag="R", doc="The read structure for the second read.") val readStructureReadTwo: Option[ReadStructure] = None,
 @arg(flag="m", doc="Per-barcode and per-lane metrics written to this file.")
 val metrics: FilePath,
 @arg(doc="Use the index and index2 columns in the sample sheet to identify R1's and read R2's inline sample barcode for each sample.")
 val useIndexColumns: Boolean = false,
 @arg(flag="u", doc="Output BAM file name for the unmatched records.")
 val unmatched: String = DemuxFromInlineSampleBarcode.UnmatchedSampleId + ".bam",
 @arg(flag="q", doc="A value describing how the quality values are encoded in the FASTQ.  Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.")
 val qualityFormat: Option[FastqQualityFormat] = None,
 @arg(doc="Minimum base quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.") val minQ: Int = 0,
 @arg(doc="Maximum base quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.") val maxQ: Int = SAMUtils.MAX_PHRED_SCORE,
 @arg(doc="Maximum mismatches for a barcode to be considered a match.") val maxMismatches: Int = 1,
 @arg(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.") val minMismatchDelta: Int = 2,
 @arg(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.") val maxNoCalls: Int = 2,
 @arg(flag="t", doc="The SAM tag for any molecular barcode.  If multiple molecular barcodes are specified, they will be concatenated and stored here.")
    val umiTag: String = "RX",
 @arg(flag="n", doc="The number of threads to use while de-multiplexing") val numThreads: Int = 1
) extends FgBioTool with LazyLogging {

  import CustomSampleSheet._
  import DemuxFromInlineSampleBarcode._

  private val readStructures = Seq(Some(readStructureReadOne), readStructureReadTwo).flatten

  /** The number of records to batch before demultiplexing in parallel. */
  private val recordBatchSize = 1e6.toInt

  validate(inputReadTwo.isEmpty || inputReadOne.size == inputReadTwo.size, "Number of read one and read two files should match")

  Io.assertReadable(inputReadOne ++ inputReadTwo)
  Io.assertReadable(sampleSheet)
  Io.assertWritableDirectory(output)

  validate(readStructures.flatMap(_.sampleBarcode).nonEmpty, s"No sample barcodes found in read structures: " + readStructures.map(_.toString).mkString(", "))

  override def execute(): Unit = {
    // Get the FASTQ quality encoding format
    val qualityFormat = this.qualityFormat.getOrElse(determineQualityFormat(inputReadOne, this.qualityFormat, Some(this.logger)))

    // Read in the sample sheet, and override the sample barcode bases if necessary to point to the correct
    // column in the sample sheet.
    val sampleSheet: SampleSheet = maybeCustomSampleSheet(this.sampleSheet, useIndexColumns=useIndexColumns)

    val metricsCollector = new SampleBarcodeMetricCollector(samples=sampleSheet.toSeq, readStructures=readStructures)

    // Create the files for reading and writing
    val r1Readers = inputReadOne.flatMap(FastqSource(_))
    val r2Readers = inputReadTwo.flatMap(FastqSource(_))
    val writers   = createSamFileWriters(sampleSheet, this.output, this.unmatched)
    val converter = new FastqDemultiplexer(
      headers = writers.map(_.getFileHeader),
      samples          = sampleSheet.toSeq,
      readStructures   = readStructures,
      umiTag           = umiTag,
      qualityFormat    = qualityFormat,
      minQ             = minQ,
      maxQ             = maxQ,
      maxMismatches    = maxMismatches,
      minMismatchDelta = minMismatchDelta,
      maxNoCalls       = maxNoCalls
    )

    // Convert the records
    import com.fulcrumgenomics.FgBioDef.ParSupport
    val demuxRecords = if (r2Readers.isEmpty) {
      r1Readers
        .grouped(recordBatchSize)
        .flatMap { records =>
          records.parWith(parallelism=numThreads).map(rec => converter.demultiplex(rec))
        }
    }
    else {
      // Add read numbers if not present
      val reader1 = r1Readers.map { r => r.copy(readNumber=Some(r.readNumber.getOrElse(1))) }
      val reader2 = r2Readers.map { r => r.copy(readNumber=Some(r.readNumber.getOrElse(2))) }
      reader1.zip(reader2)
        .grouped(recordBatchSize)
        .flatMap { records =>
          records.parWith(parallelism=numThreads).map { case (r1, r2) => converter.demultiplex(r1, r2) }
        }
    }

    // Write the records out in its own thread
    val progress = new ProgressLogger(this.logger, unit=1e6.toInt)
    // Developer Note: toStream ensures that writing and metrics collection doesn't happen asynchronously given the
    // parallelism above.
    demuxRecords.toStream.foreach { demuxRecord =>
      val writer = writers(demuxRecord.sampleIndex)
      demuxRecord.records.foreach { rec =>
        writer.addAlignment(rec)
        progress.record(rec)
      }
      metricsCollector.increment(sampleIndex=demuxRecord.sampleIndex, numMismatches=demuxRecord.numMismatches)
    }

    // Close the writer; NB: the inputs close automatically
    writers.foreach(_.close())

    // Write the metrics
    Metric.write(metrics, metricsCollector.barcodeMetrics)
  }
}

/** Utilities for handling sample barcode bases in non-standard columns in the sample sheet. */
private[fastq] object CustomSampleSheet {

  /** Custom columns in the Illumina Experiment Manager Sample Sheet for sample barcodes. */
  val R1SampleBarcodeBases = "R1_Barcode_Bases"
  val R2SampleBarcodeBases = "R2_Barcode_Bases"

  /** Reads in the sample sheet.  If `useIndexColumns` is `false`, set the sample barcodes for each sample to be
    * the inline sample barcode bases stored in the columns with names [[R1SampleBarcodeBases]]
    * and [[R2SampleBarcodeBases]].
    */
  def maybeCustomSampleSheet(sampleSheet: FilePath, useIndexColumns: Boolean): SampleSheet = {
    val originalSampleSheet = SampleSheet(sampleSheet)
    val modifiedSampleSheet = if (useIndexColumns) originalSampleSheet else new SampleSheet(originalSampleSheet.map(new CustomSample(_)).toSeq)
    // Update the sample barcode for fast access
    modifiedSampleSheet.setSampleBarcodes()
  }

  /** A custom sample class to support inline sample barcodes for read one and read two.  This will ignore the i7 and
    * i5 bases. */
  private class CustomSample(sample: Sample, r1SampleBarcodeKey: String = R1SampleBarcodeBases, r2SampleBarcodeKey: String = R2SampleBarcodeBases)
    extends Sample(
      sampleOrdinal      = sample.sampleOrdinal,
      sampleId           = sample.sampleId,
      sampleName         = sample.sampleName,
      libraryId          = sample.libraryId,
      project            = sample.project,
      description        = sample.description,
      lane               = sample.lane,
      i7IndexBases       = None,
      i5IndexBases       = None,
      sampleBarcode      = None,
      extendedAttributes = sample.extendedAttributes
    ) {
    val r1SampleBarcodeBases: Option[String] = extendedAttribute(r1SampleBarcodeKey)
    val r2SampleBarcodeBases: Option[String] = extendedAttribute(r2SampleBarcodeKey)

    /** Returns the sample barcodes in order of sequencing. */
    override def sampleBarcodeBases: Seq[Option[String]] = {
      Seq(r1SampleBarcodeBases, r2SampleBarcodeBases)
    }
  }
}


/** Collects sample barcode metrics per-sample */
private class SampleBarcodeMetricCollector(val samples: Seq[Sample],
                                           val readStructures: Seq[ReadStructure]) {
  import DemuxFromInlineSampleBarcode.UnmatchedSampleId

  require(readStructures.nonEmpty, "No read structures given.")

  private val SampleBarcodeDelimiter = "-"

  /** The barcode string for reads that do not match. */
  private val noMatchBarcode: String = {
    readStructures.flatMap(_.sampleBarcode).map("N" * _.length).mkString(SampleBarcodeDelimiter)
  }
  require(noMatchBarcode.nonEmpty, "No sample barcodes found in read structures: " + readStructures.map(_.toString).mkString(", "))

  require(samples.flatMap(_.sampleBarcode).map(_.concatenatedBarcode).sorted.distinct.length == samples.length, "Unique sample barcodes required for all samples")

  /** The sample barcode metrics in the same order as samples in the sample sheet, with an additional sample for unmatched reads. */
  private val metrics: Seq[SampleBarcodeMetric] = {
    val sampleMetrics = samples.map { sample =>
      val barcode: String = sample.sampleBarcode match {
        case Some(b) => b.concatenatedBarcode
        case None    => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a sample barcode")
      }
      val libraryId: String = sample.libraryId match {
        case Some(id) => id
        case None     => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a library id")
      }
      SampleBarcodeMetric(barcodeName=sample.sampleName, libraryName=libraryId, barcode=barcode)
    }
    sampleMetrics.toSeq ++ Seq(SampleBarcodeMetric(UnmatchedSampleId, UnmatchedSampleId, noMatchBarcode))
  }

  /** Increments the metrics for the given sample ordinal, which matched with the given number of mismatches.  The
    * sample ordinal less than zero indicates the "unmatched" sample. */
  def increment(sampleIndex: Int, numMismatches: Int): Unit = {
    if (0 <= sampleIndex && sampleIndex < this.metrics.length-1) {
      val metric = metrics(sampleIndex)
      metric.reads += 1
      metric.pf_reads += 1
      if (numMismatches == 0) {
        metric.perfect_matches += 1
        metric.pf_perfect_matches += 1
      }
      else if (numMismatches == 1) {
        metric.one_mismatch_matches += 1
        metric.pf_one_mismatch_matches += 1
      }
    }
    else {
      require(sampleIndex < this.metrics.length, s"Sample ordinal out of range: $sampleIndex")
      val metric = this.metrics.last
      metric.reads += 1
      metric.pf_reads += 1
    }
  }

  /** Finalize the barcode metrics and return them. */
  def barcodeMetrics: Iterable[SampleBarcodeMetric] = {
    SampleBarcodeMetric.finalizeMetrics(this.metrics.map { m => (m.barcode, m) }.toMap, this.noMatchBarcode)
    metrics
  }
}



private[fastq] object FastqDemultiplexer {
  /** A class to store the sample ordinal and associated demultiplexed [[SAMRecord]]s.
    * @param sampleIndex the sample index to which the records were assigned (0-based)
    * @param numMismatches the # of mismatches if it matched a sample with a sample barcode, -1 otherwise (the 'unmatched'
    *                      sample).
    * @param records the records, one for each read that have template bases.
    */
  case class DemuxRecord(sampleIndex: Int, numMismatches: Int, records: Seq[SAMRecord])

  /** Counts the nucleotide mismatches between two strings of the same length.  Ignores no calls in expectedBases.
    * Observed base qualities less than the minimum base quality are counted as mismatches if not a no call.
    */
  private[fastq] def countMismatches(observedBases: Array[Byte], expectedBases: Array[Byte]): Int = {
    var idx = 0
    var count = 0
    while (idx < observedBases.length) {
      val expectedBase = expectedBases(idx)
      val observedBase = observedBases(idx)
      if (!SequenceUtil.isNoCall(expectedBase) && !SequenceUtil.basesEqual(observedBase, expectedBase)) {
        count += 1
      }
      idx += 1
    }
    count
  }
}

/** Performs sample demultiplexes from reads from the same template/fragment.
  *
  * A [[SAMFileHeader]] should be given per sample and a [[ReadStructure]] per read from the same template/fragment.
  * Use the [[demultiplex()]] method to create a [[SAMRecord]] for each read with template bases.  Any molecular barcodes
  * will be extracted and stored in the tag specified by [[umiTag]].
  *
  * @param headers a sequence of [[SAMFileHeader]]s for each sample in the same order as samples in the sample sheet.
  * @param samples the samples ordered by sample ordinal.
  * @param readStructures the read structures, one for each read that will be given to [[demultiplex()]].
  * @param umiTag the tag to store any molecular barcodes.  The barcodes from reads will be delimited by "-".
  * @param qualityFormat the quality format for the FASTQ records.
  * @param minQ the minimum base quality to enforce.
  * @param maxQ the maximum base quality to enforce.
  * @param maxMismatches the maximum mismatches to match a sample barcode.
  * @param minMismatchDelta the minimum difference between number of mismatches in the best and second best barcodes for
  *                         a barcode to be considered a match.
  * @param maxNoCalls the maximum number of no calls in the sample barcode bases allowed for matching.
  */
private class FastqDemultiplexer(val headers: Seq[SAMFileHeader],
                                 val samples: Seq[Sample],
                                 val readStructures: Seq[ReadStructure],
                                 val umiTag: String = "RX",
                                 val qualityFormat: FastqQualityFormat = FastqQualityFormat.Standard,
                                 val minQ: Int = 0,
                                 val maxQ: Int =  SAMUtils.MAX_PHRED_SCORE,
                                 val maxMismatches: Int = 2,
                                 val minMismatchDelta: Int = 1,
                                 val maxNoCalls: Int = 2) {
  import FastqDemultiplexer._

  require(readStructures.nonEmpty, "No read structures were given")
  require(samples.map(_.sampleOrdinal).toList == samples.map(_.sampleOrdinal).sorted, "Samples not provided in ascending order by sample ordinal")
  require(samples.flatMap(_.sampleBarcode).map(_.concatenatedBarcode).sorted.distinct.length == samples.length, "Unique sample barcodes required for all samples")
  require(headers.length == samples.length + 1, s"Not enough SAMFileHeaders; found ${headers.length} but should have ${samples.length+1}.")
  require(headers.forall(_.getReadGroups.size == 1), "All headers must have one read group. Ordinals are: " + headers.zipWithIndex.filter(_._1.getReadGroups.size != 1).map(_._2).mkString(", "))

  /** The concatenation of all read structures read structure. */
  private val fullReadStructure: ReadStructure = new ReadStructure(segs = readStructures.flatten, resetOffsets = true)

  /** The quality format converter. */
  private val solexaQualityConverter: SolexaQualityConverter = SolexaQualityConverter.getSingleton

  /** The sub-reads for the sample barcodes across all read structures. */
  private val sampleBarcodeSegments = fullReadStructure.sampleBarcode

  /** Get the sample barcode bases across all bases. */
  private def sampleBarcodeBases(allBases: String): String = {
    ReadStructure.structureRead(bases=allBases, segments=sampleBarcodeSegments).map(_.bases).mkString
  }

    /** Gets the sample ordinal given the bases from all cycles and the number of mismatches between the bases and matched
      * sample barcode.  If no match is found, the ordinal for the oridinal for the 'unmatched" sample is given. */
  private def matchSampleBarcode(bases: String): (Int, Int) = {
    val observedBarcode = sampleBarcodeBases(bases).getBytes
    val numNoCalls      = observedBarcode.count(base => SequenceUtil.isNoCall(base))

    // Get the best and second best sample barcode matches.
    val (bestSampleIndex, bestMismatches, secondBestMismatches) = if (numNoCalls <= maxNoCalls) {
      samples.map { sample =>
        val expectedBarcode = sample.sampleBarcode
          .getOrElse(unreachable(s"Sample with id '${sample.sampleId}' did not have a sample barcode"))
          .barcodeBytes
        val numMismatches = countMismatches(observedBarcode, expectedBarcode)
        (sample.sampleOrdinal-1, numMismatches)
      }
      .toList.sortBy(_._2).take(2) match {
        case Nil                              => (-1,           Int.MaxValue, Int.MaxValue)
        case List(bestTuple)                  => (bestTuple._1, bestTuple._2, Int.MaxValue)
        case List(bestTuple, secondBestTuple) => (bestTuple._1, bestTuple._2, secondBestTuple._2)
      }
    }
    else {
      (-1, Int.MaxValue, Int.MaxValue)
    }

    // Make sure we are within the parameter limits and update barcode metrics if necessary
    if (maxMismatches < bestMismatches || maxNoCalls < numNoCalls || (secondBestMismatches - bestMismatches) < minMismatchDelta) {
      (this.samples.size, -1)
    }
    else {
      (bestSampleIndex, bestMismatches)
    }
  }

  /** Demultiplexes a given set of reads from the same template.  The same number of reads should be given as read
    * structures.
    *
    * The sample barcoded bases from each read are extracted and concatenated in the same order as the given reads. They
    * are matched against the sample barcode bases for each sample.  The
    * */
  def demultiplex(reads: FastqRecord*): DemuxRecord = {
    require(reads.nonEmpty, "No reads given for demultiplexing.")
    require(reads.length == readStructures.length, s"Expected the same number of reads ('${reads.length}') as read structures ('${readStructures.length}').")
    val numReadsWithTemplate = readStructures.count(_.template.nonEmpty)
    val pairedEnd            = numReadsWithTemplate == 2
    val allBases             = reads.map(_.bases).mkString

    // Get the sample ordinal
    val (sampleIndex, numMismatches) = matchSampleBarcode(allBases)

    // Get the molecular barcode bases
    val molecularBarcode = ReadStructure.structureRead(bases=allBases, segments=this.fullReadStructure.molecularBarcode).map(_.bases).filter(_.nonEmpty).mkString("-")

    // Create the SAMRecords.
    val header   = headers(sampleIndex)
    val sampleId = header.getReadGroups.get(0).getId
    val records  = reads.zip(readStructures).flatMap { case (read, readStructure) =>
      val templateSegments = readStructure.template
      if (templateSegments.isEmpty) {
        None
      }
      else {
        val subReads: Seq[SubRead] = ReadStructure.structureReadWithQualities(bases=read.bases, qualities=read.quals, segments=templateSegments)
        val bases                  = subReads.map(_.bases).mkString
        val quals                  = subReads.flatMap(_.quals).mkString

        // Build the [[SAMRecord]]
        val record  = new SAMRecord(header)
        record.setReadName(read.name)
        record.setReadString(bases)
        record.setReadUnmappedFlag(true)
        if (pairedEnd) {
          require(read.readNumber.isDefined, s"Read number not defined for :${read.name}.")
          record.setReadPairedFlag(true)
          record.setMateUnmappedFlag(true)
          if (read.readNumber.forall(_ == 1)) {
            record.setFirstOfPairFlag(true)
          }
          else {
            require(read.readNumber.forall(_ == 2), s"Expected read number 2 but found '${read.readNumber.getOrElse("None")}' for read: ${read.name}")
            record.setSecondOfPairFlag(true)
          }
        }
        record.setAttribute(ReservedTagConstants.READ_GROUP_ID, sampleId)
        if (molecularBarcode.nonEmpty) record.setAttribute("RX", molecularBarcode)

        // Have fun with setting the qualities
        val readQualsBytes = quals.getBytes
        convertQuality(readQualsBytes, qualityFormat)
        readQualsBytes.foreach { qual =>
          val uQual: Int = qual & 0xff
          require(minQ <= uQual && uQual <= maxQ, s"Base quality $uQual is not in the range $minQ ... $maxQ for read ${read.name}")
        }
        record.setBaseQualities(readQualsBytes)

        Some(record)
      }
    }

    DemuxRecord(sampleIndex=sampleIndex, numMismatches=numMismatches, records=records)
  }

  /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
  private[fastq] def convertQuality(quals: Array[Byte], version: FastqQualityFormat) = {
    version match {
      case FastqQualityFormat.Standard => SAMUtils.fastqToPhred(quals)
      case FastqQualityFormat.Solexa => solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals)
      case FastqQualityFormat.Illumina => solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals)
    }
  }
}
