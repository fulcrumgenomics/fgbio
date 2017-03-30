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

import java.io.{Closeable, Serializable}

import com.fulcrumgenomics.FgBioDef.unreachable
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.fastq.FastqDemultiplexer.{DemuxRecord, DemuxResult}
import com.fulcrumgenomics.illumina.{Sample, SampleSheet}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.ReadStructure.SubRead
import com.fulcrumgenomics.util.{ReadStructure, SampleBarcodeMetric, SampleBarcode => _, _}
import dagr.commons.CommonsDef.{DirPath, FilePath, PathPrefix, PathToBam, PathToFastq, yieldAndThen}
import dagr.commons.io.PathUtil
import dagr.commons.util.{LazyLogging, Logger}
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.util.{ProgressLogger => _, _}
import htsjdk.samtools.{SAMRecord, _}

import scala.concurrent.forkjoin.ForkJoinPool

object DemuxFastqs {

  /** The name of the metrics file if none is given by the user. */
  val DefaultDemuxMetricsFileName = "demux_barcode_metrics.txt"

  /** The name of the sample for unmatched reads. */
  val UnmatchedSampleId: String = "unmatched"

  /** The maximum # of records in RAM per writer. */
  private[fastq] val MaxRecordsInRamPerSamFileWriter: Int = 1e6.toInt

  /** Decides whether or not to use asynchronous IO for the writers. This has a big performance benefit at the cost of
    * some RAM. RAM may balloon if there is a need to sort the output. */
  private[fastq] val UseAsyncIo: Boolean  = true

  /** The number of records to batch when demultiplexing in parallel. */
  private val DemuxBatchRecordsSize = 1e5.toInt

  /** Creates the sample output BAM path for the given sample. */
  private[fastq] def outputPrefixFrom(output: DirPath, sample: Sample): PathPrefix = {
    val sampleBarcode = sample.sampleBarcodeString
    require(sampleBarcode.nonEmpty, s"Sample barcode missing for sample: ${sample.sampleName}")
    output.resolve(PathUtil.sanitizeFileName(s"${sample.sampleId}-${sample.sampleName}-$sampleBarcode"))
  }

  /** Gets the quality format of the FASTQs. */
  private def determineQualityFormat(fastqs: Seq[PathToFastq], expectedQualityFormat: Option[FastqQualityFormat] = None, logger: Option[Logger] = None): FastqQualityFormat = {
    val readers = fastqs.map { fastq => new FastqReader(fastq.toFile) }
    val detector: QualityEncodingDetector = new QualityEncodingDetector
    detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers:_*)
    readers.foreach(_.close())
    val format = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, expectedQualityFormat.orNull)
    if (detector.isDeterminationAmbiguous) {
      logger.foreach(_.warning("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities."))
    }
    logger.foreach(_.info(String.format("Auto-detected quality format as: %s.", format)))
    format
  }

  /** Create a sample with sample barcodes extracted from a custom column. */
  private[fastq] def withCustomSampleBarcode(sample: Sample, columnForSampleBarcode: String): Sample = {
    new Sample(
      sampleOrdinal      = sample.sampleOrdinal,
      sampleId           = sample.sampleId,
      sampleName         = sample.sampleName,
      libraryId          = sample.libraryId,
      project            = sample.project,
      description        = sample.description,
      lane               = sample.lane,
      i7IndexBases       = sample.i7IndexBases,
      i5IndexBases       = sample.i5IndexBases,
      extendedAttributes = sample.extendedAttributes
    ) {
      override val sampleBarcodes: Seq[Option[String]] = {
        val barcode = extendedAttribute(columnForSampleBarcode)
        require(barcode.nonEmpty, s"Sample barcode not found in column '$columnForSampleBarcode' for sample id '${sample.sampleId}'.")
        Seq(barcode)
      }
    }
  }

  /** Create the unmatched sample with the given sample ordinal. */
  private[fastq] def unmatchedSample(sampleOrdinal: Int, readStructures: Seq[ReadStructure]): Sample = {
    val noMatchBarcode: String = readStructures.flatMap(_.sampleBarcode).map("N" * _.length).mkString("-")
    require(noMatchBarcode.nonEmpty, "No sample barcodes found in read structures: " + readStructures.map(_.toString).mkString(", "))
    Sample(sampleOrdinal=sampleOrdinal, sampleId=UnmatchedSampleId, sampleName=UnmatchedSampleId, libraryId=UnmatchedSampleId, i7IndexBases=Some(noMatchBarcode))
  }

  /** A little class to group reads from the same template together. */
  private class ZippedIterator(sources: Seq[Iterator[FastqRecord]]) extends Iterator[Seq[FastqRecord]] {
    require(sources.nonEmpty, "No sources provided")
    def hasNext(): Boolean = sources.exists(_.hasNext)
    def next(): Seq[FastqRecord] = {
      if (!this.hasNext) throw new NoSuchElementException("Calling next() when hasNext() is false.")
      require(sources.forall(_.hasNext) == sources.head.hasNext, "Sources are out of sync.")
      val records = sources.map(_.next)
      // Check that the FASTQ records all have the same name
      require(records.forall(_.name == records.head.name), "Sources out of sync, found read names: \n" + records.map(_.name).mkString("\n"))
      records
    }
  }

  private[fastq] def toSampleOutputPrefix(sample: Sample, isUnmatched: Boolean, output: DirPath, unmatched: String): PathPrefix = {
    if (isUnmatched) PathUtil.removeExtension(output.resolve(unmatched)) else outputPrefixFrom(output, sample)
  }

  /** Creates a demultiplexing iterator that performs demultiplexing in parallel.
    *
    * @param sources the FASTQ sources, one per read.
    * @param demultiplexer the demultiplexer to use.  The demultiplexer's [[FastqDemultiplexer.demultiplex()]] method
    *                      expects the same number of reads as sources.
    */
  def demultiplexingIterator(sources: Seq[Iterator[FastqRecord]],
                             demultiplexer: FastqDemultiplexer,
                             threads: Int,
                             batchSize: Int = DemuxBatchRecordsSize): Iterator[DemuxResult] = {

    require(demultiplexer.expectedNumberOfReads == sources.length,
      s"The demultiplexer expects ${demultiplexer.expectedNumberOfReads} reads but ${sources.length} FASTQ sources given.")

    val zippedIterator = new ZippedIterator(sources=sources)
    if (threads > 1) {
      // Developer Note: Iterator does not support parallel operations, so we need to group together records into a
      // [[List]] or [[Seq]].  A fixed number of records are grouped to reduce memory overhead.
      import com.fulcrumgenomics.FgBioDef.ParSupport
      val pool = new ForkJoinPool(threads)
      zippedIterator
        .grouped(batchSize)
        .flatMap { batch =>
          batch
            .parWith(pool=pool)
            .map { readRecords => demultiplexer.demultiplex(readRecords: _*) }
        }.seq // Developer Note: toStream ensures that the only parallelism is within the flatMap
    }
    else {
      zippedIterator.map { readRecords => demultiplexer.demultiplex(readRecords: _*) }
    }
  }
}

@clp(
  description =
    """
      |Performs sample demultiplexing on FASTQs.
      |
      |The sample barcode for each sample in the sample sheet will be compared against the sample barcode bases extracted from
      |the FASTQs, to assign each read to a sample.  Reads that do not match any sample within the given error tolerance
      |will be placed in the 'unmatched' file.
      |
      |The output directory will contain one BAM file per sample in the sample sheet or metadata CSV file, plus a BAM for
      |reads that could not be assigned to a sample given the criteria.  The output file names will be the concatenation
      |of sample id, sample name, and sample barcode bases (expected not observed), delimited by "-".  A metrics file
      |will also be output providing analogous information to the metric described here:
      |https://broadinstitute.github.io/picard/picard-metric-definitions.html#SampleBarcodeMetric
      |
      |Alternatively, gzipped FASTQs can be written using the "--output-fastqs=true" option instead of BAMs.  For paired
      |end data, the output will have the suffix "_R1.fastq.gz" and "_R2.fastq.gz" for read one and read two respectively.
      |The sample barcode and molecular barcodes (concatenated) will be appended to the read name and delimited by a
      |colon.
      |
      |The output base qualities will be standardized to Sanger/SAM format.
      |
      |FASTQs and associated read structures for each read should be given:
      |- a single fragment read should have one FASTQ and one read structure
      |- paired end reads should have two FASTQs and two read structures
      |- a dual-index sample with paired end reads should have four FASTQs and four read structures given: two for the
      |  two index reads, and two for the template reads.
      |
      |The read structures may contain sample barcode bases ('B'), molecular identifier bases ('M'), template bases ('T'),
      |and bases to skip ('S'). Both reads must have template bases.  Any molecular identifiers will be concatenated using
      |the '-' delimiter and placed in the given SAM record tag ("RX" by default).  Similarly, the sample barcode bases
      |from the given read will be placed in the "BC" tag.
      |
      |Metadata about the samples should be given in either an Illumina Experiment Manager sample sheet or a metadata CSV
      |file.  Formats are described in detail below.
      |
      |The read structures will be used to extract the observed sample barcode, template bases, and molecular identifiers
      |from each read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in
      |the sample metadata and associated read structures.
      |
      |## Sample Sheet
      |The read group's sample id, sample name, and library id all correspond to the similarly named values in the
      |sample sheet.  Library id will be the sample id if not found, and the platform unit will be the sample name
      |concatenated with the sample barcode bases delimited by a ".".
      |
      |The sample section of the sample sheet should contain information related to each sample with the following
      |columns:
      |  - Sample_ID:   The sample identifier unique to the sample in the sample sheet.
      |  - Sample_Name: The sample name.
      |  - Library_ID:  The library Identifier.  The combination sample name and library identifier should be unique
      |                 across the samples in the sample sheet.
      |  - Description: The description of the sample, which will be placed in the description field in the output BAM's
      |                 read group.  This column may be omitted.
      |
      |Additionally, the sample barcode should be specified in a column named 'Sample_Barcode'.  The name of the column
      |containing the sample barcode can be changed using the --column-for-sample-barcode option.  If the sample barcode
      |is present across multiple reads (ex. dual-index, or inline in both reads of a pair), then the expected barcode
      |bases from each read should be concatenated and placed in the 'Sample_Barcode' column.  The concatenation should
      |be in the same order as the order of the reads' FASTQs and read structures given to this tool.
      |
      |## Metadata CSV
      |In lieu of a sample sheet, a simple CSV file may be provided with the necessary metadata.  This file should
      |contain the same columns as described above for the sample sheet (Sample_ID, Sample_Name, Library_ID, and
      |Description).
      |
      |## Example Command Line
      |
      |As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both reading a sample
      |barcode, as well as an in-line 8bp sample barcode in read one, the command line would be
      |  --inputs r1.fq i1.fq i2.fq r2.fq --read-structures 8B92T 8B 8B 100T \
      |    --sample-sheet SampleSheet.csv --metrics metrics.txt --output output_folder
    """,
  group=ClpGroups.Fastq
)
class DemuxFastqs
(@arg(flag="i", doc="One or more input fastq files each corresponding to a sub-read (ex. index read, read one, read two, fragment).") val inputs: Seq[PathToFastq],
 @arg(flag="o", doc="The output directory in which to place sample BAMs.") val output: DirPath,
 @arg(flag="x", doc="A file containing the metadata about the samples.") val metadata: FilePath,
 @arg(flag="r", doc="The read structure for each of the FASTQs.") val readStructures: Seq[ReadStructure],
 @arg(flag="m", doc="The file to which per-barcode metrics are written.  If none given, a file named 'demux_barcode_metrics.txt' will be written to the output directory.") val metrics: Option[FilePath] = None,
 @arg(flag="c", doc="The column name in the sample sheet or metadata CSV for the sample barcode.") val columnForSampleBarcode: String = "Sample_Barcode",
 @arg(flag="u", doc="Output BAM file name for the unmatched records.") val unmatched: String = DemuxFastqs.UnmatchedSampleId + ".bam",
 @arg(flag="q",
   doc="""A value describing how the quality values are encoded in the FASTQ.
      |Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66),
      |Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard
      |for phred scaled scores with a character shift of 33.  If this value
      |is not specified, the quality format will be detected automatically.
  """)
val qualityFormat: Option[FastqQualityFormat] = None,
 @arg(flag="t", doc="The number of threads to use while de-multiplexing. The performance does not increase linearly with the # of threads and seems not to improve beyond 2-4 threads.") val threads: Int = 1,
 @arg(doc="Maximum mismatches for a barcode to be considered a match.") val maxMismatches: Int = 1,
 @arg(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.") val minMismatchDelta: Int = 2,
 @arg(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.") val maxNoCalls: Int = 2,
 @arg(doc="The sort order for the output sam/bam file (typically unsorted or queryname).") val sortOrder: SortOrder = SortOrder.unsorted,
 @arg(doc="The SAM tag for any molecular barcode.  If multiple molecular barcodes are specified, they will be concatenated and stored here.") val umiTag: String = ConsensusTags.UmiBases,
 @arg(doc="The platform unit (typically '<flowcell-barcode>-<samle-barcode>.<lane>')") val platformUnit: Option[String] = None,
 @arg(doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
 @arg(doc="Predicted median insert size, to insert into the read group header") val predictedInsertSize: Option[Integer] = None,
 @arg(doc="Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)") val platformModel: Option[String] = None,
 @arg(doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
 @arg(doc="Date the run was produced, to insert into the read group header") val runDate: Option[Iso8601Date] = None,
 @arg(doc="Output FASTQs (.fastq.gz) in addition to BAM files") val outputFastqs: Boolean = false
) extends FgBioTool with LazyLogging {

  import DemuxFastqs._

  private[fastq] val metricsPath = metrics.getOrElse(output.resolve(DefaultDemuxMetricsFileName))

  validate(inputs.length == readStructures.length, "The same number of read structures should be given as FASTQs.")
  validate(readStructures.flatMap(_.sampleBarcode).nonEmpty, s"No sample barcodes found in read structures: " + readStructures.map(_.toString).mkString(", "))

  private val pairedEnd = readStructures.count(_.template.nonEmpty) match {
    case 1 => false
    case 2 => true
    case n => invalid(s"Found $n read structures with template bases but expected 1 or 2.")
  }

  Io.assertReadable(inputs)
  Io.assertReadable(metadata)
  Io.assertWritableDirectory(output)

  /** The quality format converter. */
  private val solexaQualityConverter: SolexaQualityConverter = SolexaQualityConverter.getSingleton

  /** Create a sample sheet from either the input sample sheet path or the metadata CSV. */
  private val sampleSheet: SampleSheet = {
    val lines = Io.readLines(metadata).toSeq
    if (lines.exists(_.contains("[Data]"))) {
      logger.info("Assuming input metadata file is an Illumina Experiment Manager Sample Sheet.")
      SampleSheet(lines.toIterator, lane=None)
    }
    else {
      logger.info("Assuming input metadata file is simple CSV file.")
      SampleSheet(Iterator("[Data]") ++ lines, lane=None)
    }
  }

  override def execute(): Unit = {
    // Get the FASTQ quality encoding format
    val qualityFormat = this.qualityFormat.getOrElse(determineQualityFormat(inputs, this.qualityFormat, Some(this.logger)))

    // Read in the sample sheet and create the sample information
    val samplesFromSampleSheet = sampleSheet.map(s => withCustomSampleBarcode(s, columnForSampleBarcode))
    val samples                = samplesFromSampleSheet.toSeq :+ unmatchedSample(samplesFromSampleSheet.size, this.readStructures)
    val sampleInfos            = samples.map(sample => SampleInfo(sample, sample.sampleName == UnmatchedSampleId))
    val sampleToWriter         = sampleInfos.map { info => info.sample -> toWriter(info) }.toMap

    // Validate that the # of sample barcode bases in the read structure matches the # of sample barcode in the sample sheet.
    {
      val rsNumSampleBarcodeBases = readStructures.map(_.sampleBarcode.map(_.length).sum).sum
      samples.foreach { sample =>
        val numSampleBarcodeBases = sample.sampleBarcodeBytes.length
        require(numSampleBarcodeBases == rsNumSampleBarcodeBases,
          s"The number of sample barcodes bases did not match; read structures: $rsNumSampleBarcodeBases sample (${sample.sampleId}): $numSampleBarcodeBases"
        )
      }
    }

    // Create the files for reading and writing
    val demultiplexer = new FastqDemultiplexer(
      sampleInfos      = sampleInfos,
      readStructures   = this.readStructures,
      umiTag           = umiTag,
      maxMismatches    = maxMismatches,
      minMismatchDelta = minMismatchDelta,
      maxNoCalls       = maxNoCalls
    )

    val progress = new ProgressLogger(this.logger, unit=1e6.toInt)

    // An iterator that uses the given fastq demultiplexer to convert FASTQ records from the same fragment/template to
    // DemuxRecord in parallel
    val sources   = inputs.map(FastqSource(_))
    val iterator  = demultiplexingIterator(
      sources       = sources,
      demultiplexer = demultiplexer,
      threads       = threads
    )

    // Write the records out in its own thread
    iterator.foreach { demuxRecord =>
      demuxRecord.sampleInfo.metric.increment(numMismatches=demuxRecord.numMismatches)
      val writer = sampleToWriter(demuxRecord.sampleInfo.sample)
      demuxRecord.records.foreach { rec =>

        val readQualsBytes = rec.quals.getBytes
        convertQuality(readQualsBytes, qualityFormat, rec.name)

        writer.add(rec.copy(quals=readQualsBytes.mkString))
        progress.record()
      }
    }

    // Close the writer; NB: the inputs close automatically
    sampleToWriter.values.foreach(_.close())

    // Write the metrics
    val metricsMap = sampleInfos.map { sampleInfo => (sampleInfo.metric.barcode, sampleInfo.metric) }.toMap
    val unmatchedBarcode = sampleInfos.find { sampleInfo => sampleInfo.isUnmatched }.getOrElse(unreachable("No unmatched sample."))
    SampleBarcodeMetric.finalizeMetrics(metricsMap, unmatchedBarcode.metric.barcode)
    Metric.write(metricsPath, sampleInfos.map(_.metric))
  }

  private def toWriter(sampleInfo: SampleInfo): DemuxWriter = {
    val sample = sampleInfo.sample
    val isUnmatched = sample.sampleName == UnmatchedSampleId
    val prefix = toSampleOutputPrefix(sample, isUnmatched, output, this.unmatched)

    if (this.outputFastqs) {
      new FastqRecordWriter(prefix, pairedEnd)
    }
    else {
      val readGroup = new SAMReadGroupRecord(sample.sampleId)
      readGroup.setSample(sample.sampleName)
      readGroup.setLibrary(sample.libraryId)
      readGroup.setPlatform("Illumina")
      sample.description.foreach(readGroup.setDescription)
      platformUnit.foreach(readGroup.setPlatformUnit)
      sequencingCenter.foreach(readGroup.setSequencingCenter)
      predictedInsertSize.foreach(readGroup.setPredictedMedianInsertSize)
      runDate.foreach(readGroup.setRunDate)
      platformModel.foreach(readGroup.setPlatformModel)

      val header: SAMFileHeader = new SAMFileHeader
      header.addReadGroup(readGroup)
      header.setSortOrder(sortOrder)
      comments.foreach(header.addComment)

      new SamRecordWriter(PathUtil.pathTo(prefix + ".bam"), header, this.umiTag)
    }
  }

  /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
  private[fastq] def convertQuality(quals: Array[Byte], version: FastqQualityFormat, name: String) = {
    version match {
      case FastqQualityFormat.Standard => SAMUtils.fastqToPhred(quals)
      case FastqQualityFormat.Solexa => solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals)
      case FastqQualityFormat.Illumina => solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals)
    }
    // Check that the converted qualities are same.
    quals.foreach { qual =>
      val uQual: Int = qual & 0xff
      require(0 <= uQual && uQual <= SAMUtils.MAX_PHRED_SCORE, s"Base quality $uQual is not in the range 0 ... ${SAMUtils.MAX_PHRED_SCORE} for read: $name")
    }
  }
}

/** A writer than writes [[DemuxRecord]]s */
private trait DemuxWriter extends Closeable {
  def add(rec: DemuxRecord): Unit
}

/** A writer that writes [[DemuxRecord]]s as [[SAMRecord]]s. */
private class SamRecordWriter(output: PathToBam,
                              val header: SAMFileHeader,
                              val umiTag: String) extends DemuxWriter {
  private val writer = new SAMFileWriterFactory()
    .setUseAsyncIo(DemuxFastqs.UseAsyncIo)
    .setMaxRecordsInRam(DemuxFastqs.MaxRecordsInRamPerSamFileWriter)
    .makeSAMOrBAMWriter(header, header.getSortOrder == SortOrder.unsorted, output.toFile)

  private val rgId: String = this.header.getReadGroups.get(0).getId

  def add(rec: DemuxRecord): Unit = {
    val record  = new SAMRecord(header)
    record.setReadName(rec.name)
    record.setReadString(rec.bases)
    record.setBaseQualities(rec.quals.getBytes)
    record.setReadUnmappedFlag(true)
    if (rec.pairedEnd) {
      record.setReadPairedFlag(true)
      record.setMateUnmappedFlag(true)
      if (rec.readNumber == 1) record.setFirstOfPairFlag(true) else record.setSecondOfPairFlag(true)
    }
    rec.sampleBarcode.foreach(bc => record.setAttribute("BC", bc))
    record.setAttribute(ReservedTagConstants.READ_GROUP_ID, rgId)
    rec.molecularBarcode.foreach(mb => record.setAttribute(umiTag, mb))
    writer.addAlignment(record)
  }

  override def close(): Unit = writer.close()
}

/** A writer that writes [[DemuxRecord]]s as [[FastqRecord]]s. */
private class FastqRecordWriter(prefix: PathPrefix, val pairedEnd: Boolean) extends DemuxWriter {
  private val writers: IndexedSeq[FastqWriter] = {
    (if (pairedEnd) Seq("_R1.fastq.gz", "_R2.fastq.gz") else Seq(".fastq.gz"))
    .map(ext => FastqWriter(Io.toWriter(PathUtil.pathTo(prefix + ext))))
      .toIndexedSeq
  }

  def add(rec: DemuxRecord): Unit = {
    val record = FastqRecord(
      name       = Seq(Some(rec.name), rec.sampleBarcode, rec.molecularBarcode).flatten.mkString(":"),
      bases      = rec.bases,
      quals      = rec.quals,
      comment    = None,
      readNumber = Some(rec.readNumber)
    )
    val writer = record.readNumber match {
      case Some(1) | None => writers.head
      case Some(2)        => writers.last
      case _              => throw new IllegalStateException(s"Read number was invalid: ${record.readNumber}")
    }
    writer.write(record)
  }

  override def close(): Unit = this.writers.foreach(_.close())
}


/** A class to store information about a sample. */
private[fastq] case class SampleInfo(sample: Sample, isUnmatched: Boolean = false) {
  val metric: SampleBarcodeMetric = {
    val barcode: String = this.sample.sampleBarcodeString
    require(barcode.nonEmpty, s"Sample with id '${sample.sampleId}' did not have a sample barcode")
    SampleBarcodeMetric(barcodeName=sample.sampleName, libraryName=sample.libraryId, barcode=barcode)
  }
}

private[fastq] object FastqDemultiplexer {

  /** Stores the minimal information for a single template read. */
  case class DemuxRecord(name: String, bases: String, quals: String, molecularBarcode: Option[String], sampleBarcode: Option[String], readNumber: Int, pairedEnd: Boolean)

  /** A class to store the [[SampleInfo]] and associated demultiplexed [[DemuxRecord]]s.
    * @param sampleInfo the [[SampleInfo]] for the matched sample.
    * @param numMismatches the # of mismatches if it matched a sample with a sample barcode.  This will be [[Int.MaxValue]]
    *                      for the unmatched sample.
    * @param records the records, one for each read that has template bases.
    */
  case class DemuxResult(sampleInfo: SampleInfo, numMismatches: Int, records: Seq[DemuxRecord])

  /** Counts the nucleotide mismatches between two strings of the same length.  Ignores no calls in expectedBases. */
  private[fastq] def countMismatches(observedBases: Array[Byte], expectedBases: Array[Byte]): Int = {
    require(observedBases.length == expectedBases.length, s"observedBases: ${observedBases.length} expectedBases: ${expectedBases.length}")
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

/** Assigns reads from the same fragment/template to a sample.
  *
  * A [[SampleInfo]] should be given per sample and a [[ReadStructure]] per read from the same template/fragment.
  * Use the [[demultiplex()]] method to create a [[DemuxRecord]] for each read with template bases.  Any molecular barcodes
  * will be extracted and stored in the tag specified by [[umiTag]].
  *
  * @param sampleInfos the sample information, one per sample.
  * @param readStructures the read structures, one for each read that will be given to [[demultiplex()]].
  * @param umiTag the tag to store any molecular barcodes.  The barcodes from reads will be delimited by "-".
  * @param maxMismatches the maximum mismatches to match a sample barcode.
  * @param minMismatchDelta the minimum difference between number of mismatches in the best and second best barcodes for
  *                         a barcode to be considered a match.
  * @param maxNoCalls the maximum number of no calls in the sample barcode bases allowed for matching.
  */
private class FastqDemultiplexer(val sampleInfos: Seq[SampleInfo],
                                 val readStructures: Seq[ReadStructure],
                                 val umiTag: String = ConsensusTags.UmiBases,
                                 val maxMismatches: Int = 2,
                                 val minMismatchDelta: Int = 1,
                                 val maxNoCalls: Int = 2) {
  import FastqDemultiplexer._

  require(readStructures.nonEmpty, "No read structures were given")

  {
    val samples = sampleInfos.map(_.sample)
    require(samples.map(_.sampleBarcodeString).sorted.distinct.length == samples.length, "Unique sample barcodes required for all samples")
  }

  private val sampleInfosNoUnmatched = sampleInfos.filterNot(_.isUnmatched)
  private val unmatchedSample        = sampleInfos.find(_.isUnmatched).getOrElse(throw new IllegalArgumentException("No unmatched sample provided."))

  /** The number of reads that are expected to be given to the [[demultiplex()]] method. */
  def expectedNumberOfReads: Int = readStructures.length

  /** True if the read structure implies paired end reads will be produced, false otherwise. */
  val pairedEnd: Boolean = readStructures.count(_.template.nonEmpty) == 2

  /** Converts the segments across all reads to a string, with each read delimited by the given delimiter. */
  private def toBarcode(allBases: Seq[String], delimiter: String, allSegments: Seq[Seq[ReadSegment]]): String = {
    // Get the barcode bases delimited by read, not segment
    allBases.zip(allSegments).map { case (bases, segments) =>
      // Get the bases for the given read.  For example, we may have 8B2S4B100T.  If the segments are sample barcodes,
      // we want the bases corresponding to the 8B4S.
      ReadStructure.structureRead(bases=bases, segments=segments).map(_.bases).filter(_.nonEmpty).mkString
    }.filter(_.nonEmpty).mkString(delimiter)
  }

  /** The sub-reads for the sample barcodes across all read structures. */
  private val sampleBarcodeSegments = this.readStructures.map(_.sampleBarcode)

  /** Get the sample barcodes for all reads, with reads delimited by the given delimiter. */
  private def sampleBarcode(allBases: Seq[String], delimiter: String): String = toBarcode(allBases, delimiter, sampleBarcodeSegments)

  /** The sub-reads for the molecular barcodes across all read structures. */
  private val molecularBarcodeSegments = this.readStructures.map(_.molecularBarcode)

  /** Get the molecular barcode bases for all reads, with reads delimited by the given delimiter. */
  private def molecularBarcode(allBases: Seq[String], delimiter: String): String = toBarcode(allBases, delimiter, molecularBarcodeSegments)

  /** Gets the [[SampleInfo]] and the number of mismatches between the bases and matched sample barcode.  If no match is
    * found, the unmatched sample and [[Int.MaxValue]] are returned. */
  private def matchSampleBarcode(allBases: Seq[String]): (SampleInfo, Int) = {
    val observedBarcode = sampleBarcode(allBases, delimiter="").getBytes
    val numNoCalls      = observedBarcode.count(base => SequenceUtil.isNoCall(base))

    // Get the best and second best sample barcode matches.
    val (bestSampleInfo, bestMismatches, secondBestMismatches) = if (numNoCalls <= maxNoCalls) {
      this.sampleInfosNoUnmatched.map { sampleInfo =>
        val sample          = sampleInfo.sample
        val expectedBarcode = sample.sampleBarcodeBytes
        require(expectedBarcode.nonEmpty, s"Sample with id '${sample.sampleId}' did not have a sample barcode")
        val numMismatches   = countMismatches(observedBarcode, expectedBarcode)
        (sampleInfo, numMismatches)
      }
      .toList.sortBy(_._2).take(2) match {
        case Nil                              => (this.unmatchedSample, Int.MaxValue, Int.MaxValue)
        case List(bestTuple)                  => (bestTuple._1, bestTuple._2, Int.MaxValue)
        case List(bestTuple, secondBestTuple) => (bestTuple._1, bestTuple._2, secondBestTuple._2)
      }
    }
    else {
      (this.unmatchedSample, Int.MaxValue, Int.MaxValue)
    }

    // Make sure we are within the parameter limits and update barcode metrics if necessary
    if (maxMismatches < bestMismatches || maxNoCalls < numNoCalls || (secondBestMismatches - bestMismatches) < minMismatchDelta) {
      (this.unmatchedSample, Int.MaxValue)
    }
    else {
      (bestSampleInfo, bestMismatches)
    }
  }

  /** Demultiplexes a given set of reads from the same template.  The same number of reads should be given as read
    * structures.
    *
    * The sample barcoded bases from each read are extracted and concatenated in the same order as the given reads. They
    * are matched against the sample barcode bases for each sample.
    * */
  def demultiplex(reads: FastqRecord*): DemuxResult = {
    require(reads.nonEmpty, "No reads given for demultiplexing.")
    require(reads.length == expectedNumberOfReads, s"Expected '$expectedNumberOfReads' number of reads but found '${reads.length}'.")

    // Get the sample
    val allBases                    = reads.map(_.bases)
    val (sampleInfo, numMismatches) = matchSampleBarcode(allBases)

    // Create the [[DemuxRecord]]s.
    val molecularBarcode  = this.molecularBarcode(allBases, "-")
    val sampleBarcode     = this.sampleBarcode(allBases, "-")
    var readNumber        = 1
    val records           = reads.zip(readStructures).flatMap { case (read, readStructure) =>
      val templateSegments = readStructure.template
      if (templateSegments.isEmpty) {
        None
      }
      else {
        // NB: strict=false in case folks have trimmed the reads prior can be handled gracefully.
        val subReads: Seq[SubRead] = ReadStructure.structureReadWithQualities(bases=read.bases, qualities=read.quals, segments=templateSegments, strict=false)
        val bases                  = subReads.map(_.bases).mkString
        val quals                  = subReads.flatMap(_.quals).mkString
        require(bases.length == quals.length, s"Bases and qualities have differing lengths for read: ${read.header}")
        require(readNumber == 1 || readNumber == 2, s"Invalid read number '$readNumber' for read: ${read.header}")

        def f(s: String): Option[String] = if (s.isEmpty) None else Some(s)
        yieldAndThen(Some(DemuxRecord(read.header, bases, quals, f(molecularBarcode), f(sampleBarcode), readNumber, pairedEnd)))(readNumber += 1)
      }
    }

    DemuxResult(sampleInfo=sampleInfo, numMismatches=numMismatches, records=records)
  }
}
