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

import java.io.Closeable

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, unreachable}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{DirPath, FilePath, PathPrefix, PathToFastq}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, Logger}
import com.fulcrumgenomics.fastq.FastqDemultiplexer.{DemuxRecord, DemuxResult}
import com.fulcrumgenomics.illumina.{Sample, SampleSheet}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.ReadStructure.SubRead
import com.fulcrumgenomics.util.{ReadStructure, SampleBarcodeMetric, _}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.util.{Iso8601Date, SequenceUtil}
import enumeratum.EnumEntry

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable.ListBuffer
import scala.concurrent.forkjoin.ForkJoinPool

object DemuxFastqs {

  /** The name of the metrics file if none is given by the user. */
  val DefaultDemuxMetricsFileName = "demux_barcode_metrics.txt"

  /** The name of the sample for unmatched reads. */
  val UnmatchedSampleId: String = "unmatched"

  /** The maximum # of records in RAM per SAM/BAM writer. */
  private[fastq] val MaxRecordsInRam: Int = 5e6.toInt

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

  /** Gets the quality format of the FASTQs. If a format is given, checks that the given format is compatible. */
  private def determineQualityFormat(fastqs: Seq[PathToFastq],
                                     format: Option[QualityEncoding] = None,
                                     logger: Option[Logger] = None): QualityEncoding = {
    val detector = new QualityEncodingDetector
    detector.sample(fastqs.iterator.flatMap(FastqSource(_)).map(_.quals))

    format match {
      case Some(f) =>
        require(detector.isCompatible(f), s"Fastq is not compatible with provided quality encoding: $f")
        f
      case None =>
        val encs = detector.rankedCompatibleEncodings(q=30)
        require(encs.nonEmpty, "Could not determine quality score encoding in fastq. No known encodings are valid for all observed qualities.")
        if (encs.size > 1) logger.foreach(_.warning(s"Making ambiguous determination about fastq's quality encoding; possible encodings: ${encs.mkString(", ")}."))
        logger.foreach(_.info(s"Auto-detected quality format as: ${encs.head}"))
        encs.head
    }
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
    val barcodeSegments = readStructures.flatMap(_.sampleBarcodeSegments)
    require(barcodeSegments.nonEmpty, "No sample barcodes found in read structures: " + readStructures.mkString(", "))
    require(barcodeSegments.forall(_.hasFixedLength), "Barcode segments must have fixed lengths in: " + readStructures.mkString(", "))

    val noMatchBarcode: String = barcodeSegments.map("N" * _.fixedLength).mkString
    Sample(sampleOrdinal=sampleOrdinal, sampleId=UnmatchedSampleId, sampleName=UnmatchedSampleId, libraryId=UnmatchedSampleId, i7IndexBases=Some(noMatchBarcode))
  }

  private[fastq] def toSampleOutputPrefix(sample: Sample, isUnmatched: Boolean, illuminaStandards: Boolean, output: DirPath, unmatched: String): PathPrefix = {
    (isUnmatched, illuminaStandards) match {
      case (true, true)   => output.resolve(f"${UnmatchedSampleId}_S${sample.sampleOrdinal}_L${sample.lane.getOrElse(1)}%03d")
      case (false, true)  => output.resolve(f"${sample.sampleName}_S${sample.sampleOrdinal}_L${sample.lane.getOrElse(1)}%03d")
      case (true, false)  => PathUtil.removeExtension(output.resolve(unmatched))
      case (false, false) => outputPrefixFrom(output, sample)
    }
  }

  /** Creates a demultiplexing iterator that performs demultiplexing in parallel.
    *
    * @param sources the FASTQ sources, one per read.
    * @param demultiplexer the demultiplexer to use.  The demultiplexer's [[FastqDemultiplexer.demultiplex()]] method
    *                      expects the same number of reads as sources.
    */
  def demultiplexingIterator(sources: Seq[FastqSource],
                             demultiplexer: FastqDemultiplexer,
                             threads: Int,
                             batchSize: Int = DemuxBatchRecordsSize): Iterator[DemuxResult] = {

    require(demultiplexer.expectedNumberOfReads == sources.length,
      s"The demultiplexer expects ${demultiplexer.expectedNumberOfReads} reads but ${sources.length} FASTQ sources given.")

    val zippedIterator = FastqSource.zipped(sources)
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
      |The type of output is specified with the `--output-type` option, and can be BAM (`--output-type BamOnly`),
      |gzipped FASTQ (`--output-type FastqOnly`), or both (`--output-type Both`).
      |
      |For BAM output, the output directory will contain one BAM file per sample in the sample sheet or metadata CSV file,
      |plus a BAM for reads that could not be assigned to a sample given the criteria.  The output file names will be the
      |concatenation of sample id, sample name, and sample barcode bases (expected not observed), delimited by `-`.  A
      |metrics file will also be output providing analogous information to the metric described
      |[SampleBarcodeMetric](http://fulcrumgenomics.github.io/fgbio/metrics/latest/#samplebarcodemetric).
      |
      |For gzipped FASTQ output, one or more gzipped FASTQs per sample in the sample sheet or metadata CSV file will be
      |written to the output directory. For paired end data, the output will have the suffix `_R1.fastq.gz` and
      |`_R2.fastq.gz` for read one and read two respectively.  The sample barcode and molecular barcodes (concatenated)
      |will be appended to the read name and delimited by a colon.  If the `--illumina-standards` option is given, then
      |the output read names and file names will follow the
      |[Illumina standards described here](https://help.basespace.illumina.com/articles/tutorials/upload-data-using-web-uploader/).
      |
      |The output base qualities will be standardized to Sanger/SAM format.
      |
      |FASTQs and associated read structures for each read should be given:
      |
      |- a single fragment read should have one FASTQ and one read structure
      |- paired end reads should have two FASTQs and two read structures
      |- a dual-index sample with paired end reads should have four FASTQs and four read structures given: two for the
      |  two index reads, and two for the template reads.
      |
      |Read structures are made up of `<number><operator>` pairs much like the `CIGAR` string in BAM files. Four kinds of
      |operators are recognized:
      |
      |1. `T` identifies a template read
      |2. `B` identifies a sample barcode read
      |3. `M` identifies a unique molecular index read
      |4. `S` identifies a set of bases that should be skipped or ignored
      |
      |The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote "all remaining
      |bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length. Both reads must
      |have template bases.  Any molecular identifiers will be concatenated using
      |the `-` delimiter and placed in the given SAM record tag (`RX` by default).  Similarly, the sample barcode bases
      |from the given read will be placed in the `BC` tag.
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
      |concatenated with the sample barcode bases delimited by a `.`.
      |
      |The sample section of the sample sheet should contain information related to each sample with the following columns:
      |
      |```
      |  * Sample_ID:   The sample identifier unique to the sample in the sample sheet.
      |  * Sample_Name: The sample name.
      |  * Library_ID:  The library Identifier.  The combination sample name and library identifier should be unique
      |                 across the samples in the sample sheet.
      |  * Description: The description of the sample, which will be placed in the description field in the output BAM's
      |                 read group.  This column may be omitted.
      |```
      |
      |Additionally, the sample barcode should be specified in a column named `Sample_Barcode`.  The name of the column
      |containing the sample barcode can be changed using the `--column-for-sample-barcode` option.  If the sample barcode
      |is present across multiple reads (ex. dual-index, or inline in both reads of a pair), then the expected barcode
      |bases from each read should be concatenated and placed in the `Sample_Barcode` column.  The concatenation should
      |be in the same order as the order of the reads' FASTQs and read structures given to this tool.
      |
      |## Metadata CSV
      |
      |In lieu of a sample sheet, a simple CSV file may be provided with the necessary metadata.  This file should
      |contain the same columns as described above for the sample sheet (`Sample_ID`, `Sample_Name`, `Library_ID`, and
      |`Description`).
      |
      |## Example Command Line
      |
      |As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both reading a sample
      |barcode, as well as an in-line 8bp sample barcode in read one, the command line would be
      |
      |```
      |--inputs r1.fq i1.fq i2.fq r2.fq --read-structures 8B92T 8B 8B 100T \
      |    --sample-sheet SampleSheet.csv --metrics metrics.txt --output output_folder
      |```
      |
      |## Illumina FASTQ File Naming
      |
      |The output FASTQs will be written according to Illumina naming standards when the `--illumina-standards` option
      |is specified.  This requires that the input fastq files contain Illumina format read names with comments containing
      |`<ReadNum>:<FilterFlag>:<0>:<SampleNumber>` (no angle brackets).  The output create output files with Illumina
      |names:
      |
      |1. The file extension will be `_R1_001.fastq.gz` for read one, and `_R2_001.fastq.gz` for read two (if paired end).
      |2. The per-sample output prefix will be `<SampleName>_S<SampleOrdinal>_L<LaneNumber>` (without angle brackets).
    """,
  group=ClpGroups.Fastq
)
class DemuxFastqs
(@arg(flag='i', doc="One or more input fastq files each corresponding to a sub-read (ex. index-read, read-one, read-two, fragment).") val inputs: Seq[PathToFastq],
 @arg(flag='o', doc="The output directory in which to place sample BAMs.") val output: DirPath,
 @arg(flag='x', doc="A file containing the metadata about the samples.") val metadata: FilePath,
 @arg(flag='r', doc="The read structure for each of the FASTQs.") val readStructures: Seq[ReadStructure],
 @arg(flag='m', doc="The file to which per-barcode metrics are written.  If none given, a file named `demux_barcode_metrics.txt` will be written to the output directory.") val metrics: Option[FilePath] = None,
 @arg(flag='c', doc="The column name in the sample sheet or metadata CSV for the sample barcode.") val columnForSampleBarcode: String = "Sample_Barcode",
 @arg(flag='u', doc="Output BAM file name for the unmatched records.") val unmatched: String = DemuxFastqs.UnmatchedSampleId + ".bam",
 @arg(flag='q',
   doc="""A value describing how the quality values are encoded in the FASTQ.
      |Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66),
      |Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard
      |for phred scaled scores with a character shift of 33.  If this value
      |is not specified, the quality format will be detected automatically.
    """)
 val qualityFormat: Option[QualityEncoding] = None,
 @arg(flag='t', doc="The number of threads to use while de-multiplexing. The performance does not increase linearly with the # of threads and seems not to improve beyond 2-4 threads.") val threads: Int = 1,
 @arg(doc="Maximum mismatches for a barcode to be considered a match.") val maxMismatches: Int = 1,
 @arg(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.") val minMismatchDelta: Int = 2,
 @arg(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.") val maxNoCalls: Int = 2,
 @arg(doc="The sort order for the output sam/bam file (typically unsorted or queryname).") val sortOrder: SortOrder = SortOrder.queryname,
 @arg(doc="The SAM tag for any molecular barcode.  If multiple molecular barcodes are specified, they will be concatenated and stored here.") val umiTag: String = ConsensusTags.UmiBases,
 @arg(doc="The platform unit (typically `<flowcell-barcode>-<sample-barcode>.<lane>`)") val platformUnit: Option[String] = None,
 @arg(doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
 @arg(doc="Predicted median insert size, to insert into the read group header") val predictedInsertSize: Option[Integer] = None,
 @arg(doc="Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)") val platformModel: Option[String] = None,
 @arg(doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
 @arg(doc="Date the run was produced, to insert into the read group header") val runDate: Option[Iso8601Date] = None,
 @arg(doc="The type of outputs to produce.", mutex=Array("outputFastqs")) var outputType: Option[OutputType] = None,
 @deprecated("Use outputType instead.", since="0.5.0")
 @arg(doc="*** Deprecated: use --output-type instead ***. Output gzipped FASTQs (`.fastq.gz`) instead of BAM files", mutex=Array("outputType")) val outputFastqs: Option[Boolean] = None,
 @arg(doc="Output FASTQs according to Illumina naming standards, for example, for upload to the BaseSpace Sequence Hub") val illuminaStandards: Boolean = false,
 @arg(
   doc=
     """Output all bases (i.e. all sample barcode, molecular barcode, skipped,
        and template bases) for every read with template bases (ex. read one
        and read two) as defined by the corresponding read structure(s).
     """)
 val includeAllBasesInFastqs: Boolean = false
) extends FgBioTool with LazyLogging {

  import DemuxFastqs._

  private[fastq] val metricsPath = metrics.getOrElse(output.resolve(DefaultDemuxMetricsFileName))

  // NB: remove me when outputFastqs gets removed and use outputType directly
  private val _outputType = (this.outputType, this.outputFastqs) match {
    case (None, None) => OutputType.Bam
    case (None, Some(fastqs)) => if (fastqs) OutputType.Fastq else OutputType.Bam
    case (Some(tpe), None) => tpe
    case _ => unreachable("Bug: outputType and outputFastqs should never be both defined")
  }

  validate(inputs.length == readStructures.length, "The same number of read structures should be given as FASTQs.")
  validate(readStructures.flatMap(_.sampleBarcodeSegments).nonEmpty, s"No sample barcodes found in read structures: " + readStructures.map(_.toString).mkString(", "))

  private val pairedEnd = readStructures.count(_.templateSegments.nonEmpty) match {
    case 1 => false
    case 2 => true
    case n => invalid(s"Found $n read structures with template bases but expected 1 or 2.")
  }

  if (illuminaStandards) {
    validate(this._outputType != OutputType.Bam, s"--illumina-standards may only be used with '--output-type ${OutputType.Fastq}' and '--output-type ${OutputType.BamAndFastq}'.")
    validate(pairedEnd,  "--illumina-standards may only be used with paired end data")
  }

  Io.assertReadable(inputs)
  Io.assertReadable(metadata)

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
    Io.mkdirs(this.output)
    Io.assertCanWriteFile(this.metricsPath)

    // Get the FASTQ quality encoding format
    val qualityEncoding: QualityEncoding = determineQualityFormat(fastqs=inputs, format=this.qualityFormat, logger=Some(this.logger))

    // Read in the sample sheet and create the sample information
    val samplesFromSampleSheet = sampleSheet.map(s => withCustomSampleBarcode(s, columnForSampleBarcode))
    val samples                = samplesFromSampleSheet.toSeq :+ unmatchedSample(samplesFromSampleSheet.size, this.readStructures)
    val sampleInfos            = samples.map(sample => SampleInfo(sample, sample.sampleName == UnmatchedSampleId))
    val sampleToWriter         = sampleInfos.map { info => info.sample -> toWriter(info, sampleInfos.length) }.toMap

    // Validate that the # of sample barcode bases in the read structure matches the # of sample barcode in the sample sheet.
    {
      val rsNumSampleBarcodeBases = readStructures.map(_.sampleBarcodeSegments.map(_.fixedLength).sum).sum
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
      maxNoCalls       = maxNoCalls,
      includeOriginal  = this.includeAllBasesInFastqs
    )

    val progress = ProgressLogger(this.logger, unit=1e6.toInt)

    // An iterator that uses the given fastq demultiplexer to convert FASTQ records from the same fragment/template to
    // DemuxRecord in parallel
    val sources   = inputs.map(FastqSource(_))
    val iterator  = demultiplexingIterator(
      sources       = sources,
      demultiplexer = demultiplexer,
      threads       = threads
    )

    // Write the records out in its own thread
    iterator.foreach { demuxResult =>
      demuxResult.sampleInfo.metric.increment(numMismatches=demuxResult.numMismatches)
      val writer = sampleToWriter(demuxResult.sampleInfo.sample)
      demuxResult.records.foreach { rec =>
        writer.add(rec.copy(quals=qualityEncoding.toStandardAscii(rec.quals)))
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

  private def toWriter(sampleInfo: SampleInfo, numSamples: Int): DemuxWriter = {
    val sample = sampleInfo.sample
    val isUnmatched = sample.sampleName == UnmatchedSampleId
    val prefix = toSampleOutputPrefix(sample, isUnmatched, illuminaStandards, output, this.unmatched)

    val writers = new ListBuffer[DemuxWriter]()

    if (this._outputType.producesFastq) {
      writers += new FastqRecordWriter(prefix, this.pairedEnd, illuminaStandards)
    }

    if (this._outputType.producesBam) {
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

        writers += new SamRecordWriter(prefix, header, this.umiTag, numSamples)
    }

    new DemuxWriter {
      def add(rec: DemuxRecord): Unit = writers.foreach { writer => writer.add(rec) }
      override def close(): Unit = writers.foreach(_.close())
    }
  }
}

/** A writer than writes [[DemuxRecord]]s */
private trait DemuxWriter extends Closeable {
  def add(rec: DemuxRecord): Unit
}

/** A writer that writes [[DemuxRecord]]s as [[SamRecord]]s. */
private class SamRecordWriter(prefix: PathPrefix,
                              val header: SAMFileHeader,
                              val umiTag: String,
                              val numSamples: Int) extends DemuxWriter {
  val order: Option[SamOrder] = if (header.getSortOrder == SortOrder.unsorted) None else SamOrder(header)
  private val writer = SamWriter(PathUtil.pathTo(prefix + ".bam"), header, sort=order,
    async = DemuxFastqs.UseAsyncIo,
    maxRecordsInRam = Math.max(10000,  DemuxFastqs.MaxRecordsInRam / numSamples))

  private val rgId: String = this.header.getReadGroups.get(0).getId

  def add(rec: DemuxRecord): Unit = {
    val record = SamRecord(header)
    record.name     = rec.name
    record.bases    = rec.bases
    record.quals    = rec.quals
    record.unmapped = true
    if (rec.pairedEnd) {
      record.paired = true
      record.mateUnmapped = true
      if (rec.readNumber == 1) record.firstOfPair = true else record.secondOfPair = true
    }
    rec.sampleBarcode.foreach(bc => record("BC") = bc)
    record(ReservedTagConstants.READ_GROUP_ID) =  rgId
    rec.molecularBarcode.foreach(mb => record(umiTag) = mb)
    writer += record
  }

  override def close(): Unit = writer.close()
}

private[fastq] object FastqRecordWriter {
  private[fastq] def extensions(pairedEnd: Boolean, illuminaStandards: Boolean = false): Seq[String] = (pairedEnd, illuminaStandards) match {
    case (true, true)   => Seq("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    case (true, false)  => Seq("_R1.fastq.gz", "_R2.fastq.gz")
    case (false, true)  => Seq("_R1_001.fastq.gz")
    case (false, false) => Seq(".fastq.gz")
  }
}

/** A writer that writes [[DemuxRecord]]s as [[FastqRecord]]s. */
private[fastq] class FastqRecordWriter(prefix: PathPrefix, val pairedEnd: Boolean, val illuminaStandards: Boolean = false) extends DemuxWriter {
  private val writers: IndexedSeq[FastqWriter] = FastqRecordWriter.extensions(pairedEnd=pairedEnd, illuminaStandards=illuminaStandards).map { ext =>
    FastqWriter(Io.toWriter(PathUtil.pathTo(prefix + ext)))
  }.toIndexedSeq

  private[fastq] def readName(rec: DemuxRecord): String = {
    if (illuminaStandards) {
      val comment = rec.comment.getOrElse {
        throw new IllegalArgumentException(s"Comment required with illumina-standards and read: '${rec.name}'")
      }

      require(rec.name.count(_ == ':') == 6,
        s"Expected the read name format 'Instrument:RunID:FlowCellID:Lane:Tile:X:Y', found '${rec.name}'")
      require(comment.count(_ == ':') == 3,
        s"Expected the comment format 'ReadNum:FilterFlag:0:SampleNumber', found '${rec.name}'")

      s"${rec.name} $comment"
    }
    else {
      Seq(Some(rec.name), rec.sampleBarcode, rec.molecularBarcode).flatten.mkString(":")
    }
  }

  def add(rec: DemuxRecord): Unit = {
    val record = FastqRecord(
      name       = readName(rec),
      bases      = rec.originalBases.getOrElse(rec.bases),
      quals      = rec.originalQuals.getOrElse(rec.quals),
      comment    = None,
      readNumber = if (illuminaStandards) None else Some(rec.readNumber)
    )

    val writer = this.writers.lift(rec.readNumber-1).getOrElse {
      throw new IllegalStateException(s"Read number was invalid: ${rec.readNumber}")
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
  case class DemuxRecord(name: String, bases: String, quals: String, molecularBarcode: Option[String],
                         sampleBarcode: Option[String], readNumber: Int, pairedEnd: Boolean, comment: Option[String],
                         originalBases: Option[String] = None, originalQuals: Option[String] = None)

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
  * @param includeOriginal true if to provide set the values for `originaBases` and `originalQuals` in [[DemuxResult]],
  *                        namely the bases and qualities FOR ALL bases, including molecular barcode, sample barcode,
  *                        and skipped bases.
  */
private class FastqDemultiplexer(val sampleInfos: Seq[SampleInfo],
                                 readStructures: Seq[ReadStructure],
                                 val umiTag: String = ConsensusTags.UmiBases,
                                 val maxMismatches: Int = 2,
                                 val minMismatchDelta: Int = 1,
                                 val maxNoCalls: Int = 2,
                                 val includeOriginal: Boolean = false) {
  import FastqDemultiplexer._

  require(readStructures.nonEmpty, "No read structures were given")
  private val variableReadStructures = readStructures.map(_.withVariableLastSegment)

  {
    val samples = sampleInfos.map(_.sample)
    require(samples.map(_.sampleBarcodeString).sorted.distinct.length == samples.length, "Unique sample barcodes required for all samples")
  }

  private val sampleInfosNoUnmatched = sampleInfos.filterNot(_.isUnmatched)
  private val unmatchedSample        = sampleInfos.find(_.isUnmatched).getOrElse(throw new IllegalArgumentException("No unmatched sample provided."))

  /** The number of reads that are expected to be given to the [[demultiplex()]] method. */
  def expectedNumberOfReads: Int = this.variableReadStructures.length

  /** True if the read structure implies paired end reads will be produced, false otherwise. */
  val pairedEnd: Boolean = this.variableReadStructures.count(_.templateSegments.nonEmpty) match {
    case 0 => throw new IllegalArgumentException("No template reads in any read structure.")
    case 1 => false
    case 2 => true
    case n => throw new IllegalArgumentException(s"$n template reads defined. Can't process > 2 template reads.")
  }

  /** Gets the [[SampleInfo]] and the number of mismatches between the bases and matched sample barcode.  If no match is
    * found, the unmatched sample and [[Int.MaxValue]] are returned. */
  private def matchSampleBarcode(subReads: Seq[SubRead]): (SampleInfo, Int) = {
    val observedBarcode = subReads.filter(_.kind == SegmentType.SampleBarcode).map(_.bases).mkString.getBytes
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

    // Generate the sub-reads by type
    val subReads = reads.zip(this.variableReadStructures).flatMap { case (read, rs) => rs.extract(read.bases, read.quals) }

    // Get the sample
    val (sampleInfo, numMismatches) = matchSampleBarcode(subReads)

    // Method to get all the bases of a given type
    def bases(segmentType: SegmentType): Option[String] = {
      val b = subReads.filter(_.kind == segmentType).map(_.bases).mkString("-")
      if (b.isEmpty) None else Some(b)
    }

    // Get the molecular and sample barcodes
    val molecularBarcode = bases(SegmentType.MolecularBarcode)
    val sampleBarcode    = bases(SegmentType.SampleBarcode)

    val demuxRecords = reads.zip(this.variableReadStructures)
      .filter { case (_, rs) => rs.exists(_.kind == SegmentType.Template) }
      .zipWithIndex
      .map { case ((read, rs), readIndex) =>
        val segments   = rs.extract(read.bases, read.quals)
        val readNumber = readIndex + 1
        val templates  = segments.filter(_.kind == SegmentType.Template)
        require(templates.nonEmpty, s"Bug: require at least one template in read $readIndex; read structure was ${segments.mkString}")
        DemuxRecord(
          name             = read.name,
          bases            = templates.map(_.bases).mkString,
          quals            = templates.map(_.quals).mkString,
          molecularBarcode = molecularBarcode,
          sampleBarcode    = sampleBarcode,
          readNumber       = readNumber,
          pairedEnd        = this.pairedEnd,
          comment          = read.comment,
          originalBases    = if (this.includeOriginal) Some(read.bases) else None,
          originalQuals    = if (this.includeOriginal) Some(read.quals) else None
        )
      }

    DemuxResult(sampleInfo=sampleInfo, numMismatches=numMismatches, records=demuxRecords)
  }
}

sealed trait OutputType extends EnumEntry {
  def producesBam: Boolean
  def producesFastq: Boolean
}
object OutputType extends FgBioEnum[OutputType] {
  def values: IndexedSeq[OutputType] = findValues
  case object Fastq extends OutputType { val producesBam: Boolean = false; val producesFastq: Boolean = true; }
  case object Bam extends OutputType { val producesBam: Boolean = true; val producesFastq: Boolean = false; }
  case object BamAndFastq extends OutputType { val producesBam: Boolean = true; val producesFastq: Boolean = true; }
}