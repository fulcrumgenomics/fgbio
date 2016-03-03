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
 *
 */
package com.fulcrumgenomics.fastq

import java.io.File
import java.nio.file.Path
import java.util.NoSuchElementException

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util._
import com.fulcrumgenomics.util.miseq.{Sample, SampleBarcode, SampleSheet}
import dagr.commons.CommonsDef.{DirPath, FilePath, PathToFastq, _}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.util._

import scala.collection.mutable

object DemuxFastqs {
  import SampleBarcode.SampleBarcodeDelimiter

  val UnmatchedSampleId: String = "unmatched"

  /** Custom columns in the Illumina Experiment Manager Sample Sheet. */
  object CustomSampleSheetHeaderNames {
    val R1SampleBarcodeBases = "R1_Barcode_Bases"
    val R2SampleBarcodeBases = "R2_Barcode_Bases"
  }

  /** Creates the sample output BAM for the given sample.  [setSampleBarcode] should be called first. */
  def sampleOutputBam(output: Path, sample: Sample): Path = {
    val sampleBarcode: String = sample.sampleBarcode.getOrElse(throw new IllegalArgumentException("Sample barcode missing for sample: " + sample.sampleName)).toString
    output.resolve(IOUtil.makeFileNameSafe(String.format("%s.%s.%s.bam", sample.sampleId, sample.sampleName, sampleBarcode)))
  }

  private val solexaQualityConverter: SolexaQualityConverter = SolexaQualityConverter.getSingleton

  /** A custom sample class to support inline sample barcodes for read one and read two. */
  private[fastq] class CustomSample(sample: Sample)
    extends Sample(
      sampleOrdinal      = sample.sampleOrdinal,
      sampleId           = sample.sampleId,
      sampleName         = sample.sampleName,
      libraryId          = sample.libraryId,
      project            = sample.project,
      description        = sample.description,
      lane               = sample.lane,
      i7IndexBases       = sample.i7IndexBases,
      i5IndexBases       = sample.i5IndexBases,
      sampleBarcode      = sample.sampleBarcode,
      extendedAttributes = sample.extendedAttributes
    ) {
    val r1SampleBarcodeBases: Option[String] = extendedAttribute(CustomSampleSheetHeaderNames.R1SampleBarcodeBases)
    val r2SampleBarcodeBases: Option[String] = extendedAttribute(CustomSampleSheetHeaderNames.R2SampleBarcodeBases)

    /** Returns the sample barcodes in order of sequencing. */
    override def sampleBarcodeBases: Seq[Option[String]] = {
      Seq(i7IndexBases, i5IndexBases, r1SampleBarcodeBases, r2SampleBarcodeBases)
    }
  }

  /** Counts the nucleotide mismatches between two strings of the same length.  Ignores no calls in expectedBases.
    * Observed base qualities less than the minimum base quality are counted as mismatches if not a no call. Observed
    * qualities may not be null. */
  private[fastq] def countMismatches(observedBases: Array[Byte], expectedBases: Array[Byte]): Int = {
    observedBases.zip(expectedBases).count {
      case (observedBase, expectedBase) => !SequenceUtil.isNoCall(expectedBase) && !SequenceUtil.basesEqual(observedBase, expectedBase)
    }
  }

  /** Reads multiple FASTQs at the same time, ensuring their read names match along the way.  This is
    * useful for when we don't have a single interleaved FASTQ. */
  private[fastq] class ParallelFastqReader(fastq: Seq[Path], fastqs: Seq[Path]*) extends Iterator[Seq[FastqRecord]] {

    val readers: Seq[SequentialFastqSource] = (fastq +: fastqs.toSeq).map { otherFastqs: Seq[Path] =>
      if (otherFastqs.size != fastq.size) throw new IllegalArgumentException("List of list fastqs had differing lengths.")
      SequentialFastqSource(otherFastqs)
    }

    def next: Seq[FastqRecord] = {
      val records = readers.map(_.next())
      verifySameReadNames(records)
      records
    }

    def hasNext: Boolean = {
      if (readersHaveNext(true)) true
      else if (readersHaveNext(false)) false
      else throw new IllegalArgumentException("Fastqs have differing lengths")
    }

    private def readersHaveNext(desiredValue: Boolean): Boolean = readers.forall(_.hasNext == desiredValue)

    private def verifySameReadNames(records: Seq[FastqRecord]): Unit = {
      val readName = records.head.name
      if (records.exists(_.name != readName)) {
        throw new IllegalArgumentException(String.format("Mismatching read names in FASTQS:\n%s\n%s\n", readName, records.find(_.name != readName).get))
      }
    }

    def close() = readers.foreach(_.close())
  }

  /** Helper case class to store the result of converting SAM records. */
  private case class AssignedRecords(sampleOrdinal: Int, records: Seq[SAMRecord])

  /** Helper class to convert fastq records into SAM records and write them out.
    *
    * @param input the iterator through groups of FASTQ records (r1, r2, i7, and i5).
    * @param headers the headers, one per sample
    * @param sampleSheet the sample sheet describing all the samples
    * @param readStructures the read structures info, one for each expected fastq record that will be passed
    *                           into `convertFastqsAndWriteRecords`.
    */
  private class FastqConverter(val input: Iterator[Seq[FastqRecord]],
                               val headers: Seq[SAMFileHeader],
                               val sampleSheet: SampleSheet,
                               val umiTag: String,
                               val multiUmiTags: Seq[String],
                               val qualityFormat: FastqQualityFormat,
                               val minQ: Int,
                               val maxQ: Int,
                               val maxMismatches: Int,
                               val minMismatchDelta: Int,
                               val maxNoCalls: Int,
                               val readStructures: ReadStructure*) extends Iterator[AssignedRecords] {
    /** The barcode string for reads that do not match. */
    private val noMatchBarcode: String = {
      readStructures.flatMap {
        readStructure => // get the sample barcode read structure for this read.
          val rs = readStructure.sampleBarcodeStructure(resetOffsets=true)
          if (rs.isEmpty) None
          else Some(rs)
        }
        .map { readStructure =>
          StringUtil.repeatCharNTimes('N', readStructure.map(_.length).sum)
        }.mkString(SampleBarcodeDelimiter)
    }

    /** Read structures for the sample barcodes, one read structure per read (ex. r1, r2, i7, and i5). */
    private val sampleBarcodeReadStructures = readStructures.map { readStructure =>
     readStructure.sampleBarcodeStructure(resetOffsets=false)
    }

    /** Read structures for the molecular barcodes, one read structure per read (ex. r1, r2, i7, and i5). */
    private val molecularBarcodeReadStructures = readStructures.map { readStructure =>
      readStructure.molecularBarcodeStructure(resetOffsets=false)
    }

    /** Sample barcode to its associate metric */
    private val barcodeToMetrics: Map[String, SampleBarcodeMetric] = {
      val metricMap = new mutable.LinkedHashMap[String, SampleBarcodeMetric]()
      sampleSheet.foreach { sample =>
        val barcode: String = sample.sampleBarcode match {
          case Some(b) => b.concatenatedBarcode
          case None    => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a sampleBarcode")
        }
        val libraryId: String = sample.libraryId match {
          case Some(id) => id
          case None     => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a library id")
        }
        metricMap.put(barcode, SampleBarcodeMetric(barcodeName=sample.sampleName, libraryName=libraryId, barcode=barcode))
      }
      metricMap.put(noMatchBarcode, SampleBarcodeMetric(UnmatchedSampleId, UnmatchedSampleId, noMatchBarcode))
      metricMap.toMap
    }

    def hasNext: Boolean = {
      false
    }

    def next(): AssignedRecords = {
      if (!hasNext) throw new NoSuchElementException("Calling next() when hasNext is false")

      val fastqRecords                        = input.next()
      val sampleBarcodeReadBases: String      = getSampleBarcode(fastqRecords).mkString
      val sampleOrdinal: Int                  = getSampleOrdinalFromSampleBarcode(sampleBarcodeReadBases, sampleSheet).getOrElse(headers.length)
      val header: SAMFileHeader               = headers(sampleOrdinal-1)
      val sampleId: String                    = header.getReadGroups.get(0).getId
      val molecularBarcodes: Seq[String]      = getMolecularBarcodes(fastqRecords)
      val pairedEnd: Boolean                  = readStructures.count(_.template.nonEmpty) > 1

      val records = fastqRecords.zipWithIndex.flatMap { case (record, recordIndex) =>
        val readStructure = readStructures(recordIndex)
        if (readStructure.template.isEmpty) {
          None
        }
        else {
          if (1 < recordIndex) throw new IllegalArgumentException("Found more than two ends with template bases")
          val samRecord: SAMRecord = makeSAMRecord(header, record, readStructure, sampleId, pairedEnd, recordIndex == 0)
          
          // set attributes: sample and molecular barcodes
          samRecord.setAttribute("BC", sampleBarcodeReadBases)
          val umiValue = molecularBarcodes.mkString(SampleBarcodeDelimiter)
          if (umiValue.nonEmpty) samRecord.setAttribute(umiTag, molecularBarcodes.mkString(SampleBarcodeDelimiter))

          if (multiUmiTags.nonEmpty) {
            val numMolecularBarcodeSegments = getNumMolecularBarcodes(fastqRecords)
            if (numMolecularBarcodeSegments == multiUmiTags.length) {
              multiUmiTags.zip(molecularBarcodes).foreach { case (tag, barcode) => samRecord.setAttribute(tag, barcode) }
            }
            else if (multiUmiTags.length == 1) {
              samRecord.setAttribute(multiUmiTags.head, molecularBarcodes.mkString(SampleBarcodeDelimiter))
            }
            else { // more than one tag, but not the same as the number of molecular barcodes
              unreachable(s"Did not find only one SAM tag or the same # of SAM tags for the molecular barcodes (${multiUmiTags.length} != $numMolecularBarcodeSegments)")
            }
          }
          Some(samRecord)
        }
      }
      AssignedRecords(sampleOrdinal=sampleOrdinal, records=records)
    }

    def getBarcodeMetrics: Iterable[SampleBarcodeMetric] = {
      SampleBarcodeMetric.finalizeMetrics(this.barcodeToMetrics, this.noMatchBarcode)
      this.barcodeToMetrics.values
    }

    private def noMatchBarcodeMetric: SampleBarcodeMetric = this.barcodeToMetrics(noMatchBarcode)

    /** Searches for a matching sample given the observed barcode.  Returns -1 if no match is found, given the
      * sample sheet and matching parameters */
    private def getSampleOrdinalFromSampleBarcode(sampleBarcodeReadBases: String, sampleSheet: SampleSheet): Option[Int] = {
      val numNoCalls: Int = sampleBarcodeReadBases.getBytes.count(base => SequenceUtil.isNoCall(base))
      val (bestSampleOrdinal, bestMismatches, secondBestMismatches) = if (numNoCalls <= maxNoCalls) {
        sampleSheet.map { sample =>
          val barcodeBytes = sample.sampleBarcode match {
            case Some(barcode) => barcode.barcodeBytes
            case None => throw new IllegalArgumentException(s"Sample barcode required for sample with id '${sample.sampleId}'")
          }
          (sample.sampleOrdinal, countMismatches(sampleBarcodeReadBases.getBytes, barcodeBytes))
        }.toList.sortBy(_._2).take(2) match {
          case Nil                              => (-1,           Integer.MAX_VALUE, Integer.MAX_VALUE)
          case List(bestTuple)                  => (bestTuple._1, bestTuple._2,      Integer.MAX_VALUE)
          case List(bestTuple, secondBestTuple) => (bestTuple._1, bestTuple._2,      secondBestTuple._2)
        }
      }
      else {
        (-1, Integer.MAX_VALUE, Integer.MAX_VALUE)
      }
      // Make sure we are within the parameter limits and update barcode metrics if necessary
      if (maxMismatches < bestMismatches || maxNoCalls < numNoCalls || (secondBestMismatches - bestMismatches) < minMismatchDelta) {
        noMatchBarcodeMetric.reads +=  1
        noMatchBarcodeMetric.pf_reads += 1
        None
      }
      else {
        val sampleBarcode = sampleSheet.get(bestSampleOrdinal - 1).sampleBarcode.getOrElse(unreachable(s"Sample barcode required for sample"))
        val bestBarcodeMetric: SampleBarcodeMetric = this.barcodeToMetrics.getOrElse(sampleBarcode.concatenatedBarcode, unreachable("No metric for sample"))
        bestBarcodeMetric.reads += 1
        bestBarcodeMetric.pf_reads += 1
        if (bestMismatches == 0) {
          bestBarcodeMetric.perfect_matches += 1
          bestBarcodeMetric.pf_perfect_matches += 1
        }
        else if (bestMismatches == 1) {
          bestBarcodeMetric.one_mismatch_matches += 1
          bestBarcodeMetric.pf_one_mismatch_matches += 1
        }
        Some(bestSampleOrdinal)
      }
    }

    /** Creates a SAM record (always paired) */
    private def makeSAMRecord(header: SAMFileHeader, fastqRecord: FastqRecord, readStructure: ReadStructure,
                              readGroupId: String, pairedEnd: Boolean, firstOfPair: Boolean): SAMRecord = {
      val (readbases, readQualities) = readStructure.structureReadWithQualities(
        bases     = fastqRecord.bases,
        qualities = fastqRecord.quals,
        strict=false
      ).map { case (b, q, rs) => (b, q) }.unzip
      val record: SAMRecord = new SAMRecord(header)
      val readBases = readbases.mkString.getBytes
      val baseQualities = readQualities.mkString.getBytes

      record.setReadName(fastqRecord.name)
      record.setReadBases(readBases)
      record.setReadUnmappedFlag(true)
      if (pairedEnd) {
        record.setReadPairedFlag(true)
        record.setMateUnmappedFlag(true)
        if (firstOfPair) record.setFirstOfPairFlag(true)
        else record.setSecondOfPairFlag(true)
      }
      record.setAttribute(ReservedTagConstants.READ_GROUP_ID, readGroupId)

      convertQuality(baseQualities, qualityFormat)
      baseQualities.foreach {
        qual =>
          val uQual: Int = qual & 0xff
          if (uQual < minQ || uQual > maxQ) {
            throw new IllegalStateException(s"Base quality $uQual is not in the range $minQ ... $maxQ for read ${fastqRecord.header}")
          }
      }
      record.setBaseQualities(baseQualities)

      record
    }

    /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
    private[fastq] def convertQuality(quals: Array[Byte], version: FastqQualityFormat) {
      version match {
        case FastqQualityFormat.Standard => SAMUtils.fastqToPhred(quals)
        case FastqQualityFormat.Solexa => solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals)
        case FastqQualityFormat.Illumina => solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals)
      }
    }

    /** Gets the sample barcode bases */
    private def getSampleBarcode(records: Seq[FastqRecord]): Seq[String] = {
      records.zip(sampleBarcodeReadStructures).map {
        case (rec, rs) => rs.structureRead(rec.bases).map(_._1).mkString
      }
    }

    private def getMolecularBarcodes(records: Seq[FastqRecord]): Seq[String] = {
      records.zip(molecularBarcodeReadStructures).map {
        case (rec, rs) => rs.structureRead(rec.bases).map(_._1).mkString
      }
    }

    private def getNumMolecularBarcodes(records: Seq[FastqRecord]): Int = {
      records.zip(readStructures).map {
        case (rec, rs) => rs.molecularBarcode.length
      }.sum
    }
  }
}

/**
  * This tool de-multiplexes a set of FASTQs based on the given Illumina Experiment Manager Sample Sheet for dual-indexed
  * sequencing runs.
  *
  * See the USAGE for a detailed description.
  *
  * Possible Future Improvements:
  * - adapter trimming
  * - more metadata stored in the SampleSheet.csv
  */
@clp(
  description = """|Demultiplexes FASTQs based on the given Illumina Experiment Manager Sample Sheet.
                   |
                   |Fastqs and read structures for read one, read two, i7 read, and i5 read should be given.  The read structures
                   |may contain sample barcode bases ('B'), molecular identifier bases ('M'), template bases ('T'), and bases to
                   |skip ('S'). Template bases should only be found on read one or read two.  Any molecular identifiers will be
                   |concatenated using the '-' delimiter and placed in the given SAM record tag ("RX" by default).  The order of
                   |concatenation will be read one, read two, i7 read, and i5 read, only considering reads with molecular identifiers.
                   |Similarly, the sample barcode bases from the given read will be placed in the "BC" tag, using the same rules as
                   |molecular identifiers, but applied to sample barcodes.
                   |
                   |The sample barcode for each sample in the sample sheet will be compared against the sample barcode bases in
                   |each read, to assign each read to a sample.  Reads that to not match any sample within the given error tolerance
                   |will be placed in the 'unmatched' file.
                   |
                   |The output directory will contain one BAM file per sample in the sample sheet, plus a BAM for reads that could
                   |not be assigned to a sample given the criteria.  The output files will be the concatenation of sample id, sample
                   |name, and sample barcode bases (expected not observed), delimited by ".".  A metrics file will also be output
                   |for sample barcodes.  More information about these metrics can be found here:
                   |  https://broadinstitute.github.io/picard/picard-metric-definitions.html#SampleBarcodeMetric
                   |
                   |The read group's sample id, sample name, and library id all correspond to the similarly named values in the sample sheet.
                   |Library id will be the sample id if not found, and the platform unit will be the sample name concatenated with the sample
                   |barcode bases delimited by a ".".
                   |
                   |The sample section of the sample sheet should contain information related to each sample with the following keys and values:
                   |  - Sample Identifier:  Sample_ID
                   |  - Sample Name:        Sample_Name
                   |  - Library Identifier: Library_ID
                   |  - Sample Project:     Sample_Project
                   |  - Description:        Description
                   |
                   |The following are optional values may included in the sample sheet to specify sample barcode bases in any read:
                   |  - Read One Inline Sample Barcode Bases*: R1_Barcode_Bases
                   |  - Read Two Inline Sample Barcode Bases*: R2_Barcode_Bases
                   |  - i7 Sample Barcode Bases:              Index
                   |  - i5 Sample Barcode Bases:              Index2
                   |* Note that these are custom extensions and not standard.
                   |
                   |The read structures will be used to extract the observed sample barcode and molecular identifiers from each
                   |read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in the sample sheet
                   |and associated read structures.  Please see the following link for the read structure:
                   |  https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/parser/ReadStructure.html
                   |
                   |""",
  group = ClpGroups.Utilities)
class DemuxFastqs
(
  @arg(flag="f1", doc="Input fastq file (optionally gzipped) for the first read of paired end data.")
  val fastq1: List[PathToFastq],
  @arg(flag="f2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", minElements = 0)
  val fastq2: List[PathToFastq] = List.empty,
  @arg(flag="i7", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i7 sequencing primer. This is typically the I1 FASTQ (index read one for dual-indexed data).", minElements = 0)
  val fastq7: List[PathToFastq] = List.empty,
  @arg(flag="i5", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i5 sequencing primer. This is typically the I2 FASTQ (index read two for dual-indexed data).", minElements = 0)
  val fastq5: List[PathToFastq] = List.empty,
  @arg(doc="Read structure for the first read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
  val rs1: String = "1000T",
  @arg(doc="Read structure for the second read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
  val rs2: Option[String] = None,
  @arg(doc="Read structure for the index read of the Illumina i7 sequencing primer.  The total bases must match the length of the index read. See the DemuxFastqs help message for more details.")
  val rs7: Option[String] = None,
  @arg(doc="Read structure for the index read of the Illumina i5 sequencing primer.  The total bases must match the length of the index read. See the DemuxFastqs help message for more details.")
  val rs5: Option[String] = None,
  @arg(flag="ss", doc="The Sample Sheet (SampleSheet.csv).")
  val sampleSheet: FilePath,
  @arg(flag="o", doc="The output directory, with a BAM per sample.")
  val output: DirPath,
  @arg(flag="m", doc="Per-barcode and per-lane metrics written to this file.")
  val metrics: FilePath,
  @arg(flag="u", doc="Output BAM file name for the unmatched records.")
  val unmatched: String = DemuxFastqs.UnmatchedSampleId + ".bam",
  @arg(flag="l", doc="The lane number to analyze.") val lane: Option[Int] = None,
  @arg(flag="t", doc="The SAM tag for the molecular barcode.  If multiple molecular barcodes are given, they will be concatenated and stored here.")
  val umiTag: String = "RX",
  @arg(flag="mt", doc="A SAM tag per multiple molecular barcode if multiple barcodes are present. Either one tag can be given for all molecular barcodes or the same number as individual molecular barcode segments.", minElements = 0)
  val multiUmiTags: Seq[String] = Seq.empty,
  @arg(flag="q", doc="A value describing how the quality values are encoded in the FASTQ.  Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.")
  var qualityFormat: Option[FastqQualityFormat] = None,
  @arg(flag="pl", doc="The platform type (e.g. illumina, solid) to insert into the read group header") val platform: Option[String] = Some("ILLUMINA"),
  @arg(flag="cn", doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
  @arg(flag="pi", doc="Predicted median insert size, to insert into the read group header") val predictedInsertSize: Option[Integer] = None,
  @arg(flag="pg", doc="Program group to insert into the read group header.") val programGroup: Option[String] = None,
  @arg(flag="pm", doc="Platform model to insert into the group header (free-form text providing further details of the platform/technology used)") val platformModel: Option[String] = None,
  @arg(flag="co", doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
  @arg(flag="ds", doc="Inserted into the read group header") val description: Option[String] = None,
  @arg(flag="dt", doc="Date the run was produced, to insert into the read group header") val runDate: Option[Iso8601Date] = None,
  @arg(flag="so", doc="The sort order for the output sam/bam file.") val sortOrder: SortOrder = SortOrder.queryname,
  @arg(doc="Minimum base quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.") val minQ: Int = 0,
  @arg(doc="Maximum base quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.") val maxQ: Int = SAMUtils.MAX_PHRED_SCORE,
  @arg(doc="Allow (and ignore) empty lines") val allowAndIgnoreEmptyLines: Boolean = false,
  @arg(doc="Maximum mismatches for a barcode to be considered a match.") val maxMismatches: Int = 1,
  @arg(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.") val minMismatchDelta: Int = 2,
  @arg(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.") val maxNoCalls: Int = 2
) extends FgBioTool with LazyLogging {

  import DemuxFastqs._

  // Check all input and output files are readable or writable respectively.
  Io.assertReadable(fastq1 ++ fastq2 ++ fastq7 ++ fastq5)
  Io.assertReadable(this.sampleSheet)
  Io.assertWritableDirectory(this.output)
  Io.assertCanWriteFile(this.metrics)

  override def execute(): Unit = {
    // Make some data structures to help us on our way
    val r1ReadStructure          = ReadStructure(rs1)
    val r2ReadStructure          = rs2.map(ReadStructure(_))
    val i7ReadStructure          = rs7.map(ReadStructure(_))
    val i5ReadStructure          = rs5.map(ReadStructure(_))

    // Verify that the template bases are only on read one or read two
    if (r1ReadStructure.template.isEmpty)           throw new IllegalArgumentException("No template bases found in the read structure for read one.")
    if (r2ReadStructure.exists(_.template.isEmpty)) throw new IllegalArgumentException("No template bases found in the read structure for read two.")
    if (i7ReadStructure.exists(_.template.nonEmpty)) throw new IllegalArgumentException("Template bases not allowed in the i7 read.")
    if (i5ReadStructure.exists(_.template.nonEmpty)) throw new IllegalArgumentException("Template bases not allowed in the i5 read.")

    val fastqs = Seq(fastq1, fastq2, fastq7, fastq5)
    val readStructures: Seq[ReadStructure] = Seq(Some(r1ReadStructure), r2ReadStructure, i7ReadStructure, i5ReadStructure).flatten
    if (readStructures.length != fastqs.count(_.nonEmpty)) throw new IllegalArgumentException("The # of read structures did not match the number of fastqs.")

    // check the molecular barcodes and associated tags, if any exist.
    val molecularBarcodes = readStructures.flatMap { rs => rs.molecularBarcode }
    if (molecularBarcodes.isEmpty && multiUmiTags.nonEmpty) {
      throw new IllegalArgumentException("No molecular barcodes found but molecular barcode SAM tags were given.")
    }
    if (multiUmiTags.length > 1 && multiUmiTags.length != molecularBarcodes.length) {
      throw new IllegalArgumentException("Either one SAM tag or the same number as the molecular barcode segments should be given " +
        s"(found ${molecularBarcodes.length} molecular barcodes but given ${multiUmiTags.length} tags."
      )
    }

    // Read the sample sheet, convert to a custom sample type, and set the barcode bases for each sample
    val sampleSheet: SampleSheet = new SampleSheet(
      SampleSheet(this.sampleSheet, lane=lane)
        .map(new CustomSample(_))
        .toSeq
    ).sampleBarcodes(readStructures:_*)

    // Determine the quality format
    determineQualityFormat()

    // Create the files for writing
    val writers = createSamFileWriters(sampleSheet)

    // Open the input FASTQs for reading
    val nonEmptyFastqs = fastqs.filter(_.nonEmpty)
    val reader = new ParallelFastqReader(nonEmptyFastqs.head, nonEmptyFastqs.tail:_*)

    // Read in the quadruples of FASTQ records
    val fastqConverter = new FastqConverter(
      reader,
      writers.map(_.getFileHeader),
      sampleSheet,
      umiTag,
      multiUmiTags,
      qualityFormat.getOrElse(unreachable("Quality format not set")),
      minQ,
      maxQ,
      maxMismatches,
      minMismatchDelta,
      maxNoCalls,
      readStructures:_*)

    val progress = new com.fulcrumgenomics.util.ProgressLogger(logger)
    fastqConverter.foreach { assignedRec =>
      val writer = writers(assignedRec.sampleOrdinal-1)
      assignedRec.records.foreach { rec => writer.addAlignment(rec) }
      progress.record(assignedRec.records:_*)
    }

    // close them all
    reader.close()
    writers.foreach(_.close())

    // write the metrics
    Metric.write[SampleBarcodeMetric](metrics=fastqConverter.getBarcodeMetrics.toSeq, path=metrics)
  }

  private def determineQualityFormat(): Unit = {
    val readers = fastq1.map { fastq => new FastqReader(fastq.toFile, allowAndIgnoreEmptyLines) }
    val detector: QualityEncodingDetector = new QualityEncodingDetector
    detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers:_*)
    readers.foreach(_.close())
    val format = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, qualityFormat.orNull)
    if (detector.isDeterminationAmbiguous) {
      logger.warning("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.")
    }
    logger.info(String.format("Auto-detected quality format as: %s.", format))
    this.qualityFormat = Some(format)
  }

  /** Creates SAM file writers for each sample in the Sample Sheet, returning them in the same order as specified
    * by each sample's ordinal.  An extra writer is appended for records that do not match a sample barcode. */
  private def createSamFileWriters(sampleSheet: SampleSheet): List[SAMFileWriter] = {
    val writers: List[SAMFileWriter] = sampleSheet.map { sample =>
      val outputBam: Path       = DemuxFastqs.sampleOutputBam(this.output, sample)
      val header: SAMFileHeader = createSamFileHeader(sample)
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputBam.toFile)
      writer
    }.toList
    val unmatchedWriter = {
      val outputFile: File      = new File(this.output.toFile, this.unmatched)
      val header: SAMFileHeader = createSamFileHeader(new Sample(writers.length - 1, DemuxFastqs.UnmatchedSampleId, DemuxFastqs.UnmatchedSampleId, Some(DemuxFastqs.UnmatchedSampleId), None, None))
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile)
       writer
    }
    writers ++ List(unmatchedWriter)
  }

  /** Creates a read group for the SAM header with the values provided on the command line.. */
  private def createSamReadGroupRecord(readGroupId: String, sampleName: String, libraryName: String, platformUnit: String): SAMReadGroupRecord = {
    val rgroup: SAMReadGroupRecord = new SAMReadGroupRecord(readGroupId)
    rgroup.setSample(sampleName)
    rgroup.setLibrary(libraryName)
    rgroup.setPlatformUnit(platformUnit)
    platform.foreach(rgroup.setPlatform)
    sequencingCenter.foreach(rgroup.setSequencingCenter)
    predictedInsertSize.foreach(rgroup.setPredictedMedianInsertSize)
    description.foreach(rgroup.setDescription)
    runDate.foreach(rgroup.setRunDate)
    platformModel.foreach(rgroup.setPlatformModel)
    programGroup.foreach(rgroup.setProgramGroup)
    rgroup
  }

  /** Creates a simple header with the values provided on the command line. */
  private def createSamFileHeader(sample: Sample): SAMFileHeader = {
    val header: SAMFileHeader = new SAMFileHeader
    val libraryId: String     = sample.libraryId.getOrElse(sample.sampleName)
    val platformUnit: String  = sample.sampleBarcode.map(sb => s"${sample.sampleName}.$sb").getOrElse(sample.sampleName)
    header.addReadGroup(createSamReadGroupRecord(sample.sampleId, sample.sampleName, libraryId, platformUnit))
    this.comments.foreach(header.addComment)
    header.setSortOrder(this.sortOrder)
    header
  }
}