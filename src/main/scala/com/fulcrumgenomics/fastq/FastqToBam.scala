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

package com.fulcrumgenomics.fastq

import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.{Io, ProgressLogger, ReadStructure}
import dagr.commons.CommonsDef.PathToFastq
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.Iso8601Date
import htsjdk.samtools.{ReservedTagConstants, SAMFileHeader, SAMFileWriter, SAMFileWriterFactory, SAMReadGroupRecord, SAMRecord}
@clp(group=ClpGroups.Fastq, description=
  """
    |DOES STUFF!!
  """)
class FastqToBam
(
  @arg(flag="i", doc="Fastq files corresponding to each sequencing read (e.g. R1, I1, etc.).") val input: Seq[PathToFastq],
  @arg(flag="o", doc="The output SAM or BAM file to be written.")                              val output: PathToBam,
  @arg(flag="r", doc="Read structures, one for each of the FASTQs.")                           val readStructures: Seq[ReadStructure],
  @arg(flag="s", doc="Sort order for the output sam/bam file (e.g. unsorted or queryname).")   val sortOrder: SortOrder = SortOrder.unsorted,
  @arg(flag="u", doc="Tag in which to store molecular barcodes/UMIs.")                         val umiTag: String = ConsensusTags.UmiBases,
  @arg(          doc="Read group ID to use in the file header.")                               val readGroupId: String = "A",
  @arg(          doc="Sequencing Platform.")                                                   val platform: String = "illumina",
  @arg(doc="Platform unit (e.g. '<flowcell-barcode>.<lane>.<sample-barcode>')")                val platformUnit: Option[String] = None,
  @arg(doc="The sequencing center from which the data originated")                             val sequencingCenter: Option[String] = None,
  @arg(doc="Predicted median insert size, to insert into the read group header")               val predictedInsertSize: Option[Integer] = None,
  @arg(doc="Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)")    val platformModel: Option[String] = None,
  @arg(doc="Comment(s) to include in the merged output file's header.", minElements = 0)       val comment: List[String] = Nil,
  @arg(doc="Date the run was produced, to insert into the read group header")                  val runDate: Option[Iso8601Date] = None
)
  extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  validate(input.length == readStructures.length, "input and read-structure must be supplied the same number of times.")
  validate(Range.inclusive(1,2).contains(readStructures.flatMap(_.template).size), "read structures must contain 1-2 template reads total.")

  override def execute(): Unit = {
    val encoding = qualityEncoding
    val writer   = makeSamWriter()
    val readers  = this.input.map(FastqSource(_))

    val progress = this.sortOrder match {
      case SortOrder.unsorted => ProgressLogger(logger, verb="written")
      case _ =>
        writer.setProgressLogger(ProgressLogger(logger, verb="writter"))
        ProgressLogger(logger, verb="read")
    }

    while (readers.forall(_.hasNext)) {
      val recs = makeSamRecords(readers.map(_.next()), readStructures, writer.getFileHeader, encoding)
      recs.foreach(writer.addAlignment)
      recs.foreach(progress.record)
    }

    writer.close()
    if (!readers.forall(r => !r.hasNext)) fail("Fastq files appear to have different numbers of records.")
  }

  /** Makes the SAMFileWriter we'll use to output the file. */
  protected def makeSamWriter(): SAMFileWriter = {
    val header = new SAMFileHeader
    header.setSortOrder(this.sortOrder)
    if (this.sortOrder == SortOrder.unsorted) header.setGroupOrder(GroupOrder.query)
    header.setComments(util.Arrays.asList(this.comment:_*))

    val rg = new SAMReadGroupRecord(this.readGroupId)
    rg.setPlatform(this.platform)
    this.platformUnit.foreach(pu => rg.setPlatformUnit(pu))
    this.sequencingCenter.foreach(cn => rg.setSequencingCenter(cn))
    this.predictedInsertSize.foreach(isize => rg.setPredictedMedianInsertSize(isize))
    this.platformModel.foreach(pm => rg.setPlatformModel(pm))
    this.runDate.foreach(date => rg.setRunDate(date))
    header.addReadGroup(rg)

    val factory =  new SAMFileWriterFactory().setCreateIndex(false)
    factory.makeWriter(header, this.sortOrder == SortOrder.unsorted, this.output.toFile, null)
  }

  protected def makeSamRecords(fqs: Seq[FastqRecord],
                               rss: Seq[ReadStructure],
                               header: SAMFileHeader,
                               encoding: QualityEncoding
                              ): Seq[SAMRecord] = {
    val subs = fqs.iterator.zip(rss.iterator).flatMap { case(fq, rs) => rs.structureReadWithQualities(fq.bases, fq.quals, strict=false)}.toIndexedSeq
    val sampleBarcode = subs.iterator.filter(_.segment.symbol == 'B').map(_.bases).mkString("-")
    val umi           = subs.iterator.filter(_.segment.symbol == 'M').map(_.bases).mkString("-")
    val templates     = subs.iterator.filter(_.segment.symbol == 'T').toList

    templates.zipWithIndex.map { case (read, index) =>
      val rec = new SAMRecord(header)
      rec.setAttribute(ReservedTagConstants.READ_GROUP_ID, this.readGroupId)
      rec.setReadName(fqs.head.name)
      rec.setReadString(read.bases)
      rec.setBaseQualities(encoding.toStandardNumeric(read.quals.getOrElse(unreachable("Quals must be here!"))))
      rec.setReadUnmappedFlag(true)
      if (templates.size == 2) {
        rec.setReadPairedFlag(true)
        rec.setMateUnmappedFlag(true)
        if (index == 1) rec.setFirstOfPairFlag(true) else rec.setSecondOfPairFlag(true)
      }

      if (sampleBarcode.nonEmpty) rec.setAttribute("BC", sampleBarcode)
      if (umi.nonEmpty) rec.setAttribute(this.umiTag, umi)

      rec
    }
  }


  /** Determine the quality encoding of the incoming fastq files. */
  protected def qualityEncoding: QualityEncoding = {
    val readers  = input.map { fastq => FastqSource(fastq) }
    val iterator = readers.tail.foldLeft(readers.head.iterator) {(a,b) => a ++ b }.map(_.quals)
    val detector = new QualityEncodingDetector
    detector.sample(iterator)
    readers.foreach(_.safelyClose())
    detector.rankedCompatibleEncodings(q=30) match {
      case Nil        => fail("Quality scores in FASTQ file do not match any known encoding.")
      case enc :: Nil => yieldAndThen(enc) { logger.info("Detected fastq quality encoding: ", enc) }
      case enc :: xs  =>
        logger.info(s"Could not uniquely determine quality encoding. Using $enc, other valid encodings: ${xs.mkString(", ")}")
        enc
    }
  }
}
