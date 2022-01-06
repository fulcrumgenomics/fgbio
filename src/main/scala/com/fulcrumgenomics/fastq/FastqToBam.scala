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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.PathToFastq
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.SegmentType._
import com.fulcrumgenomics.util.{Io, ProgressLogger, ReadStructure}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.Iso8601Date
import htsjdk.samtools.{ReservedTagConstants, SAMFileHeader, SAMReadGroupRecord}

@clp(group=ClpGroups.Fastq, description=
  """
    |Generates an unmapped BAM (or SAM or CRAM) file from fastq files.  Takes in one or more fastq files (optionally
    |gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read
    |structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, or to
    |designate bases to be skipped over.
    |
    |Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files. Four kinds of
    |operators are recognized:
    |
    |1. `T` identifies a template read
    |2. `B` identifies a sample barcode read
    |3. `M` identifies a unique molecular index read
    |4. `S` identifies a set of bases that should be skipped or ignored
    |
    |The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote "all remaining
    |bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length.  For example
    |to convert a paired-end run with an index read and where the first 5 bases of R1 are a UMI and the second
    |five bases are monotemplate you might specify:
    |
    |```
    |--input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B
    |```
    |
    |Alternative if you know your reads are of fixed length you could specify:
    |
    |```
    |--input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B
    |```
    |
    |For more information on read structures see the
    |[Read Structure Wiki Page](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)
    |
    |The same number of input files and read structures must be provided, with one exception: if supplying exactly
    |1 or 2 fastq files, both of which are solely template reads, no read structures need be provided.
    |
    |The output file can be sorted by queryname using the `--sort-order` option; the default is to produce a BAM
    |with reads in the same order as they appear in the fastq file.
  """)
class FastqToBam
(
  @arg(flag='i', doc="Fastq files corresponding to each sequencing read (e.g. R1, I1, etc.).") val input: Seq[PathToFastq],
  @arg(flag='o', doc="The output SAM or BAM file to be written.")                              val output: PathToBam,
  @arg(flag='r', doc="Read structures, one for each of the FASTQs.", minElements=0)            val readStructures: Seq[ReadStructure] = Seq(),
  @arg(flag='s', doc="If true, queryname sort the BAM file, otherwise preserve input order.")  val sort: Boolean = false,
  @arg(flag='u', doc="Tag in which to store molecular barcodes/UMIs.")                         val umiTag: String = ConsensusTags.UmiBases,
  @arg(flag='q', doc="Tag in which to store molecular barcode/UMI qualities.")                 val umiQualTag: Option[String] = None,
  @arg(          doc="Read group ID to use in the file header.")                               val readGroupId: String = "A",
  @arg(          doc="The name of the sequenced sample.")                                      val sample: String,
  @arg(          doc="The name/ID of the sequenced library.")                                  val library: String,
  @arg(          doc="Sequencing Platform.")                                                   val platform: String = "illumina",
  @arg(doc="Platform unit (e.g. '<flowcell-barcode>.<lane>.<sample-barcode>')")                val platformUnit: Option[String] = None,
  @arg(doc="Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)")    val platformModel: Option[String] = None,
  @arg(doc="The sequencing center from which the data originated")                             val sequencingCenter: Option[String] = None,
  @arg(doc="Predicted median insert size, to insert into the read group header")               val predictedInsertSize: Option[Integer] = None,
  @arg(doc="Description of the read group.")                                                   val description: Option[String] = None,
  @arg(doc="Comment(s) to include in the output file's header.", minElements = 0)              val comment: List[String] = Nil,
  @arg(doc="Date the run was produced, to insert into the read group header")                  val runDate: Option[Iso8601Date] = None
)
  extends FgBioTool with LazyLogging {

  // If no read structures are provided and we only have 1-2 fastqs, assume that they are just template reads
  private val actualReadStructures = if (readStructures.isEmpty && (1 to 2 contains input.length)) input.map(_ => ReadStructure("+T")) else readStructures

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  validate(input.length == actualReadStructures.length, "input and read-structure must be supplied the same number of times.")
  validate(1 to 2 contains actualReadStructures.flatMap(_.templateSegments).size, "read structures must contain 1-2 template reads total.")

  override def execute(): Unit = {
    val encoding = qualityEncoding
    val writer   = makeSamWriter()
    val iterator = FastqSource.zipped(this.input.map(FastqSource(_)))

    iterator.foreach { fqs => writer ++= makeSamRecords(fqs, actualReadStructures, writer.header, encoding) }
    writer.close()
  }

  /** Makes the SAMFileWriter we'll use to output the file. */
  protected def makeSamWriter(): SamWriter = {
    val header = new SAMFileHeader
    header.setSortOrder(if (sort) SortOrder.queryname else SortOrder.unsorted)
    header.setGroupOrder(GroupOrder.query)
    header.setComments(util.Arrays.asList(this.comment:_*))

    val rg = new SAMReadGroupRecord(this.readGroupId)
    rg.setSample(sample)
    rg.setLibrary(library)
    rg.setPlatform(this.platform)
    this.platformUnit.foreach(pu => rg.setPlatformUnit(pu))
    this.sequencingCenter.foreach(cn => rg.setSequencingCenter(cn))
    this.predictedInsertSize.foreach(isize => rg.setPredictedMedianInsertSize(isize))
    this.platformModel.foreach(pm => rg.setPlatformModel(pm))
    this.runDate.foreach(date => rg.setRunDate(date))
    this.description.foreach(desc => rg.setDescription(desc))
    header.addReadGroup(rg)

    SamWriter(this.output, header, sort = if (sort) Some(SamOrder.Queryname) else None)
  }

  /** Generates SamRecords for each of the template reads across the read structures. */
  protected def makeSamRecords(fqs: Seq[FastqRecord],
                               structures: Seq[ReadStructure],
                               header: SAMFileHeader,
                               encoding: QualityEncoding
                              ): Seq[SamRecord] = {
    // Make the SamRecords inside a try so we can provide more informative error messages
    try {
      val subs = fqs.iterator.zip(structures.iterator).flatMap { case(fq, rs) => rs.extract(fq.bases, fq.quals) }.toIndexedSeq
      val sampleBarcode = subs.iterator.filter(_.kind == SampleBarcode).map(_.bases).mkString("-")
      val umi           = subs.iterator.filter(_.kind == MolecularBarcode).map(_.bases).mkString("-")
      val umiQual       = subs.iterator.filter(_.kind == MolecularBarcode).map(_.quals).mkString(" ")
      val templates     = subs.iterator.filter(_.kind == Template).toList

      templates.zipWithIndex.map { case (read, index) =>
        // If the template read had no bases, we'll substitute in a single N @ Q2 below to keep htsjdk happy
        val empty = read.bases.length == 0

        val rec = SamRecord(header)
        rec(ReservedTagConstants.READ_GROUP_ID) = this.readGroupId
        rec.name  = fqs.head.name
        rec.bases = if (empty) "N" else read.bases
        rec.quals = if (empty) Array[Byte](2) else encoding.toStandardNumeric(read.quals)
        rec.unmapped = true
        if (templates.size == 2) {
          rec.paired = true
          rec.mateUnmapped = true
          if (index == 0) rec.firstOfPair = true else rec.secondOfPair = true
        }

        if (sampleBarcode.nonEmpty) rec("BC") = sampleBarcode
        if (umi.nonEmpty) {
          rec(this.umiTag) = umi
          this.umiQualTag.foreach(rec(_) = umiQual)
        }

        rec
      }
    }
    catch {
      case ex: Exception => fail(s"Failed to process record(s) with name ${fqs.map(_.name).head} due to: ${ex.getMessage}")
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
