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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.PathToBam
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fastq.{FastqRecord, FastqWriter}
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.SequenceUtil

object SamToFastq {
  /** The tag that stores bases to prepend to SEQ. */
  val FivePrimeBasesTag: String      = "S5"
  /** The tag that stores qualities to prepend to SEQ. */
  val FivePrimeQualitiesTag: String  = "Q5"
  /** The tag that stores bases to append to SEQ. */
  val ThreePrimeBasesTag: String     = "S3"
  /** The tag that stores qualities to append to SEQ. */
  val ThreePrimeQualitiesTag: String = "Q3"
}

@clp(description =
  """
    |Converts a SAM or BAM to FASTQ.
    |
    |Three outputs will be created:
    |1. Fragment reads will be written to the file `<output>.fragments.fq.gz`
    |2. First of pair reads will be written to the file `<output>.R1.fq.gz`
    |3. Second of pair reads will be written to the file `<output>.R2.fq.gz`
    |
    |Bases and qualities stored in the five prime and three prime bases/qualities tags will be prepended and appended
    |to the bases/qualities found in the SEQ/QUAL field respectively.  The SEQ field will be reverse complemented if
    |mapped to the reverse strand.
  """,
  group = ClpGroups.SamOrBam)
class SamToFastq
( @arg(flag='i', doc="Input BAM file.")                               val input: PathToBam,
  @arg(flag='o', doc="Output FASTQ file prefix.")                     val output: FilePath,
  @arg(doc="The SAM tag storing the 5' bases to prepend to SEQ.")     val fivePrimeBasesTag: String = SamToFastq.FivePrimeBasesTag,
  @arg(doc="The SAM tag storing the 5' qualities to prepend to SEQ.") val fivePrimeQualitiesTag: String = SamToFastq.FivePrimeQualitiesTag,
  @arg(doc="The SAM tag storing the 3' bases to append to SEQ.")     val threePrimeBasesTag: String = SamToFastq.ThreePrimeBasesTag,
  @arg(doc="The SAM tag storing the 3' qualities to append to SEQ.") val threePrimeQualitiesTag: String = SamToFastq.ThreePrimeQualitiesTag
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  validate(fivePrimeBasesTag.length      == 2, s"--five-prime-bases-tag must be length two, was '${fivePrimeBasesTag.length}'")
  validate(fivePrimeQualitiesTag.length  == 2, s"--five-prime-qualities-tag must be length two, was '${fivePrimeQualitiesTag.length}'")
  validate(threePrimeBasesTag.length     == 2, s"--three-prime-bases-tag must be length two, was '${threePrimeBasesTag.length}'")
  validate(threePrimeQualitiesTag.length == 2, s"--three-prime-qualities-tag must be length two, was '${threePrimeQualitiesTag.length}'")

  override def execute(): Unit = {
    val writers: Map[ReadType, FastqWriter] = ReadType.values.map {
      case _: ReadType.Fragment.type => ReadType.Fragment -> FastqWriter(PathUtil.pathTo(output + ".fragments.fq.gz"))
      case _: ReadType.ReadOne.type  => ReadType.ReadOne  -> FastqWriter(PathUtil.pathTo(output + ".R1.fq.gz"))
      case _: ReadType.ReadTwo.type  => ReadType.ReadTwo  -> FastqWriter(PathUtil.pathTo(output + ".R2.fq.gz"))
    }.toMap
    val iterator = buildIterator()
    val progress = ProgressLogger(logger, verb="written", unit=5e6.toInt)

    iterator.foreach { rec =>
      val readType = ReadType.from(rec)
      val writer   = writers(readType)
      val fastq    = toFastqRecord(rec, readType)
      writer += fastq
      progress.record()
    }

    writers.values.foreach(_.safelyClose)
  }

  /** Builds an iterator ensuring the records are returned in queryname order.  If the input is not queryname sorted,
    * the records are sorted. */
  private def buildIterator(): Iterator[SamRecord] = {
    val in       = SamSource(input)
    val order = SamOrder(in.header).getOrElse(SamOrder.Unsorted)

    if (order == SamOrder.Queryname) {
      new SelfClosingIterator(in.iterator.bufferBetter, () => in.safelyClose())
    }
    else {
      logger.info("Sorting into queryname order.")
      val progress = ProgressLogger(this.logger, "Queryname sorted")
      val sort     = Bams.sorter(SamOrder.Queryname, in.header)
      in.foreach { rec =>
        sort += rec
        progress.record(rec)
      }

      new SelfClosingIterator(sort.iterator.bufferBetter, () => sort.close())
    }
  }

  /** Converts a [[SamRecord]] into a [[FastqRecord]].  Prepends and appends to the bases the sequence stored in the
    * five prime and three prime bases tags respectively.  Similarly, prepends and appends to the qualities the qualities
    * stored in the five prime and three prime qualities tags respectively.  */
  private def toFastqRecord(rec: SamRecord, readType: ReadType): FastqRecord = {
    val bases = rec.bases
    val quals = if (rec.positiveStrand) rec.qualsString else rec.qualsString.reverse

    if (rec.negativeStrand) SequenceUtil.reverseComplement(bases)

    val fivePrimeBases  = rec.get(this.fivePrimeBasesTag).getOrElse("")
    val fivePrimeQuals  = rec.get(this.fivePrimeQualitiesTag).getOrElse("")
    val threePrimeBases = rec.get(this.threePrimeBasesTag).getOrElse("")
    val threePrimeQuals = rec.get(this.threePrimeQualitiesTag).getOrElse("")

    val readNumber = readType match {
      case _: ReadType.Fragment.type => None
      case _: ReadType.ReadOne.type  => Some(1)
      case _: ReadType.ReadTwo.type  => Some(2)
    }

    FastqRecord(
      name       = rec.name,
      bases      = fivePrimeBases + new String(bases) + threePrimeBases,
      quals      = fivePrimeQuals + quals + threePrimeQuals,
      readNumber = readNumber
    )
  }
}

