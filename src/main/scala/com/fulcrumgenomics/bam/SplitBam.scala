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

import java.io.Closeable
import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{PathPrefix, PathToBam, SafelyClosable, javaIterableToIterator, javaIteratorAsScalaIterator}
import com.fulcrumgenomics.commons.io.{Io, PathUtil}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import enumeratum.EnumEntry
import htsjdk.samtools._

import scala.collection.immutable.IndexedSeq

sealed trait SplitType extends EnumEntry
object SplitType extends FgBioEnum[SplitType] {
  def values: IndexedSeq[SplitType] = findValues
  case object Library extends SplitType
  case object ReadGroup extends SplitType
}

@clp(group = ClpGroups.SamOrBam, description=
  """
    |Splits a BAM into multiple BAMs, one per-read group (or library).
    |
    |The resulting BAMs will be named `<output-prefix>.<read-group-id>.bam`, or `<output-prefix>.<library-name>.bam`
    |when splitting by the library.  All reads without a read group, or without a library when splitting by library,
    |will be written to `<output-prefix>.unknown.bam`.  If no such reads exist, then no such file will exist.
  """)
class SplitBam
( @arg(flag='i', doc="Input SAM or BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output prefix for all SAM or BAM files (ex. output/sample-name).") val output: PathPrefix,
  @arg(flag='s', doc="Split by library instead of read group") val splitBy: SplitType = SplitType.ReadGroup,
  @arg(flag='u', doc="The name to use for the unknown file") val unknown: String = "unknown"
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in       = SamSource(input)
    val progress = ProgressLogger(logger)
    val unknownBamAndWriter = {
      val unknownBam: PathToBam = toOutput(unknown)
      val unknownWriter: SamWriter = SamWriter(toOutput(unknown), in.header)
      //new SAMFileWriterFactory().makeWriter(in.h, true, toOutput(unknown).toFile, null)
      WriterInfo(name=unknown, bam=unknownBam, writer=unknownWriter)
    }
    val writers  = createWriters(header=in.header, splitBy=splitBy).withDefaultValue(unknownBamAndWriter)
    val counter  = new SimpleCounter[WriterInfo]()

    in.foreach { rec =>
      val info = writers(rec.readGroup)
      info.writer += rec
      counter.count(info)
      progress.record(rec)
    }

    writers.values.foreach(_.close())
    unknownBamAndWriter.close()
    in.safelyClose()

    counter.toSeq sortBy { _._1.name } foreach { case (info, count) =>
      logger.info(s"Wrote $count records to ${info.bam.getFileName}")
    }

    if (counter.countOf(unknownBamAndWriter) == 0) {
      Files.delete(unknownBamAndWriter.bam)
    }
  }

  /** Stores the path to the BAM and the associated writer */
  private case class WriterInfo(name: String, bam: PathToBam, writer: SamWriter) extends Closeable{
    def close(): Unit = writer.close()
  }

  /** Initializes the writers and returns a map from each read group to a writer. */
  private def createWriters(header: SAMFileHeader, splitBy: SplitType): Map[SAMReadGroupRecord, WriterInfo] = {
    splitBy match {
      case SplitType.Library =>
        header.getReadGroups.toSeq.groupBy { rg => rg.getLibrary }
          .flatMap { case (library, readGroups) =>
            val bam = toOutput(library)
            val writer = SamWriter(bam, header)
            readGroups.map { rg => rg -> WriterInfo(name=library, bam=bam, writer=writer) }
          }
      case SplitType.ReadGroup =>
        header.getReadGroups.map { rg =>
          val bam = toOutput(rg.getId)
          val writer = SamWriter(bam, header)
          rg -> WriterInfo(name=rg.getId, bam=bam, writer=writer)
        }.toMap
    }
  }

  private[bam] def toOutput(name: String): PathToBam = {
    val outputDir = output.getParent
    val prefix    = output.getFileName
    outputDir.resolve(PathUtil.sanitizeFileName(s"$prefix.$name.bam"))
  }
}

