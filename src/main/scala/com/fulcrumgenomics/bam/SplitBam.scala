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
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.{IteratorToJavaCollectionsAdapter, PathPrefix, PathToBam, SafelyClosable, javaIterableToIterator}
import com.fulcrumgenomics.commons.io.{Io, PathUtil, Writer}
import com.fulcrumgenomics.commons.util.LazyLogging
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
    |
    |For splitting a BAM with many read groups, use `--no-async-writing` in case too many threads or memory are used.
  """)
class SplitBam
( @arg(flag='i', doc="Input SAM or BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output prefix for all SAM or BAM files (ex. output/sample-name).") val output: PathPrefix,
  @arg(flag='s', doc="Split by library instead of read group") val splitBy: SplitType = SplitType.ReadGroup,
  @arg(flag='u', doc="The name to use for the unknown file") val unknown: String = "unknown",
  @arg(doc="Do not write the records asynchronously. Use this if with a large numbers of read groups to reduce memory usage.")
  val noAsyncWriting: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)

  private val outputIsDirectory: Boolean = Files.isDirectory(output)
  if (outputIsDirectory) Io.assertWritableDirectory(output)
  else Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in       = SamSource(input)
    val progress = ProgressLogger(logger=logger)
    val unknownBamAndWriter: WriterInfo = {
      val unknownBam: PathToBam = toOutput(name=unknown)
      val unknownWriter: SamWriter = SamWriter(toOutput(unknown), toOutputHeader(in.header))
      WriterInfo(name=unknown, bam=unknownBam, writer=unknownWriter)
    }
    val writerInfoMap: Map[SAMReadGroupRecord, WriterInfo] = {
      createWriters(header=in.header, splitBy=splitBy).withDefaultValue(unknownBamAndWriter)
    }

    in.foreach { rec =>
      val info = writerInfoMap(rec.readGroup)
      info.count += 1
      info.writer.write(rec)
      progress.record(rec)
    }
    progress.logLast()

    in.safelyClose()

    writerInfoMap.values.toSeq.sortBy(_.name).distinct.foreach { info =>
      info.writer.close()
      logger.info(f"Wrote ${info.count}%,d records to ${info.bam.getFileName}")
    }

    if (unknownBamAndWriter.count == 0) {
      Files.delete(unknownBamAndWriter.bam)
    }
  }

  /** Gets the output path for the writer with a given name. */
  private[bam] def toOutput(name: String): PathToBam = {
    if (outputIsDirectory) {
      output.resolve(f"$name.bam")
    }
    else {
      val outputDir = output.getParent
      val prefix    = output.getFileName
      outputDir.resolve(PathUtil.sanitizeFileName(s"$prefix.$name.bam"))
    }
  }

  /** Initializes the writers and returns a map from each read group to a writer. */
  private def createWriters(header: SAMFileHeader, splitBy: SplitType): Map[SAMReadGroupRecord, WriterInfo] = {
    val groups: Map[String, Seq[SAMReadGroupRecord]]  = splitBy match {
      case SplitType.Library   => header.getReadGroups.toSeq.groupBy { rg => rg.getLibrary }
      case SplitType.ReadGroup => header.getReadGroups.map(rg => rg.getId -> Seq(rg)).toMap
    }
    logger.info(f"Outputting to ${groups.size}%,d files ($splitBy).")
    groups.flatMap { case (library, readGroups) =>
      val bam    = toOutput(library)
      val writer = SamWriter(
        path   = bam,
        header = toOutputHeader(header=header, readGroups:_*),
        async  = !noAsyncWriting
      )
      readGroups.map { rg => rg -> WriterInfo(name=library, bam=bam, writer=writer) }
    }
  }

  /** Stores the path to the BAM and the associated writer */
  private case class WriterInfo(name: String, bam: PathToBam, writer: Writer[SamRecord], var count: Long = 0) extends Closeable {
    def close(): Unit = writer.close()
  }

  /** Efficiently creates an output header with a single read group from the input header.
    *
    * Developer note: `SamFileHeader.clone` is very memory inefficient.  It creates a `String` from the header to clone,
    * the decodes it!
    * */
  private def toOutputHeader(header: SAMFileHeader, readGroup: SAMReadGroupRecord*): SAMFileHeader = {
    val outputHeader = new SAMFileHeader()
    // Sort order and group order
    outputHeader.setSortOrder(header.getSortOrder)
    outputHeader.setGroupOrder(header.getGroupOrder)
    // HD
    outputHeader.setAttribute(SAMFileHeader.VERSION_TAG, header.getVersion)
    // SD
    outputHeader.setSequenceDictionary(header.getSequenceDictionary)
    // RG
    outputHeader.setReadGroups(readGroup.iterator.toJavaList)
    // PG
    outputHeader.setProgramRecords(header.getProgramRecords)
    // Comment
    outputHeader.setComments(header.getComments)
    outputHeader
  }
}

