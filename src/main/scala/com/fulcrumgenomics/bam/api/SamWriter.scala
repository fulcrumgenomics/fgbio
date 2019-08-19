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

package com.fulcrumgenomics.bam.api

import java.io.Closeable

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.util.{ProgressLogger, Sorter}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._

object SamWriter extends LazyLogging {
  var DefaultCompressionLevel: Int = 5
  var DefaultUseAsyncIo: Boolean      = true
  var DefaultMaxRecordsInRam: Int  = 1e6.toInt
  var DefaultCreateIndex: Boolean  = true

  /**
    * Constructs a SamWriter for writing out [[SamRecord]]s.
    *
    * @param path the path at which the data should be written
    * @param header the header for the SAM/BAM
    * @param sort an optional SamOrder into which to sort the reads. If no order is provided the reads are
    *             assumed to be in the order described in the header. If an order is provided the reads are
    *             _always_ sorted using that SamOrder before being emitted.
    * @param ref an optional reference sequence for use in writing CRAM
    * @param async if true use multiple threads to increase writing throughput
    * @param buffer the buffer size, in bytes, to use for the writer
    * @param compression the GZIP compression level to use when compressing data
    * @param tmp an optional tmp directory to use when sorting
    * @param maxRecordsInRam the maximum number of records to keep in RAM while sorting
    * @param index if true, index the SAM/BAM file on the fly while writing
    * @param md5 if true, create an MD5 file for the output file while writing
    * @param sortProgress an optional ProgressLogger to use if/when putting records into a Sorter
    * @param writeProgress an optional ProgressLogger to use when writing records out
    */
  def apply(path: PathToBam,
            header: SAMFileHeader,
            sort: Option[SamOrder]   = None,
            ref: Option[PathToFasta] = None,
            async: Boolean           = DefaultUseAsyncIo,
            buffer: Int              = Defaults.NON_ZERO_BUFFER_SIZE,
            compression: Int         = DefaultCompressionLevel,
            tmp: Option[DirPath]     = None,
            maxRecordsInRam: Int     = 1e6.toInt,
            index: Boolean           = true,
            md5: Boolean             = false,
            sortProgress:  Option[ProgressLogger] = Some(ProgressLogger(logger, noun="records", verb="Sorted", unit=2e6.toInt)),
            writeProgress: Option[ProgressLogger] = Some(ProgressLogger(logger, noun="records", verb="Wrote ", unit=2e6.toInt))
           ): SamWriter = {

    val sorter = sort match {
      case None => None
      case Some(so) if so == SamOrder.Unsorted || so == SamOrder.Unknown =>
        so.applyTo(header)
        None
      case Some(so) =>
        so.applyTo(header)
        Some(new Sorter(maxRecordsInRam, new SamRecordCodec(header), so.sortkey))
    }

    // Only index the output file if the sort order is right, indexing is requested, and no chrom is too long
    val actuallyIndex = if (header.getSortOrder != SortOrder.coordinate || !index) false else {
      val hasChromsTooLong = header.getSequenceDictionary.getSequences.exists(_.getSequenceLength > GenomicIndexUtil.BIN_GENOMIC_SPAN)
      if (hasChromsTooLong) logger.warning(s"Cannot index $path as one or more chromosomes is too long.")
      !hasChromsTooLong
    }

    val factory = new SAMFileWriterFactory()
    factory.setUseAsyncIo(async)
    factory.setBufferSize(buffer)
    factory.setCompressionLevel(compression)
    factory.setCreateIndex(actuallyIndex)
    factory.setCreateMd5File(md5)
    tmp.foreach(dir => factory.setTempDirectory(dir.toFile))

    new SamWriter(writer=factory.makeWriter(header, true, path.toFile, ref.map(_.toFile).orNull),
      sorter=sorter, sortProgress=sortProgress, writeProgress=writeProgress)
  }
}

/** Provides the ability to write [[SamRecord]]s to an output Path. */
final class SamWriter private (private val writer: SAMFileWriter,
                               private val sorter: Option[Sorter[SamRecord,_]],
                               private val sortProgress: Option[ProgressLogger],
                               private val writeProgress: Option[ProgressLogger]
                              ) extends Closeable with Writer[SamRecord] with HeaderHelper {

  /** The associated [[htsjdk.samtools.SAMFileHeader]]. */
  override def header: SAMFileHeader = writer.getFileHeader

  /** Writes an individual item. */
  override def write(rec: SamRecord): Unit = sorter match {
    case Some(s) =>
      s += rec
      sortProgress.foreach(_.record(rec))
    case None    =>
      this.writer.addAlignment(rec.asSam)
      writeProgress.foreach(_.record(rec))
  }

  /** Closes the writer. */
  override def close(): Unit = {
    this.sorter.foreach { s =>
      s.iterator.foreach { rec =>
        this.writer.addAlignment(rec.asSam)
        writeProgress.foreach(_.record(rec))
      }

      s.close()
    }

    this.writer.close()
  }
}
