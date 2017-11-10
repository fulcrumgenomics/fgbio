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
import com.fulcrumgenomics.bam.api.QueryType.QueryType
import htsjdk.samtools._
import htsjdk.samtools.util.{Interval, Locatable}

import scala.collection.IterableView

/** Companion to the [[SamSource]] class that provides factory methods for sources. */
object SamSource {
  var DefaultUseAsyncIo: Boolean = false
  var DefaultValidationStringency: ValidationStringency = ValidationStringency.STRICT

  /**
    * Constructs a [[SamSource]] to read from the provided path.
    *
    * @param path the path to read the SAM/BAM/CRAM from
    * @param index an optional path to read the index from
    * @param ref an optional reference sequencing for decoding CRAM files
    * @param async if true use extra thread(s) to speed up reading
    * @param stringency the validation stringency to apply when reading the data
    * @param factory a SAMRecordFactory; MUST return classes that mix in [[SamRecord]]
    */
  def apply(path: PathToBam,
            index: Option[FilePath] = None,
            ref: Option[PathToFasta] = None,
            async: Boolean = DefaultUseAsyncIo,
            stringency: ValidationStringency = DefaultValidationStringency,
            factory: SAMRecordFactory = SamRecord.Factory): SamSource = {
    // Configure the factory
    val fac = SamReaderFactory.make()
    fac.samRecordFactory(factory)
    fac.setUseAsyncIo(async)
    fac.validationStringency(stringency)
    ref.foreach(r => fac.referenceSequence(r.toFile))

    // Open the input(s)
    val input = SamInputResource.of(path)
    index.foreach(i => input.index(i))
    new SamSource(fac.open(input))
  }
}

/** Describes the two types of queries that can be performed. */
object QueryType extends Enumeration {
  val Overlapping, Contained = Value
  type QueryType = Value
}

/**
  * A source class for reading SAM/BAM/CRAM files and for querying them.
  * @param reader the underlying [[SamReader]]
  */
class SamSource private(private val reader: SamReader) extends IterableView[SamRecord, SamSource] with HeaderHelper with Closeable {
  /** The [[htsjdk.samtools.SAMFileHeader]] associated with the source. */
  override val header: SAMFileHeader = reader.getFileHeader

  /** Required method for [[IterableView]]. */
  override protected def underlying: SamSource = this

  /** True if an index exists and query() calls can be made, false otherwise. */
  def indexed: Boolean = reader.hasIndex

  /** Returns an iterator over all the records in the source. */
  override def iterator: SamIterator = new SamIterator(reader.getFileHeader, reader.iterator())

  /** Returns an iterator over the records in the regions provided. */
  def query(regions: TraversableOnce[Locatable], queryType: QueryType = QueryType.Overlapping): SamIterator = {
    val queries = regions.map(l => new QueryInterval(dict.getSequenceIndex(l.getContig), l.getStart, l.getEnd)).toArray
    val contained = queryType == QueryType.Contained
    new SamIterator(header, reader.query(queries, contained))
  }

  /** Returns an iterator over the records in the region provided. */
  def query(chrom: String, start: Int, end: Int, queryType: QueryType): SamIterator = {
    query(List(new Interval(chrom, start, end)), queryType)
  }

  /** Returns an iterator over all the unmapped reads, without positions, at the end of the source. */
  def unmapped: SamIterator = new SamIterator(header, reader.queryUnmapped())

  /** Provides a string that shows where the source is reading from. */
  override def toString: String = s"SamReader(${reader.getResourceDescription})"

  override def close(): Unit = this.reader.close()

  /**
    * Returns the underlying SamReader. This should be avoided as much as possible, and the
    * SamSource should not be used again after calling [[toSamReader]].
    */
  def toSamReader: SamReader = reader
}
