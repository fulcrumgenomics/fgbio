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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.QueryType.QueryType
import com.fulcrumgenomics.cmdline.FgBioCommonArgs
import htsjdk.samtools._
import htsjdk.samtools.util.{Interval, Locatable}

import java.io.{Closeable, InputStream}
import scala.collection.compat._

/** Companion to the [[SamSource]] class that provides factory methods for sources. */
object SamSource {
  var DefaultUseAsyncIo: Boolean = false
  var DefaultValidationStringency: ValidationStringency = ValidationStringency.STRICT

  /** Configure a SAM Record Factory with a variety of parameters. */
  private def buildSamRecordFactory(
    factory: SAMRecordFactory,
    ref: Option[PathToFasta],
    async: Boolean,
    stringency: ValidationStringency,
  ): SamReaderFactory = {
    val fac = SamReaderFactory.make()
    fac.samRecordFactory(factory)
    fac.setUseAsyncIo(async)
    fac.validationStringency(stringency)
    ref.foreach(fac.referenceSequence)
    fac
  }

  /** Constructs a [[SamSource]] to read from the provided path.
    *
    * @param path the path to read the SAM/BAM/CRAM from
    * @param index an optional path to read the index from
    * @param ref an optional reference sequence for decoding CRAM files
    * @param async if true use extra thread(s) to speed up reading
    * @param stringency the validation stringency to apply when reading the data
    * @param factory a [[SamRecord.Factory]]; MUST return classes that mix in [[SamRecord]]
    */
  def apply(path: PathToBam,
            index: Option[FilePath] = None,
            ref: Option[PathToFasta] = FgBioCommonArgs.args.cramRefFasta,
            async: Boolean = DefaultUseAsyncIo,
            stringency: ValidationStringency = DefaultValidationStringency,
            factory: SAMRecordFactory = SamRecord.Factory): SamSource = {
    val fac   = buildSamRecordFactory(factory = factory, ref = ref, async = async, stringency = stringency)
    val input = SamInputResource.of(path)
    index.foreach(input.index)
    new SamSource(fac.open(input))
  }

  /** Constructs a [[SamSource]] to read from the provided input stream.
    *
    * @param stream the input stream of SAM/BAM/CRAM bytes
    * @param ref an optional reference sequence for decoding CRAM files
    * @param async if true use extra thread(s) to speed up reading
    * @param stringency the validation stringency to apply when reading the data
    * @param factory a [[SamRecord.Factory]]; MUST return classes that mix in [[SamRecord]]
    */
  def apply(
    stream: InputStream,
    ref: Option[PathToFasta],
    async: Boolean,
    stringency: ValidationStringency,
    factory: SAMRecordFactory,
  ): SamSource = {
    val fac = buildSamRecordFactory(factory = factory, ref = ref, async = async, stringency = stringency)
    new SamSource(fac.open(SamInputResource.of(stream)), closer = Some(() => stream.close()))
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
class SamSource private(private val reader: SamReader, private val closer: Option[Closeable] = None)
  extends View[SamRecord] with HeaderHelper with Closeable {

  /** The [[htsjdk.samtools.SAMFileHeader]] associated with the source. */
  override val header: SAMFileHeader = reader.getFileHeader

  /** Required for 2.12 compatibility. */
  def underlying: Iterable[SamRecord] = this

  /** True if an index exists and query() calls can be made, false otherwise. */
  def indexed: Boolean = reader.hasIndex

  /** Returns an iterator over all the records in the source. */
  override def iterator: SamIterator = new SamIterator(reader.getFileHeader, reader.iterator())

  private def newQueryInterval(region: Locatable): QueryInterval = {
    val contig = dict.get(region.getContig).getOrElse(
      throw new NoSuchElementException(s"Contig '${region.getContig}' not in SAM/BAM header.")
    )
    val start  = Math.max(region.getStart, 1);
    val end    = Math.min(region.getEnd, contig.length)
    new QueryInterval(contig.index, start, end)
  }

  /** Returns an iterator over the records in the regions provided. */
  def query(regions: IterableOnce[Locatable], queryType: QueryType = QueryType.Overlapping): SamIterator = {
    val queries = QueryInterval.optimizeIntervals(regions.iterator.map(l => newQueryInterval(l)).toArray)
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

  /** Close an optional wrapped closeable and release the SAM reader. */
  override def close(): Unit = {
    this.closer.foreach(_.close())
    this.reader.close()
  }

  /**
    * Returns the underlying SamReader. This should be avoided as much as possible, and the
    * SamSource should not be used again after calling [[toSamReader]].
    */
  def toSamReader: SamReader = reader
}
