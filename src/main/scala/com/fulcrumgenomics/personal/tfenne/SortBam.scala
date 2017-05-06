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

package com.fulcrumgenomics.personal.tfenne

import java.io.{ByteArrayOutputStream, InputStream}
import java.nio.file.{Files, Path}
import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.personal.tfenne.Sorter.{MergingIterator, SortEntry}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.util.{SortingCollection, TempStreamFactory}
import htsjdk.samtools.{BAMRecordCodec, SAMFileWriterFactory, SAMRecord, SamReaderFactory}

import scala.collection.mutable.ArrayBuffer

@clp(group=ClpGroups.Personal, description =
  """
    |Sorts a BAM into coordinate order.
  """)
class SortBam
( @arg(flag="i", doc="Input SAM or BAM.")   val input: PathToBam,
  @arg(flag="o", doc="Output SAM or BAM.")  val output: PathToBam,
  @arg(flag="m", doc="Max records in RAM.") val maxRecordsInRam: Int
) extends FgBioTool with LazyLogging {
  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)

    val in     = SamReaderFactory.make().open(input)
    val header = in.getFileHeader.clone()
    header.setSortOrder(SortOrder.coordinate)
    val sorter = new Sorter(maxRecordsInRam, new BAMRecordCodec(header), BamCoordinateKey.apply)
    val out    = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(header, true, output.toFile, null)

    val sortProgress  = new ProgressLogger(logger, verb="sorted", unit=5e6.toInt)

    in.foreach { rec =>
      sorter += rec
      sortProgress.record(rec)
    }

    logger.info("Writing.")
    val writeProgress = new ProgressLogger(logger, verb="wrote", unit=5e6.toInt)
    sorter.foreach { rec =>
      out.addAlignment(rec)
      writeProgress.record(rec)
    }

    in.safelyClose()
    out.close()
    sorter.close()
  }
}

object BamCoordinateKey {
  def apply(rec: SAMRecord): BamCoordinateKey = {
    val ref = rec.getReferenceIndex.toInt
    val idx = if (ref < 0) Int.MaxValue else ref
    new BamCoordinateKey(idx, rec.getAlignmentStart, rec.getFlags)
  }
}

case class BamCoordinateKey(refIndex: Int, pos: Int, flags: Int) extends Comparable[BamCoordinateKey] {
  override def compareTo(that: BamCoordinateKey): Int = {
    var retval = this.refIndex - that.refIndex
    if (retval == 0) retval = this.pos - that.pos
    if (retval == 0) retval = this.flags - that.flags
    retval
  }
}

object Sorter {
  case class SortEntry[B <: Comparable[B]](key: B, bytes: Array[Byte]) extends Comparable[SortEntry[B]] {
    override def compareTo(that: SortEntry[B]): Int = this.key.compareTo(that.key)
  }

  case class MergeEntry[A,B <: Comparable[B]](rec: A, key: B) extends Comparable[MergeEntry[A,B]] {
    override def compareTo(that: MergeEntry[A, B]): Int = this.key.compareTo(that.key)
  }

  /** An iterator that consumes from a single stream for a single tmp sort file. */
  private class SingleStreamIterator[A,B <: Comparable[B]](stream: InputStream,
                                                           codec: SortingCollection.Codec[A],
                                                           keyfunc: A => B
                                                          ) extends Iterator[MergeEntry[A,B]] with Comparable[SingleStreamIterator[A,B]] {
    codec.setInputStream(stream)
    private var nextEntry = decodeNext

    private def decodeNext: MergeEntry[A,B] = {
      val rec = this.codec.decode()
      if (rec == null) null else MergeEntry(rec, keyfunc(rec))
    }

    override def hasNext: Boolean = nextEntry != null

    override def next(): MergeEntry[A, B] = {
      if (nextEntry == null) throw new NoSuchElementException("next called on empty FileIterator.")
      val retval = this.nextEntry
      this.nextEntry = decodeNext
      retval
    }

    override def compareTo(that: SingleStreamIterator[A, B]): Int = this.nextEntry.key.compareTo(that.nextEntry.key)

    def close(): Unit = this.stream.safelyClose()
  }

  /** An iterator that merges sorted streams on the fly. */
  private class MergingIterator[A,B <: Comparable[B]](val streams: Seq[InputStream],
                                                      codec: SortingCollection.Codec[A],
                                                      val keyfunc: A => B) extends Iterator[A] {
    private val iterators = new util.TreeSet[SingleStreamIterator[A,B]]()
    streams.foreach { stream =>
      val iter = new SingleStreamIterator[A,B](stream, copyOf(codec), keyfunc)
      if (iter.hasNext) iterators.add(iter)
    }

    /* This method exists because the scala compiler doesn't like the public clone() override in codec. */
    private def copyOf(codec: SortingCollection.Codec[A]):  SortingCollection.Codec[A] = {
      codec.getClass.getMethod("clone").invoke(codec).asInstanceOf[SortingCollection.Codec[A]]
    }

    override def hasNext: Boolean = !this.iterators.isEmpty

    override def next(): A = {
      if (!hasNext) throw new NoSuchElementException("next called on empty iterator.")
      val iter = this.iterators.pollFirst()
      val entry = iter.next()
      if (iter.hasNext) this.iterators.add(iter)
      else iter.close()

      entry.rec
    }
  }
}

/** An implementation of a disk-backed sorting system. */
class Sorter[A,B <: Comparable[B]](val maxRecordsInRam: Int,
                                   private val codec: SortingCollection.Codec[A],
                                   private val keyfunc: A => B) extends Iterable[A] {
  private val stash = new Array[SortEntry[B]](maxRecordsInRam)
  private val files = new ArrayBuffer[Path]()
  private var recordsInMemory: Int = 0
  private val byteStream = new ByteArrayOutputStream(128 * 1024) // TODO: how to size this?
  private val _tmpStreamFactory = new TempStreamFactory

  codec.setOutputStream(byteStream)

  def +=(rec: A): Unit = {
    val key = keyfunc(rec)
    this.codec.encode(rec)
    val bytes = this.byteStream.toByteArray
    this.byteStream.reset()
    stash(recordsInMemory) = SortEntry(key, bytes)
    recordsInMemory += 1

    if (recordsInMemory == maxRecordsInRam) spill()
  }

  private def spill(): Unit = {
    if (recordsInMemory > 0) {
      util.Arrays.parallelSort(stash, 0, recordsInMemory)
      val path = Io.makeTempFile("sorter.", ".tmp")
      val out = _tmpStreamFactory.wrapTempOutputStream(Io.toOutputStream(path), Io.bufferSize)
      forloop(from = 0, until = recordsInMemory) { i =>
        out.write(stash(i).bytes)
        stash(i) = null
      }
      out.close()
      path.toFile.deleteOnExit()
      this.files += path
      this.recordsInMemory = 0
    }
  }

  def iterator: Iterator[A] = {
    spill()
    val streams = files.map(f => _tmpStreamFactory.wrapTempInputStream(Io.toInputStream(f), Io.bufferSize))
    new MergingIterator(streams, codec, keyfunc)
  }

  def close(): Unit = {
    this.files.foreach(Files.deleteIfExists)
  }
}
