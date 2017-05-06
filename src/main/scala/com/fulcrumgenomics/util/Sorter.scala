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

package com.fulcrumgenomics.util

import java.io.OutputStream
import java.util.Comparator

import com.fulcrumgenomics.FgBioDef._
import dagr.commons.reflect.ReflectionUtil
import htsjdk.samtools.util.SortingCollection
import htsjdk.samtools.util.SortingCollection.Codec

import scala.collection.mutable.ListBuffer
import scala.reflect.runtime.universe.TypeTag

object Sorter {
  val DefaultTmpDirs         = Seq(Io.defaultTempDir())
  val DefaultSampleSize      = 5000
  val DefaultInflationFactor = 6.0d
  val DefaultMemory          = Runtime.getRuntime.totalMemory() / 2

  /** Constructs a Sorter that attempts to limit the amount of memory used. */
  def apply[A](codec: Codec[A],
               comparator: Comparator[A],
               memory: Long            = DefaultMemory,
               tmpDirs: Seq[DirPath]   = DefaultTmpDirs,
               sampleSize: Int         = DefaultSampleSize,
               inflationFactor: Double = DefaultInflationFactor)(implicit tag: TypeTag[A]): Sorter[A] = {
    new Sorter(codec, comparator, Some(memory), None, tmpDirs, sampleSize, inflationFactor)
  }

  /** Constructs a Sorter that attempts limits based on the number of records. */
  def apply[A](codec: Codec[A],
               comparator: Comparator[A],
               records: Int,
               tmpDirs: Seq[DirPath])(implicit tag: TypeTag[A]): Sorter[A] = {
    val tmp = if (tmpDirs.isEmpty) DefaultTmpDirs else tmpDirs
    new Sorter(codec, comparator, None, Some(records), tmp, 0, 1)
  }
}


class Sorter[A] private (private val codec:   Codec[A],
                         private val comparator: Comparator[A],
                         private val memory:  Option[Long] = None,
                         private val records: Option[Int] = None,
                         val tmpDirs: Seq[DirPath] = Seq(Io.defaultTempDir()),
                         val sampleSize: Int = 1000,
                         val inflationFactor: Double = 5)(implicit tag: TypeTag[A]) extends Iterable[A] {

  /** Class to help size up how big the records are. */
  private class Sizer {
    val buffer = new ListBuffer[A]()
    val stream = new OutputStream {
      var bytes = 0L
      override def write(b: Int): Unit = bytes += 1
    }
    /** Calculates the average serialized size in bytes of the records sampled. */
    def bytesPerRecord: Int = Math.round(Math.ceil(stream.bytes / buffer.size.toDouble)).toInt
  }

  /** Use reflection to grab the Class[A] for the record type. */
  private val recordType = ReflectionUtil.typeTagToClass[A]

  /** Use the supplied memory limit, or else half of the heap! */
  private val memoryLimit = memory.getOrElse(Runtime.getRuntime.maxMemory() / 2)

  /** Only one of these should be defined at a time. */
  private var collection: Option[SortingCollection[A]] = None
  private var sizer: Option[Sizer] = None

  (memory, records) match {
    case (Some(bs),     None) => this.sizer = Some(new Sizer)
    case (None,     Some(rs)) => this.collection = Some(SortingCollection.newInstance(recordType, codec, comparator, rs))
    case (_, _)               => throw new IllegalArgumentException("Must supply either memory or records but not both.")
  }

  override def iterator: Iterator[A] = (this.sizer, this.collection) match {
    case (Some(s), None) => s.buffer.iterator
    case (None, Some(c)) => c.iterator()
    case _ => unreachable("Only sizer or collection should be defined at any one time.")
  }

  def +=(a: A): Unit = (this.sizer, this.collection) match {
    case (None, Some(c)) =>
      c.add(a)
    case (Some(s), None) =>
      s.buffer += a
      this.codec.setOutputStream(s.stream)
      this.codec.encode(a)

      if (s.buffer.size >= this.sampleSize) {
        System.gc()
        val bytesPerRecord = s.bytesPerRecord * this.inflationFactor
        val records = (this.memoryLimit / bytesPerRecord).toInt
        this.collection = Some(SortingCollection.newInstance(recordType, codec, comparator, records))
        for (c <- this.collection; r <- s.buffer) c.add(r)
        this.sizer = None
      }
    case _ => unreachable("Only sizer or collection should be defined at any one time.")
  }

  def append(as: A*): Unit = as.foreach(+=)
}
