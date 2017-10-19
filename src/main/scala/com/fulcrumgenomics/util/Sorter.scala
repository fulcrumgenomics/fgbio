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

import java.io._
import java.nio.file.{Files, Path}
import java.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.Writer
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.util.Sorter.{Codec, SortEntry}
import htsjdk.samtools.util.TempStreamFactory

import scala.collection.mutable.ArrayBuffer

/**
  * Companion object for [[Sorter]] that contains various types used.
  */
object Sorter {
  /**
    * A specialized tuple that contains a key object that is [[Comparable]] and a serialized version
    * of the object being sorted in an [[Array]].  The entry itself is comparable but simply
    * delegates comparison to the key object.
    *
    * @param key a comparable key that contains necessary information from the object being sorted
    * @param bytes a serialized form of the object read to be written to disk
    * @tparam B the key type
    */
  case class SortEntry[B <: Comparable[B]](key: B, bytes: Array[Byte]) extends Comparable[SortEntry[B]] {
    override def compareTo(that: SortEntry[B]): Int = this.key.compareTo(that.key)
  }

  /**
    * Trait for classes/objects that are able to encode and object into bytes and decode it
    * from bytes.
    *
    * @tparam A the type of object that is encoded/decoded
    */
  trait Codec[A] {
    /** Encode the object into an array of bytes. */
    def encode(a: A): Array[Byte]
    /** Decode an object from an array of bytes. */
    def decode(bs: Array[Byte], start: Int, length: Int): A
  }
}

/**
  * An implementation of a disk-backed sorting system.  The implementation requires two things:
  *
  * 1. An implementation of Codec that can serialize and deserialize objects
  * 2. A function that creates a [[Ordered]] key object for each object being sorted
  *
  * Both must be thread-safe as they may be invoked across threads without external synchronization
  */
class Sorter[A,B <: Ordered[B]](val maxObjectsInRam: Int,
                                private val codec: Codec[A],
                                private val keyfunc: A => B,
                                private val tmpDir: DirPath = Io.tmpDir) extends Iterable[A] with Writer[A] {
  require(maxObjectsInRam > 1, "Max records in RAM must be at least 2, and probably much higher!")
  private val stash = new Array[SortEntry[B]](maxObjectsInRam)
  private val files = new ArrayBuffer[Path]()
  private var recordsInMemory: Int = 0
  private val _tmpStreamFactory = new TempStreamFactory

  /**
    * An iterator that consumes data from a single tmp file of sorted data and produces
    * an iterator of sorted records.
    */
  private class SortedIterator(stream: InputStream, codec: Codec[A]) extends Iterator[A] with Comparable[SortedIterator] with Closeable  {
    private val dataStream = new DataInputStream(stream)
    private var bytes = new Array[Byte](0)
    private var nextValue : A = _
    private var nextKey   : B = _
    private var closed: Boolean = false
    advance()

    /** Reads the next item from the stream, if one is available */
    private def advance(): Unit = {
      if (closed) {
        clear()
      }
      else {
        val length = try { this.dataStream.readInt() } catch { case ex: EOFException => -1 }
        if (length <= 0) {
          close()
          clear()
        }
        else {
          if (this.bytes.length < length) this.bytes = new Array[Byte](length)
          val read  = dataStream.read(bytes, 0, length)
          this.nextValue = this.codec.decode(bytes, 0, length)
          this.nextKey   = keyfunc(this.nextValue)
        }
      }
    }

    /** Clears the current key and value. */
    private def clear(): Unit = {
        this.nextValue = null.asInstanceOf[A]
        this.nextKey   = null.asInstanceOf[B]
    }

    /** True if a call to [[next()]] will yield a value.  */
    override def hasNext: Boolean = nextValue != null

    /** Retrieve the next value from the iterator. */
    override def next(): A = {
      if (!hasNext) throw new NoSuchElementException("next called on empty Iterator")
      val retval = this.nextValue
      advance()
      retval
    }

    /** The way the iterators are used it is guaranteed that we'll never compare an empty iterator. */
    override def compareTo(that: SortedIterator): Int = this.nextKey.compareTo(that.nextKey)

    /** Closes the underlying stream. */
    def close(): Unit = if (!closed) { this.stream.safelyClose(); closed = true }
  }


  /**
    * An iterator that merges records from [[SortedIterator]]s and maintains ordering.
    */
  private class MergingIterator(val streams: Seq[InputStream]) extends Iterator[A] with Closeable {
    private val iterators = new util.PriorityQueue[SortedIterator]()
    streams.foreach { stream =>
      val iter = new SortedIterator(stream, codec)
      if (iter.hasNext) iterators.add(iter)
    }

    /** True if a call to [[next()]] will yield a value.  */
    override def hasNext: Boolean = !this.iterators.isEmpty

    /** Retrieve the next value from the iterator. */
    override def next(): A = {
      if (!hasNext) throw new NoSuchElementException("next called on empty iterator.")
      val iter = this.iterators.poll()
      val entry = iter.next()
      if (iter.hasNext) this.iterators.add(iter)
      else iter.close()

      entry
    }

    /** Closes all underlying iterators/streams and frees the iterators. */
    override def close(): Unit = {
      this.iterators.foreach(i => i.safelyClose())
      this.iterators.clear()
    }
  }

  /**
    * Adds a record to the sorter.  This is an amortized constant time operation, but the individual times
    * will vary wildly as the accrued objects are written to disk once the max in memory is reached.
    */
  override def write(item: A): Unit = {
    val key   = keyfunc(item)
    val bytes = this.codec.encode(item)
    stash(recordsInMemory) = SortEntry(key, bytes)
    recordsInMemory += 1

    if (recordsInMemory == maxObjectsInRam) spill()
  }

  /**
    * Writes out a temporary file containing all the accrued objects and releases the objects so
    * that they can be garbage collected.
    */
  private def spill(): Unit = {
    if (recordsInMemory > 0) {
      util.Arrays.parallelSort(stash, 0, recordsInMemory)
      val path = Io.makeTempFile("sorter.", ".tmp", dir=Some(this.tmpDir))
      val out = new DataOutputStream(_tmpStreamFactory.wrapTempOutputStream(Io.toOutputStream(path), Io.bufferSize))
      forloop(from = 0, until = recordsInMemory) { i =>
        val bytes = stash(i).bytes
        out.writeInt(bytes.length)
        out.write(bytes)
        stash(i) = null
      }
      out.close()
      path.toFile.deleteOnExit()
      this.files += path
      this.recordsInMemory = 0
    }
  }

  /**
    * Generates an iterator over the data accumulated so far.  May be called multiple times to iterate
    * over the data multiple times, and may even be invoked between adding more objects. Note however
    * that each call to iterator will cause any accumulated objects to be written to disk, so it should
    * not be invoked too frequently!
     */
  def iterator: SelfClosingIterator[A] = {
    spill()
    val streams = files.map(f => _tmpStreamFactory.wrapTempInputStream(Io.toInputStream(f), Io.bufferSize))
    val mergingIterator = new MergingIterator(streams)
    new SelfClosingIterator(mergingIterator, () => mergingIterator.close())
  }

  /** Deletes any temporary files that have been created. */
  def close(): Unit = {
    this.files.foreach(Files.deleteIfExists)
  }
}
