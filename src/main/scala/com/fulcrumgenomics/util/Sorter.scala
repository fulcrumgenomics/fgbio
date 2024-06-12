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
import java.{lang, util}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.util.Sorter.{Codec, SortEntry}
import htsjdk.samtools.util.TempStreamFactory

import java.util.concurrent.{Executors, TimeUnit}
import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag
import scala.tools.nsc.doc.html.HtmlTags.B

/**
  * Companion object for [[Sorter]] that contains various types used.
  */
object Sorter {
  /**
    * A specialized tuple that contains a key object that is [[java.lang.Comparable]] and a serialized version
    * of the object being sorted in an [[scala.Array]].  The entry itself is comparable but simply
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
  * 2. A function that creates a [[scala.Ordered]] key object for each object being sorted
  *
  * Both must be thread-safe as they may be invoked across threads without external synchronization
  */
class Sorter[A : ClassTag, B <: Ordered[B]](val maxObjectsInRam: Int,
                                 private val codec: Codec[A],
                                 private val keyfunc: A => B,
                                 private val tmpDir: DirPath = Io.tmpDir,
                                 private val threads: Int = 10) extends Iterable[A] with Writer[A] with LazyLogging {
  require(maxObjectsInRam > 1, "Max records in RAM must be at least 2, and probably much higher!")
  private val l1StashSize: Int             = (maxObjectsInRam / threads.toDouble).toInt / 10
  private val l2StashSize: Int             = l1StashSize * 10
  private var l1Stash: Array[A]            = newL1Stash
  private var l2Stash: Array[SortEntry[B]] = newL2Stash
  private var l1StashCount: Int            = 0  // Number of records in current L1 stash
  private var l2StashCount: Int            = 0  // Number of records in current L2 stash
  private val stashLock: Object            = new Object
  private val files                        = new ArrayBuffer[Path]()
  private val _tmpStreamFactory            = new TempStreamFactory
  private val executor                     = {
    if (threads <= 1) {
      Executors.newSingleThreadExecutor()
    }

    Executors.newFixedThreadPool(threads)
  }

  @inline private def newL1Stash: Array[A]            = new Array[A](l1StashSize)
  @inline private def newL2Stash: Array[SortEntry[B]] = new Array[SortEntry[B]](l2StashSize)

  private class StashDrainer(val l1: Array[A], val n: Int, forceFlush: Boolean = false) extends Runnable {
    override def run(): Unit = {
      // Encode all the L1 stash items
      val encodeStart = System.currentTimeMillis()
      val items = new Array[SortEntry[B]](n)
      forloop (from=0, until=n) { i =>
        val item  = l1(i)
        val key   = keyfunc(item)
        val bytes = codec.encode(item)
        items(i) = SortEntry(key, bytes)
      }
      val encodeEnd = System.currentTimeMillis()

      // Now lock the L2 stash and write into it, returning the taken l2 stash if it also needs flushing
      val copyStart = System.currentTimeMillis()
      val l2ToDrain = Sorter.this.stashLock.synchronized {
        Array.copy(items, 0, l2Stash, l2StashCount, items.length)
        l2StashCount += items.length

        if (forceFlush || l2StashCount == l2StashSize) {
          val l2       = l2Stash
          val l2Count  = l2StashCount
          l2Stash      = newL2Stash
          l2StashCount = 0

          // Make the file we're going to write to and add it to the list of files
          val path = Io.makeTempFile("sorter.", ".tmp", dir=Some(tmpDir))
          files += path
          path.toFile.deleteOnExit()

          Some(l2, l2Count, path)
        }
        else {
          None
        }
      }
      val copyEnd = System.currentTimeMillis()

      logger.debug(s"Encoding took ${encodeEnd-encodeStart}ms, copying took ${copyEnd-copyStart}ms on thread ${Thread.currentThread().getName}.")

      l2ToDrain.foreach { case (l2, n, path) =>
        val drainStart = System.currentTimeMillis()
        util.Arrays.parallelSort(l2, 0, n)
        val out = new DataOutputStream(_tmpStreamFactory.wrapTempOutputStream(Io.toOutputStream(path), Io.bufferSize))
        forloop(from = 0, until = n) { i =>
          val bytes = l2(i).bytes
          out.writeInt(bytes.length)
          out.write(bytes)
        }
        out.close()
        val drainEnd = System.currentTimeMillis()
        logger.debug(s"Sorting and draining to file took ${drainEnd-drainStart}ms.")
      }
    }
  }

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

  //////////////////////////////////////////////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////////////////////////////////////////////

  /**
    * Adds a record to the sorter.  This is an amortized constant time operation, but the individual times
    * will vary wildly as the accrued objects are written to disk once the max in memory is reached.
    */
  override def write(item: A): Unit = {
    this.l1Stash(l1StashCount) = item
    l1StashCount += 1
    if (l1StashCount == l1StashSize) drain(force=false)
  }

  /**
    * Writes out a temporary file containing all the accrued objects and releases the objects so
    * that they can be garbage collected.
    */
  private def drain(force: Boolean): Unit = {
    if (force || l1StashCount > 0 || l2StashCount > 0) {
      val l1            = this.l1Stash
      val l1Count       = this.l1StashCount
      this.l1Stash      = newL1Stash
      this.l1StashCount = 0

      this.executor.submit(new StashDrainer(l1, l1Count, force))
    }
  }

  /**
    * Generates an iterator over the data accumulated so far.  May be called multiple times to iterate
    * over the data multiple times, and may even be invoked between adding more objects. Note however
    * that each call to iterator will cause any accumulated objects to be written to disk, so it should
    * not be invoked too frequently!
     */
  def iterator: SelfClosingIterator[A] = {
    drain(force=true)
    if (!this.executor.isShutdown()) {
      this.executor.shutdown()
      this.executor.awaitTermination(600, TimeUnit.SECONDS)
    }

    val streams = files.iterator.map(f => _tmpStreamFactory.wrapTempInputStream(Io.toInputStream(f), Io.bufferSize)).toSeq
    val mergingIterator = new MergingIterator(streams)
    new SelfClosingIterator(mergingIterator, () => mergingIterator.close())
  }

  /** Deletes any temporary files that have been created. */
  def close(): Unit = {
    this.files.foreach(Files.deleteIfExists)
  }
}
