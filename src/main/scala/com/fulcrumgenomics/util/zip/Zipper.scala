/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.util.zip

import java.io.IOException
import java.lang.Thread.UncaughtExceptionHandler
import java.nio.{ByteBuffer, ByteOrder}
import java.nio.channels.{Channels, ReadableByteChannel, WritableByteChannel}
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicBoolean
import java.util.zip
import java.util.zip.CRC32

import com.fulcrumgenomics.util.CloseableBlockingQueue
import dagr.commons.util.LazyLogging
import htsjdk.samtools.SAMFormatException
import htsjdk.samtools.util.{BinaryCodec, BlockCompressedStreamConstants, BlockGunzipper}
import htsjdk.samtools.util.BlockCompressedStreamConstants._
import htsjdk.samtools.util.zip.DeflaterFactory
import oshi.SystemInfo
import oshi.hardware.platform.mac.MacHardwareAbstractionLayer

import scala.collection.mutable
import scala.collection.mutable.ListBuffer


/** Allows a flag to be set to stop the given runnable. Implementers may choose to stop if the flag is set, or not. */
trait StopabbleRunnable extends Runnable {
  private val _isStopped = new AtomicBoolean(false)
  def stopIt(): Unit = _isStopped.set(true)
  def isStopped: Boolean = _isStopped.get()
}

private[zip] object Zipper {
  val TimeoutInMilliseconds: Int = 1
  val TimeoutUnit: TimeUnit = TimeUnit.MILLISECONDS

  private val hal = new SystemInfo().getHardware

  /** Total number of cores in the system */
  val systemCores : Int =
    if (hal.isInstanceOf[MacHardwareAbstractionLayer])
      this.hal.getProcessor.getPhysicalProcessorCount
    else
      this.hal.getProcessor.getLogicalProcessorCount

  /** The heap size of the JVM. */
  val heapSize: Long = Runtime.getRuntime.maxMemory

  val DefaultBlockScalingFactor: Double = 0.75

  /** Calculates the number of blocks we should use when running `ZipBlocks` */
  def numBlocks(blockScalingFactor: Double = DefaultBlockScalingFactor): Int = {
    val blockSize = Blocks.sizeInBytes
    val availableRamForBlocks = heapSize * blockScalingFactor
    if (availableRamForBlocks < blockSize) throw new IllegalStateException(s"Available memory was too small ($availableRamForBlocks < $blockSize)")
    (availableRamForBlocks / blockSize).toInt
  }
}

/** A runnable that orchestrates the multi-threaded compression or decompression of data, reading from an input stream
  * and writing to an output stream.
  *
  * @param decompress true to decompress, false othewise
  * @param in the input stream from which to read
  * @param out the output stream to write compressed or uncompressed data.
  * @param numThreads the number of threads to use to compress or decompress data.
  * @param numBlocks the number of blocks to store in memory
  */
class Zipper(decompress: Boolean,
             in: ReadableByteChannel,
             out: WritableByteChannel,
             numThreads: Int = 1,
             numBlocks: Int = Zipper.numBlocks()) extends StopabbleRunnable with LazyLogging {

  // Developer note: each block potentially stores 2x the amount of bytes we really need, so we could instead have a
  // pool of byte arrays from which we draw when need.  Ignoring for now.

  // TODO:
  // - output the compression ratio, raw bytes in and out.
  // - output timing stats

  /** The queue of blocks that can be used to store new data.  The reader will retrieve blocks from this queue before
    * reading data into those blocks. The writer will put blocks into the queue after data contained has been written. */
  private val unusedQueue: CloseableBlockingQueue[Blocks] = new CloseableBlockingQueue[Blocks]()
  /** The queue of blocks that have input data, and are ready to be compressed or decompressed.  The reader will put
    * read data into blocks and store them in this queue, and consumers will consume blocks from this queue. */
  private val toConsumeQueue: CloseableBlockingQueue[Blocks] = new CloseableBlockingQueue[Blocks]()
  /** The queue of blocks that have output data, and are ready to be written to the output.  The consumers will put
    * data that has been compressed or decompressed in this queue, and the writer will consume blocks from this queue. */
  private val toWriteQueue: CloseableBlockingQueue[Blocks] = new CloseableBlockingQueue[Blocks]()

  private val reader = new Reader(
    in                = in,
    decompress        = decompress,
    inQueue           = unusedQueue,
    outQueue          = toConsumeQueue,
    numBlocksToBuffer = numBlocks
  )

  private val consumers: Seq[Consumer] = Seq.range(0, numThreads, 1).map { i =>
    val consumer    = Consumer(
      compress      = !decompress,
      inQueue       = toConsumeQueue,
      outQueue      = toWriteQueue,
      id = i,
      numThreads    = numThreads)
    consumer
  }

  private val writer = Writer(
    out       = out,
    compress  = !decompress,
    inQueue   = toWriteQueue,
    outQueue  = unusedQueue
  )

  /** Some exception from a child thread. */
  private var threadException: Option[Throwable] = None

  /** The exception handler for child threads. Any exception will stop this thread and store the exception in
    * `threadException`. */
  private val threadExceptionHandler = new UncaughtExceptionHandler {
    def uncaughtException(th: Thread, ex: Throwable): Unit = {
      threadException = Some(ex)
      stopIt()
    }
  }

  /** Create all the blocks upfront that we will use throughout. */
  Stream.range(0, numBlocks, 1).foreach { i => unusedQueue.put(new Blocks(blockId = i, compress = !decompress)) }

  /** Returns true if `isStopped` has been set, or all data has been read in, consumed, and written out. */
  private def isDone: Boolean = {
    val ret = if (isStopped)                    true // we've been stopped
    else {
      if (!reader.isStopped)                    false // reader has not completed
      else if (toConsumeQueue.nonEmpty)         false // more blocks to consume
      else if (toWriteQueue.nonEmpty)           false // more blocks to write
      else if (unusedQueue.size() != numBlocks) false // not all blocks written out (since they should be returned)
      else                                      true  // all done!
    }
    ret
  }

  def run(): Unit = {
    // Create the threads
    val threads = (Seq(reader, writer) ++ consumers).map { runnable => new Thread(runnable) }
    threads.foreach { thread => thread.setUncaughtExceptionHandler(threadExceptionHandler) }

    // start your engines
    threads.foreach { thread => thread.start() }

    while (!isDone) {
      // let others do stuff for a bit
      Thread.sleep(Zipper.TimeoutInMilliseconds)

      // if the reader reaches an EOF, close the toConsumeQueue only if it is empty
      if (!toConsumeQueue.isClosed && reader.isStopped && toConsumeQueue.isEmpty) toConsumeQueue.close()

      // if all the consumers are done, close the toWriteQueue
      if (!toWriteQueue.isClosed && !consumers.exists(c => !c.isStopped)) toWriteQueue.close()

      /*
      val total = toConsumeQueue.size() + toWriteQueue.size() + unusedQueue.size() + reader.blocksBuffered
      logger.debug(" (U=" + unusedQueue.size() + ", C=" + toConsumeQueue.size() + ", W=" + toWriteQueue.size() +
        ", R=" + reader.blocksBuffered + ")"
      )
      */
    }

    // throw the exception if one was raised in a thread
    threadException.map { thr => throw thr }

    // stop all runnables
    require(reader.isStopped) // should be done, so need to stop it
    consumers.foreach { consumer => consumer.stopIt() }
    writer.stopIt()

    // close the queues
    unusedQueue.close()
    toConsumeQueue.close()
    toWriteQueue.close()

    // join
    threads.foreach { thread => thread.join() } // TODO: timeout?

    // forces the writer to perform any last steps (ex. write the last empty GZIP block).
    writer.finish()

    // stop itself
    stopIt()
  }

  /** Closes the input and output channels, and any other streams opened using them. */
  def close(): Unit = {
    reader.close()
    writer.close()
  }
}

/** Companion methods for input and output blocks of data. */
private[zip] object Blocks {
  /** Returns the approximate size in bytes for Blocks */
  val sizeInBytes: Long = {
    SequentialBlocks.sizeInBytes(true) + SequentialBlocks.sizeInBytes(false)
  }

  /** The default block length in bytes of uncompressed data. */
  val DefaultUncompressedFixedBlockLength = DEFAULT_UNCOMPRESSED_BLOCK_SIZE

  /** The maximum block length in bytes of compressed data (ignores the GZIP header). */
  val MaxCompressedFixedBlockLengthNoHeader = MAX_COMPRESSED_BLOCK_SIZE - BLOCK_HEADER_LENGTH

  /** The maximum block length in bytes of compressed data (with the GZIP header and footer). */
  val MaxCompressedFixedBlockLength = MAX_COMPRESSED_BLOCK_SIZE - 2 // subtract two for the BGZF extra field

  /** Companion definitions and methods for sequential blocks. */
  private[zip] object SequentialBlocks {

    /** The default # of blocks in a sequential block. */
    private[zip] var DefaultCapacity: Int = 128

    /** Returns the approximate size in bytes for Blocks. */
    def sizeInBytes(compress: Boolean, capacity: Int = DefaultCapacity): Long = {
      if (compress) MaxCompressedFixedBlockLength * capacity
      else DefaultUncompressedFixedBlockLength * capacity
    }
  }

  /**
    * Stores multiple blocks to compress or decompress.  This facilitates more efficient reading and writing
    * of the data, since more data can be read or written at any one time.  It also allows for multiple blocks to
    * be compressed or decompressed by a single consumer.
    *
    * @param capacity the capacity of fixed-size blocks.
    * @param compress true if the data will be compressed, false if it is to be decompressed.
    */
  private class SequentialBlocks private(capacity: Int, compress: Boolean) {
    def this(compress: Boolean) = {
      this(capacity=SequentialBlocks.DefaultCapacity, compress=compress)
    }

    /** The buffer storing the data. */
    val byteBuffer = {
      if (compress) ByteBuffer.allocate(MaxCompressedFixedBlockLength * capacity)
      else ByteBuffer.allocate(DefaultUncompressedFixedBlockLength * capacity)
    }
    byteBuffer.order(ByteOrder.LITTLE_ENDIAN)

    /** The number of blocks stored. */
    var numStoredBlocks: Int = 0

    /** The size of each stored block. */
    var blockSizes = new ListBuffer[Int]() // NB: for uncompressed data, only the last size can be different

    /** Gets the underlyng array of bytes. */
    def array: Array[Byte] = this.byteBuffer.array()

    /** Clears the underlying bytes, so this object can be re-used. */
    def clear(): Unit = {
      this.byteBuffer.clear()
      blockSizes.clear()
      numStoredBlocks = 0
    }
  }
}

/** Stores the input data to be processed, and the output data after processing.
  *
  * @param blockId the block identifier, typically the index of the block as it was read.
  * @param compress true if the input data will be compressed, false otherwise.
  */
private class Blocks(var blockId: Long, compress: Boolean) extends Ordered[Blocks] {
  import Blocks._

  private val inData = new SequentialBlocks(compress = !compress)
  private val outData = new SequentialBlocks(compress = compress)

  // convenience methods for the input data
  def inByteBuffer: ByteBuffer = inData.byteBuffer
  def inBlockSizes: ListBuffer[Int] = inData.blockSizes
  def inArray: Array[Byte] = inData.array
  def incrementNumInBlocks(): Unit = inData.numStoredBlocks += 1

  // convenience methods for the output data
  def outByteBuffer: ByteBuffer = outData.byteBuffer
  def outBlockSizes: ListBuffer[Int] = outData.blockSizes
  def outArray: Array[Byte] = outData.array
  def incrementNumOutBlocks(): Unit = outData.numStoredBlocks += 1

  /** Clears the underling input and output data. */
  def clear(): Unit = {
    this.inData.clear()
    this.outData.clear()
  }

  /** Compares two blocks by their id. */
  def compare(that: Blocks): Int = this.blockId.compare(that.blockId)
}

/** Reads data from an input byte channel, obtains blocks ] from the input queue to store the read data, and
  * subsequently stores those blocks in the output queue.
  *
  * @param in the input byte channel from which to read data
  * @param inQueue the input queue from which to blocks are obtain to store data read in.
  * @param outQueue the output queue to which to enqueue blocks read in.
  * @param numBlocksToBuffer the number of blocks to internally buffer before requiring they be enqueued in the output
  *                          queue.
  */
private class Reader(val in: ReadableByteChannel,
                     decompress: Boolean,
                     inQueue: CloseableBlockingQueue[Blocks],
                     outQueue: CloseableBlockingQueue[Blocks],
                     numBlocksToBuffer: Int
                    ) extends StopabbleRunnable with LazyLogging {
  private val blocksToFill = new ListBuffer[Blocks]()
  private var blocksFilled = new ListBuffer[Blocks]()

  //def blocksBuffered: Int = blocksToFill.length

  /** The identifier of the next block that will be read in. */
  private var nextBlockId: Long = 0L

  def run(): Unit = {
    var numRead = 0
    while (0 <= numRead && !isStopped) {
      // grab as many blocks from the input queue as we can
      inQueue.takeToMax(blocksToFill, numBlocksToBuffer - blocksToFill.length + blocksFilled.length)

      // read into the buffer of blocks until we reach the end of the file or all blocks have been filled with data.
      var blocksIdx = 0
      numRead = 1 // so we can enter the loop
      // while we have more blocks to read into, the last time we read we returned data, and we are not to stop
      while (blocksIdx < blocksToFill.length && 0 < numRead && !isStopped) {
        // get a block into which to read bytes
        val blocks = blocksToFill(blocksIdx)
        blocks.blockId = nextBlockId
        nextBlockId += 1
        // read into the block if it has bytes remaining in the block
        if (blocks.inByteBuffer.hasRemaining) {
          // keep reading while we data was read and there are bytes remaining in the block
          while (0 < numRead && blocks.inByteBuffer.hasRemaining) {
            // try and read more bytes into this block
            numRead = in.read(blocks.inByteBuffer)
            // NB: do not call `stopIt` since that is meant for extraordinary measures, not normal things like reaching
            // the end of the stream.
          }
          // only move to the next block if we filled the data in this block
          if (0 < numRead) require(!blocks.inByteBuffer.hasRemaining)
        }
        blocksIdx += 1
      }

      // add blocks that have been filled to `blocksFilled`
      if (0 < blocksIdx) {
        blocksFilled ++= blocksToFill.take(blocksIdx)
        if (blocksIdx == blocksToFill.length) blocksToFill.clear()
        else blocksToFill.remove(0, blocksIdx)
      }

      // try and enqueue blocks from `blocksFilled`
      if (blocksFilled.nonEmpty) {
        var added = 0
        if (numBlocksToBuffer == blocksFilled.length) { // we must add some
          // stop if we were not able to add at least one
          added = outQueue.putSome(blocksFilled, min=1)
          if (added == 0) stopIt()
        }
        else {
          added = outQueue.putSome(blocksFilled, min=0)
        }
        if (added == blocksFilled.length) blocksFilled.clear()
        else blocksFilled.remove(0, added)
      }
    }

    // if we were stopped, put all the blocks filled back into the input queue, otherwise we reached the end of the file,
    // so move all the filled blocks to the output queue.
    if (isStopped) {
      if (blocksFilled.nonEmpty && !inQueue.putAll(blocksFilled)) {
        throw new IllegalStateException("Could not add all the blocks read into the input queue")
      }
    }
    else {
      if (blocksFilled.nonEmpty && !outQueue.putAll(blocksFilled)) {
        throw new IllegalStateException("Could not add all the blocks read into the output queue")
      }
      blocksFilled.clear()
    }

    // return all blocks we were buffering into the input queue. Our parent thread uses this to detect completion.
    if (blocksToFill.nonEmpty && !inQueue.putAll(blocksToFill)) {
      throw new IllegalStateException("Could not add all the blocks read into the input queue")
    }
    blocksToFill.clear()

    stopIt()
  }

  /** Closes the input channel. */
  def close(): Unit = {
    in.close()
  }
}

private object Consumer {
  def apply(compress: Boolean,
            inQueue: CloseableBlockingQueue[Blocks],
            outQueue: CloseableBlockingQueue[Blocks],
            id: Int,
            numThreads: Int): Consumer = {
    if (compress) new Deflater(inQueue=inQueue, outQueue=outQueue, id=id, numThreads=numThreads)
    else new Inflator(inQueue=inQueue, outQueue=outQueue, id=id, numThreads=numThreads)
  }

  /** The maximum number of blocks to consume in a batch. */
  private val maxNumBlocks = 8
}

/** Takes blocks from the input queue, compresses or decompresses them, and puts them in the output queue.
  *
  * @param inQueue the queue of blocks to consume (compress or decompress).
  * @param outQueue the queue of blocks to which consumed blocks are enqueued.
  * @param id the unique id of this thread (typically its zero-based index).
  * @param numThreads the total number of consumer threads created.
  */
private abstract class Consumer(inQueue: CloseableBlockingQueue[Blocks],
                                outQueue: CloseableBlockingQueue[Blocks],
                                id: Int,
                                numThreads: Int) extends StopabbleRunnable with LazyLogging {
  /** The number of blocks consumed. */
  var numBlocksConsumed: Int = 0

  def run(): Unit = {
    val buffer = new ListBuffer[Blocks]()

    while (!isStopped && !inQueue.isClosed) {
      // get a block from the input queue
      val numTaken = inQueue.takeToMax(buffer, Consumer.maxNumBlocks)
      if (numTaken == 0) {
        // the input queue must have been closed by our parent, so lets just stop.
        stopIt()
      }
      else {
        // consume the blocks
        buffer.foreach { blocks =>
          consume(blocks)
          numBlocksConsumed += 1
        }
        // output the consumed blocks
        outQueue.putAll(buffer)
        buffer.clear()
      }
    }
    stopIt()
  }

  /** Consume (compress or decompress) a block of data. */
  protected def consume(blocks: Blocks): Unit
}

/** Deflates blocks from the input queue and places them in the output queue.
  *
  * A given block will first be deflated using a deflater at the given compression level.  In some cases, the
  * compressed data may be larger than the uncompressed data.  In this case, a deflater with no compression will be
  * used.  We may still produce compressed data larger than the uncompressed data.  In this case, split the bytes to
  * be compressed into two blocks.  This last case occurs very infrequently, and typically only on random data.
  *
  * @param inQueue the queue of blocks to consume (compress or decompress).
  * @param outQueue the queue of blocks to which consumed blocks are enqueued.
  * @param id the unique id of this thread (typically its zero-based index).
  * @param numThreads the total number of consumer threads created.
  * @param compressionLevel the compression level used to defaulte.
  * @param deflaterFactory the factory to create deflaters.
  */
private class Deflater(inQueue: CloseableBlockingQueue[Blocks],
                       outQueue: CloseableBlockingQueue[Blocks],
                       id: Int,
                       numThreads: Int,
                       compressionLevel: Int = BlockCompressedStreamConstants.DEFAULT_COMPRESSION_LEVEL,
                       deflaterFactory: DeflaterFactory = new DeflaterFactory
                      ) extends Consumer(inQueue, outQueue, id, numThreads) {
  import Blocks._

  /** The main deflater. */
  private val deflater: zip.Deflater = deflaterFactory.makeDeflater(compressionLevel, true)
  /** The deflater to fall back to inf the main deflater compressed data larger than the uncompressed data. */
  private val noCompressionDeflater: zip.Deflater = new zip.Deflater(java.util.zip.Deflater.NO_COMPRESSION, true)

  protected def consume(blocks: Blocks): Unit = {
    // compress the input data, and store it in the output
    var inOffset = 0
    // consumes one block's worth of data at a time.
    while (inOffset < blocks.inByteBuffer.position()) {
      val remaining = Math.min(blocks.inByteBuffer.position() - inOffset, DefaultUncompressedFixedBlockLength)
      consumeBlock(inOffset, remaining, blocks)
      inOffset += remaining
    }
  }

  /** Consumes one blocks worth of data at a time.
    *
    * @param inputOffset the offset into the input data.
    * @param bytesToCompress the number of bytes to compress.
    * @param blocks the block storing the input data, and storing where the compressed data will be stored.
    */
  private def consumeBlock(inputOffset: Int, bytesToCompress: Int, blocks: Blocks): Unit = {
    val input = blocks.inArray
    val output = blocks.outArray

    // first attempt at compression
    deflater.reset()
    deflater.setInput(input, inputOffset, bytesToCompress)
    deflater.finish()
    var numBytesCompressed = deflater.deflate(output, blocks.outByteBuffer.position, MaxCompressedFixedBlockLengthNoHeader)

    // in the unlikely case that we increased the size of the data, try without compression
    if (!deflater.finished) {
      // compression failed, just go uncompressed, it should fit
      noCompressionDeflater.reset()
      noCompressionDeflater.setInput(input, inputOffset, bytesToCompress)
      noCompressionDeflater.finish()
      numBytesCompressed = noCompressionDeflater.deflate(output, blocks.outByteBuffer.position, MaxCompressedFixedBlockLengthNoHeader)
      // if it still doesn't fit, try splitting the compressed data into two blocks.
      if (!noCompressionDeflater.finished) {
        // compress the first and second half of the data
        consumeBlock(inputOffset, bytesToCompress/2, blocks)
        consumeBlock(inputOffset + bytesToCompress/2, bytesToCompress - bytesToCompress/2, blocks)
        // NB: do not update position since recursive calls will handle it
      }
      else {
        updatePosition(bytesToCompress, numBytesCompressed, blocks)
      }
    }
    else {
      updatePosition(bytesToCompress, numBytesCompressed, blocks)
    }
  }

  /**
    * Updates the number of bytes stored in the output byte buffer, increments the number of blocks stored in output
    * byte buffer, and stores the size of the uncompressed and compressed bytes for this block.
    *
    * @param bytesToCompress the number of bytes to compress.
    * @param numBytesCompressed the number of bytes after compression.
    * @param blocks stores the input and output data.
    */
  private def updatePosition(bytesToCompress: Int, numBytesCompressed: Int, blocks: Blocks): Unit = {
    blocks.outByteBuffer.position(blocks.outByteBuffer.position + numBytesCompressed)
    blocks.incrementNumOutBlocks()
    blocks.inBlockSizes.append(bytesToCompress)
    blocks.outBlockSizes.append(numBytesCompressed)
  }
}

/** Inflates blocks from the input queue and places them in the output queue.
  *
  * @param inQueue the queue of blocks to consume (compress or decompress).
  * @param outQueue the queue of blocks to which consumed blocks are enqueued.
  * @param id the unique id of this thread (typically its zero-based index).
  * @param numThreads the total number of consumer threads created.
  */
private class Inflator(inQueue: CloseableBlockingQueue[Blocks],
                       outQueue: CloseableBlockingQueue[Blocks],
                       id: Int,
                       numThreads: Int
                      ) extends Consumer(inQueue, outQueue, id, numThreads) {
  private val blockGunzipper: BlockGunzipper = new BlockGunzipper

  protected def consume(blocks: Blocks): Unit = {
    var index = 0
    var compressedOffset   = 0
    var uncompressedOffset = 0
    val compressedBlocksLength = blocks.inByteBuffer.position()
    while (compressedOffset < compressedBlocksLength) {
      val compressedArray   = blocks.inArray
      val uncompressedArray = blocks.outArray

      // get the compressed size
      require(compressedOffset + BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET <= compressedBlocksLength)
      val compressedSize: Int = unpackInt16(compressedArray, compressedOffset + BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET) + 1

      if (compressedSize < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || compressedSize + compressedOffset > compressedBlocksLength) {
        throw new IOException("Unexpected compressed block length: " + compressedSize)
      }
      blocks.inBlockSizes.append(compressedSize)

      // get the uncompressed size
      val uncompressedSize: Int = unpackInt32(compressedArray, compressedOffset + compressedSize - 4) // uncompressed length
      if (uncompressedSize < 0) throw new IllegalStateException("BGZF file has invalid uncompressedLength: " + uncompressedSize)
      blocks.outBlockSizes.append(uncompressedSize)

      // inflate to the output
      blockGunzipper.unzipBlock(uncompressedArray, uncompressedOffset, compressedArray, compressedOffset, compressedSize)

      compressedOffset   += compressedSize
      uncompressedOffset += uncompressedSize
      index += 1
    }
    require(uncompressedOffset == blocks.outBlockSizes.sum)
    blocks.outByteBuffer.position(uncompressedOffset)
  }

  private def unpackInt16(buffer: Array[Byte], offset: Int): Int = {
    (buffer(offset) & 0xFF) | ((buffer(offset + 1) & 0xFF) << 8)
  }

  private def unpackInt32(buffer: Array[Byte], offset: Int): Int = {
    (buffer(offset) & 0xFF) | ((buffer(offset + 1) & 0xFF) << 8) | ((buffer(offset + 2) & 0xFF) << 16) | ((buffer(offset + 3) & 0xFF) << 24)
  }
}

/** Companion methods for the writer. */
private object Writer {
  def apply(out: WritableByteChannel, compress: Boolean, inQueue: CloseableBlockingQueue[Blocks], outQueue: CloseableBlockingQueue[Blocks]): Writer = {
    if (compress) new GzipWriter(out = out, inQueue = inQueue, outQueue = outQueue)
    else new UncompressedWriter(out = out, inQueue = inQueue, outQueue = outQueue)
  }
}

/** Writes data to an output stream.  Blocks are retrieved from the input queue, and after written, placed in the output
  * queue for re-use by the reader.
  *
  * The writer must write blocks in the order they are read, so it assumes that the id of blocks are sequential, and
  * buffers blocks until the block with the next expected id is encountered.
  *
  * @param out the output channel to which data is written.
  * @param inQueue the input queue of blocks with data to write.
  * @param outQueue the output queue to enqueue blocks that have been written.
  */
private abstract class Writer(protected val out: WritableByteChannel,
                              inQueue: CloseableBlockingQueue[Blocks],
                              outQueue: CloseableBlockingQueue[Blocks]
                             ) extends StopabbleRunnable with LazyLogging {
  /** The buffer of blocks to write, sorted by their id. */
  private val blocksToWrite: mutable.TreeSet[Blocks] = new mutable.TreeSet[Blocks]()
  /** Buffer of blocks that can be placed on the output queue. */
  private val blocksToReturn = new ListBuffer[Blocks]()
  /** The identifier of the next block to be written. */
  private var nextBlockId: Long = 0

  def run(): Unit = {
    while (!isStopped && !inQueue.isClosed) {
      val buffer = new ListBuffer[Blocks]()
      // get some blocks, at least 4, but it may be fewer if the queue closes
      inQueue.takeAtLeast(buffer, 4)
      blocksToWrite ++= buffer

      // write some blocks to the output channel
      while (writeBlocks()) {}

      // return blocks to the output queue
      returnBlocks()
    }

    while (writeBlocks()) {}
    returnBlocks()

    stopIt()
  }

  /** Writes a block from `blocksToWrite` if the next block's id matches the expected id.
    *
    * @return true if a block was written, false otherwise.
    */
  private def writeBlocks(): Boolean = {
    blocksToWrite.headOption match {
      case None => false
      case Some(b) if b.blockId > nextBlockId => false
      case Some(b) =>
        if (b.blockId < nextBlockId) throw new IllegalStateException("Out of order may occur") // more than Long.MaxValue!
        write(b)
        blocksToReturn.append(b)
        this.blocksToWrite -= b
        nextBlockId += 1
        true
    }
  }

  /** Enqueues blocks that have been written into the output queue. */
  private def returnBlocks(): Unit = {
    // return blocks
    blocksToReturn.foreach { b => b.clear() }
    if (blocksToReturn.nonEmpty && !outQueue.putAll(blocksToReturn)) {
      throw new IllegalStateException("Could not add all the blocks written into the output queue (isStopped: " + isStopped + " isClosed: " + outQueue.isClosed + ")")
    }
    blocksToReturn.clear()
  }

  /** Writes a block to the output stream. */
  protected def write(blocks: Blocks): Unit

  /** Perform any last operations (ex. write the ending empty GZIP block and closes the output. */
  def finish(): Unit = {}

    /** Close the output channel. */
  def close(): Unit = out.close()
}

/** Writes uncompressed data to the output stream.
  *
  * @param out the output channel to which data is written.
  * @param inQueue the input queue of blocks with data to write.
  * @param outQueue the output queue to enqueue blocks that have been written.
  */
private class UncompressedWriter(out: WritableByteChannel,
                                 inQueue: CloseableBlockingQueue[Blocks],
                                 outQueue: CloseableBlockingQueue[Blocks]
                                ) extends Writer(out, inQueue, outQueue) {

  var numWritten: Long = 0
  protected def write(blocks: Blocks): Unit = {

    // set limit to the current position, and position to zero for writing
    blocks.outByteBuffer.limit(blocks.outByteBuffer.position())
    blocks.outByteBuffer.position(0)
    numWritten += blocks.outByteBuffer.limit()
    out.write(blocks.outByteBuffer)
  }
}

/** Writes gzip blocks to the output stream.
  *
  * @param out the output channel to which data is written.
  * @param inQueue the input queue of blocks with data to write.
  * @param outQueue the output queue to enqueue blocks that have been written.
  */
private class GzipWriter(out: WritableByteChannel,
                         inQueue: CloseableBlockingQueue[Blocks],
                         outQueue: CloseableBlockingQueue[Blocks]
                        ) extends Writer(out, inQueue, outQueue) {

  /** Computes the crc code for a GZIP block. */
  private val crc32: CRC32 = new CRC32
  /** The code to write binary data. */
  private val codec = new BinaryCodec(Channels.newOutputStream(out))

  protected def write(blocks: Blocks): Unit = {
    var blockIndex = 0

    // uncompressed dataa
    var uncompressedOffset     = 0
    val uncompressedArray      = blocks.inArray

    // compressed data
    var compressedOffset       = 0
    val compressedArray        = blocks.outArray
    val compressedByteBuffer   = blocks.outByteBuffer
    val compressedBlocksLength = compressedByteBuffer.position()

    // compress each block of data
    while (compressedOffset < compressedBlocksLength) {
      val uncompressedBlockSize = blocks.inBlockSizes(blockIndex)
      val compressedBlockSize   = blocks.outBlockSizes(blockIndex)

      // compute the crc code for the uncompressed data.
      crc32.reset()
      crc32.update(uncompressedArray, uncompressedOffset, uncompressedBlockSize)

      // write the block
      writeGzipBlock(
        compressedArray  = compressedArray,
        compressedOffset = compressedOffset,
        compressedSize   = compressedBlockSize,
        uncompressedSize = uncompressedBlockSize,
        crc              = crc32.getValue
      )

      // increment the number of blocks deflated, and move to the next block
      blocks.incrementNumOutBlocks()
      uncompressedOffset += uncompressedBlockSize
      compressedOffset   += compressedBlockSize
      blockIndex         += 1
    }
  }

  /**
    * Writes the entire gzip block, assuming the compressed data is stored in compressedBuffer.
    *
    * @param compressedArray the compressed array of bytes.
    * @param compressedOffset the offset of the compressed data in the array.
    * @param compressedSize the number of bytes for the compressed data.
    * @param uncompressedSize the number of bytes for the uncompressed data.
    * @param crc the crc code for the uncompressed data.
    * @return the size of gzip block that was written.
    */
  private def writeGzipBlock(compressedArray: Array[Byte], compressedOffset: Int, compressedSize: Int, uncompressedSize: Int, crc: Long): Int = {
    codec.writeByte(BlockCompressedStreamConstants.GZIP_ID1)
    codec.writeByte(BlockCompressedStreamConstants.GZIP_ID2)
    codec.writeByte(BlockCompressedStreamConstants.GZIP_CM_DEFLATE)
    codec.writeByte(BlockCompressedStreamConstants.GZIP_FLG)
    codec.writeInt(0)
    codec.writeByte(BlockCompressedStreamConstants.GZIP_XFL)
    codec.writeByte(BlockCompressedStreamConstants.GZIP_OS_UNKNOWN)
    codec.writeShort(BlockCompressedStreamConstants.GZIP_XLEN)
    codec.writeByte(BlockCompressedStreamConstants.BGZF_ID1)
    codec.writeByte(BlockCompressedStreamConstants.BGZF_ID2)
    codec.writeShort(BlockCompressedStreamConstants.BGZF_LEN)
    val totalBlockSize: Int = compressedSize + BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH + BlockCompressedStreamConstants.BLOCK_FOOTER_LENGTH
    codec.writeShort((totalBlockSize - 1).toShort)
    codec.writeBytes(compressedArray, compressedOffset, compressedSize)
    codec.writeInt(crc.toInt)
    codec.writeInt(uncompressedSize)
    totalBlockSize
  }

  /** Perform any last operations (ex. write the ending empty GZIP block. */
  override def finish(): Unit = {
    codec.writeBytes(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK)
    super.finish()
  }

  override def close(): Unit = {
    codec.close()
  }
}
