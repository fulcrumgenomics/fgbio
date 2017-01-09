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

import java.io.{ByteArrayInputStream, ByteArrayOutputStream}
import java.nio.channels.Channels

import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.{BlockCompressedInputStream, BlockCompressedOutputStream, BlockCompressedStreamConstants}
import org.scalatest.BeforeAndAfterAll

/**
  * Tests for Zipper.
  */
class ZipperTest extends UnitSpec with BeforeAndAfterAll {

  var previousDefaultCapacity = -1
  override def beforeAll(): Unit = {
    previousDefaultCapacity = Blocks.SequentialBlocks.DefaultCapacity
    // keep this small, as we don't need many to test having multiple blocks in one byte array, and we want to
    // test having many `Blocks` (with a capital B).
    Blocks.SequentialBlocks.DefaultCapacity = 4
  }

  override def afterAll(): Unit = {
    Blocks.SequentialBlocks.DefaultCapacity = previousDefaultCapacity
  }

  /** The inputBytesSize of the data to test. */
  private val inputBytesSize = 100000000

  "Zipper" should "compress data" in {
    // create the random input data
    val inputBytes  = new Array[Byte](inputBytesSize)
    scala.util.Random.nextBytes(inputBytes)

    // run the Zipper (compress)
    val outStream = new ByteArrayOutputStream(inputBytesSize)
    val zipper = new Zipper(
      decompress = false,
      in         = Channels.newChannel(new ByteArrayInputStream(inputBytes)),
      out        = Channels.newChannel(outStream),
      numThreads = 1,
      numBlocks  = 4
    )
    zipper.run()
    zipper.close()

    // grab the compressed bytes
    val outputBytes = outStream.toByteArray

    // decompress it with a known good tool
    val decompressInputStream = new BlockCompressedInputStream(new ByteArrayInputStream(outputBytes))
    val decompressedBytes = new Array[Byte](inputBytesSize)
    decompressInputStream.read(decompressedBytes) shouldBe inputBytesSize
    decompressInputStream.close()

    // make sure we get back the original bytes
    java.util.Arrays.equals(inputBytes, decompressedBytes) shouldBe true
  }

  it  should "decompress data" in {
    // the # of bytes to compress

    // create the random input data
    val inputBytes  = new Array[Byte](inputBytesSize)
    scala.util.Random.nextBytes(inputBytes)

    // compress it with a known good tool
    val compressedBytesStream = new ByteArrayOutputStream(inputBytesSize)
    val compressOutputStream = new BlockCompressedOutputStream(compressedBytesStream, null)
    compressOutputStream.write(inputBytes)
    compressOutputStream.close()

    // run the Zipper (decompress)
    val outStream = new ByteArrayOutputStream(inputBytesSize)
    val zipper = new Zipper(
      decompress = true,
      in         = Channels.newChannel(new ByteArrayInputStream(compressedBytesStream.toByteArray)),
      out        = Channels.newChannel(outStream),
      numThreads = 1,
      numBlocks  = 4
    )
    zipper.run()
    zipper.close()

    // grab the decompressed bytes
    val outputBytes = outStream.toByteArray

    // make sure we get back the original bytes
    java.util.Arrays.equals(inputBytes, outputBytes) shouldBe true
  }
}
