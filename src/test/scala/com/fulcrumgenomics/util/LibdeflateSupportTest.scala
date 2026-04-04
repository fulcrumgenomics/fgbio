/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.testing.{ErrorLogLevel, SamBuilder, UnitSpec}
import htsjdk.samtools.BAMRecordCodec
import htsjdk.samtools.util.{BlockCompressedInputStream, BlockCompressedOutputStream}
import htsjdk.samtools.util.zip.{DeflaterFactory, InflaterFactory}
import org.scalatest.Retries
import org.scalatest.tagobjects.Retryable

import java.nio.file.Path
import java.util.zip.{Deflater, Inflater}
import scala.util.Random


/** Tests for [[LibdeflateDeflater]], [[LibdeflateInflater]], [[LibdeflateDeflaterFactory]],
  * [[LibdeflateInflaterFactory]], and [[LibdeflateSupport]], covering unit-level correctness,
  * cross-compatibility with JDK compression, end-to-end BAM round-trips, and benchmarks.
  */
class LibdeflateSupportTest extends UnitSpec with ErrorLogLevel with Retries {
  private val supported = LibdeflateSupport.isSupported
  private val levels    = Seq(2, 5, 9)

  override def withFixture(test: NoArgTest) = {
    if (isRetryable(test))
      withRetry { super.withFixture(test) }
    else
      super.withFixture(test)
  }

  /** Generates random byte data of the given size with some compressibility (mix of random and repeated bytes). */
  private def randomData(size: Int, seed: Int = 42): Array[Byte] = {
    val random = new Random(seed)
    val data = new Array[Byte](size)
    for (i <- data.indices) {
      data(i) = if (i % 3 == 0) random.nextInt(256).toByte else 65.toByte
    }
    data
  }

  /** Decompresses using the JDK Inflater for cross-compatibility verification. */
  private def jdkInflate(compressed: Array[Byte], compressedLen: Int, expectedSize: Int): Array[Byte] = {
    val inflater = new Inflater(true) // nowrap=true for raw DEFLATE
    val output = new Array[Byte](expectedSize)
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, expectedSize)
    inflater.end()
    produced shouldBe expectedSize
    output
  }

  /** Compresses using the JDK Deflater for cross-compatibility verification. Returns (compressed, compressedLen). */
  private def jdkDeflate(input: Array[Byte], level: Int = 5): (Array[Byte], Int) = {
    val deflater = new Deflater(level, true) // nowrap=true for raw DEFLATE
    val output = new Array[Byte](input.length + 256)
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val written = deflater.deflate(output, 0, output.length)
    deflater.finished() shouldBe true
    deflater.end()
    (output, written)
  }

  /** Builds a BAM file with programmatically generated reads for use in benchmark and round-trip tests.
    * Creates 500 read pairs (1000 records) spread across two contigs.
    */
  private def buildTestBam(): Path = {
    val builder = new SamBuilder(readLength = 100)
    val random  = new Random(42)
    Range.inclusive(1, 500).foreach { _ =>
      builder.addPair(contig = random.nextInt(2), start1 = random.nextInt(10000) + 1, start2 = random.nextInt(10000) + 1)
    }
    val bam = makeTempFile("benchmark", ".bam")
    builder.write(bam)
    bam
  }

  // ---------------------------------------------------------------------------
  // LibdeflateDeflater unit tests
  // ---------------------------------------------------------------------------

  Seq(0, 1, 5, 6, 9).foreach { level =>
    s"LibdeflateDeflater at level $level" should "compress data that JDK Inflater can decompress" in {
      if (!supported) cancel("jlibdeflate not available on this platform")
      val input = randomData(4096)
      val deflater = new LibdeflateDeflater(level, nowrap = true)
      val compressed = new Array[Byte](input.length + 1024)

      deflater.reset()
      deflater.setInput(input, 0, input.length)
      deflater.finish()
      val compressedLen = deflater.deflate(compressed, 0, compressed.length)

      deflater.finished() shouldBe true
      compressedLen should be > 0

      val decompressed = jdkInflate(compressed, compressedLen, input.length)
      decompressed shouldBe input
      deflater.end()
    }
  }

  "LibdeflateDeflater" should "handle empty input" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val deflater = new LibdeflateDeflater(5, nowrap = true)
    val compressed = new Array[Byte](64)

    deflater.reset()
    deflater.setInput(new Array[Byte](0), 0, 0)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)

    deflater.finished() shouldBe true
    val decompressed = jdkInflate(compressed, compressedLen, 0)
    decompressed shouldBe empty
    deflater.end()
  }

  it should "report not finished when output buffer is too small" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(4096)
    val deflater = new LibdeflateDeflater(5, nowrap = true)
    val tinyOutput = new Array[Byte](2) // way too small

    deflater.reset()
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val written = deflater.deflate(tinyOutput, 0, tinyOutput.length)

    written shouldBe 0
    deflater.finished() shouldBe false
    deflater.end()
  }

  it should "be reusable after reset" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val deflater = new LibdeflateDeflater(5, nowrap = true)

    for (i <- 1 to 3) {
      val input = randomData(1024, seed = i)
      val compressed = new Array[Byte](input.length + 256)

      deflater.reset()
      deflater.setInput(input, 0, input.length)
      deflater.finish()
      val compressedLen = deflater.deflate(compressed, 0, compressed.length)
      deflater.finished() shouldBe true

      val decompressed = jdkInflate(compressed, compressedLen, input.length)
      decompressed shouldBe input
    }

    deflater.end()
  }

  it should "handle a large input block" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(65536)
    val deflater = new LibdeflateDeflater(5, nowrap = true)
    val compressed = new Array[Byte](input.length + 1024)

    deflater.reset()
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)
    deflater.finished() shouldBe true

    val decompressed = jdkInflate(compressed, compressedLen, input.length)
    decompressed shouldBe input
    deflater.end()
  }

  // ---------------------------------------------------------------------------
  // LibdeflateInflater unit tests
  // ---------------------------------------------------------------------------

  Seq(1, 5, 9).foreach { level =>
    s"LibdeflateInflater at level $level" should "decompress JDK-compressed data" in {
      if (!supported) cancel("jlibdeflate not available on this platform")
      val input = randomData(4096)
      val (compressed, compressedLen) = jdkDeflate(input, level)

      val inflater = new LibdeflateInflater(nowrap = true)
      val output = new Array[Byte](input.length)

      inflater.reset()
      inflater.setInput(compressed, 0, compressedLen)
      val produced = inflater.inflate(output, 0, input.length)

      produced shouldBe input.length
      output shouldBe input
      inflater.end()
    }
  }

  "LibdeflateInflater" should "round-trip with LibdeflateDeflater" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(8192)

    val deflater = new LibdeflateDeflater(5, nowrap = true)
    val compressed = new Array[Byte](input.length + 1024)
    deflater.reset()
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)
    deflater.finished() shouldBe true
    deflater.end()

    val inflater = new LibdeflateInflater(nowrap = true)
    val output = new Array[Byte](input.length)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, input.length)

    produced shouldBe input.length
    output shouldBe input
    inflater.end()
  }

  it should "handle empty input" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val (compressed, compressedLen) = jdkDeflate(new Array[Byte](0))

    val inflater = new LibdeflateInflater(nowrap = true)
    val output = new Array[Byte](0)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, 0)

    produced shouldBe 0
    inflater.end()
  }

  it should "be reusable after reset" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val inflater = new LibdeflateInflater(nowrap = true)

    for (i <- 1 to 3) {
      val input = randomData(1024, seed = i)
      val (compressed, compressedLen) = jdkDeflate(input)
      val output = new Array[Byte](input.length)

      inflater.reset()
      inflater.setInput(compressed, 0, compressedLen)
      val produced = inflater.inflate(output, 0, input.length)

      produced shouldBe input.length
      output shouldBe input
    }

    inflater.end()
  }

  it should "handle a large input block" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(65536)
    val (compressed, compressedLen) = jdkDeflate(input)

    val inflater = new LibdeflateInflater(nowrap = true)
    val output = new Array[Byte](input.length)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, input.length)

    produced shouldBe input.length
    output shouldBe input
    inflater.end()
  }

  // ---------------------------------------------------------------------------
  // Zlib-wrapped (nowrap=false) cross-compatibility tests
  // ---------------------------------------------------------------------------

  "LibdeflateDeflater with nowrap=false" should "produce zlib-wrapped output that JDK Inflater can decompress" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(4096)
    val deflater = new LibdeflateDeflater(5, nowrap = false)
    val compressed = new Array[Byte](input.length + 1024)

    deflater.reset()
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)
    deflater.finished() shouldBe true
    deflater.end()

    // Decompress with JDK Inflater (nowrap=false expects zlib wrapper)
    val jdkInflater = new Inflater(false)
    val output = new Array[Byte](input.length)
    jdkInflater.setInput(compressed, 0, compressedLen)
    val produced = jdkInflater.inflate(output, 0, input.length)
    jdkInflater.end()

    produced shouldBe input.length
    output shouldBe input
  }

  it should "produce zlib-wrapped output at multiple compression levels" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    for (level <- Seq(0, 1, 5, 9)) {
      val input = randomData(2048, seed = level)
      val deflater = new LibdeflateDeflater(level, nowrap = false)
      val compressed = new Array[Byte](input.length + 1024)

      deflater.reset()
      deflater.setInput(input, 0, input.length)
      deflater.finish()
      val compressedLen = deflater.deflate(compressed, 0, compressed.length)
      deflater.finished() shouldBe true
      deflater.end()

      val jdkInflater = new Inflater(false)
      val output = new Array[Byte](input.length)
      jdkInflater.setInput(compressed, 0, compressedLen)
      val produced = jdkInflater.inflate(output, 0, input.length)
      jdkInflater.end()

      produced shouldBe input.length
      output shouldBe input
    }
  }

  "LibdeflateInflater with nowrap=false" should "decompress JDK zlib-wrapped data" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(4096)

    // Compress with JDK Deflater (nowrap=false produces zlib wrapper)
    val jdkDeflater = new Deflater(5, false)
    val compressed = new Array[Byte](input.length + 256)
    jdkDeflater.setInput(input, 0, input.length)
    jdkDeflater.finish()
    val compressedLen = jdkDeflater.deflate(compressed, 0, compressed.length)
    jdkDeflater.finished() shouldBe true
    jdkDeflater.end()

    // Decompress with LibdeflateInflater (nowrap=false)
    val inflater = new LibdeflateInflater(nowrap = false)
    val output = new Array[Byte](input.length)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, input.length)
    inflater.end()

    produced shouldBe input.length
    output shouldBe input
  }

  it should "decompress JDK zlib-wrapped data at multiple compression levels" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    for (level <- Seq(1, 5, 9)) {
      val input = randomData(2048, seed = level)

      val jdkDeflater = new Deflater(level, false)
      val compressed = new Array[Byte](input.length + 256)
      jdkDeflater.setInput(input, 0, input.length)
      jdkDeflater.finish()
      val compressedLen = jdkDeflater.deflate(compressed, 0, compressed.length)
      jdkDeflater.finished() shouldBe true
      jdkDeflater.end()

      val inflater = new LibdeflateInflater(nowrap = false)
      val output = new Array[Byte](input.length)
      inflater.reset()
      inflater.setInput(compressed, 0, compressedLen)
      val produced = inflater.inflate(output, 0, input.length)
      inflater.end()

      produced shouldBe input.length
      output shouldBe input
    }
  }

  "Libdeflate zlib-wrapped round-trip" should "round-trip with nowrap=false" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val input = randomData(8192)

    val deflater = new LibdeflateDeflater(5, nowrap = false)
    val compressed = new Array[Byte](input.length + 1024)
    deflater.reset()
    deflater.setInput(input, 0, input.length)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)
    deflater.finished() shouldBe true
    deflater.end()

    val inflater = new LibdeflateInflater(nowrap = false)
    val output = new Array[Byte](input.length)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, input.length)
    inflater.end()

    produced shouldBe input.length
    output shouldBe input
  }

  it should "handle empty input with nowrap=false" in {
    if (!supported) cancel("jlibdeflate not available on this platform")
    val deflater = new LibdeflateDeflater(5, nowrap = false)
    val compressed = new Array[Byte](64)
    deflater.reset()
    deflater.setInput(new Array[Byte](0), 0, 0)
    deflater.finish()
    val compressedLen = deflater.deflate(compressed, 0, compressed.length)
    deflater.finished() shouldBe true
    deflater.end()

    val inflater = new LibdeflateInflater(nowrap = false)
    val output = new Array[Byte](0)
    inflater.reset()
    inflater.setInput(compressed, 0, compressedLen)
    val produced = inflater.inflate(output, 0, 0)
    inflater.end()

    produced shouldBe 0
  }

  // ---------------------------------------------------------------------------
  // End-to-end BAM round-trip tests
  // ---------------------------------------------------------------------------

  "Libdeflate compression" should "produce a valid BAM that can be read back with JDK inflater" in {
    if (!supported) cancel("jlibdeflate is not available on this platform")

    val testBam = buildTestBam()
    val source  = SamSource(testBam)
    val records = source.toList
    val header  = source.header
    source.close()

    // Write with libdeflate
    val output = makeTempFile("test", ".bam")
    val os     = new BlockCompressedOutputStream(output, 5, new LibdeflateDeflaterFactory)
    val codec  = new BAMRecordCodec(header)
    codec.setOutputStream(os)
    records.foreach(rec => codec.encode(rec.asSam))
    os.close()

    // Read back with JDK inflater
    val is        = new BlockCompressedInputStream(output.toFile, new InflaterFactory)
    val readCodec = new BAMRecordCodec(header)
    readCodec.setInputStream(is, output.toFile.toString)
    val readBack = Iterator.continually(readCodec.decode()).takeWhile(_ != null).toList
    is.close()

    readBack.size shouldBe records.size
    readBack.zip(records).foreach { case (read, original) =>
      read.getReadName shouldBe original.asSam.getReadName
      read.getReadBases shouldBe original.asSam.getReadBases
      read.getBaseQualities shouldBe original.asSam.getBaseQualities
    }
  }

  it should "read a JDK-compressed BAM with libdeflate inflater" in {
    if (!supported) cancel("jlibdeflate is not available on this platform")

    val testBam = buildTestBam()
    val source  = SamSource(testBam)
    val records = source.toList
    val header  = source.header
    source.close()

    // Write with JDK
    val output = makeTempFile("test", ".bam")
    val os     = new BlockCompressedOutputStream(output, 5, new DeflaterFactory)
    val codec  = new BAMRecordCodec(header)
    codec.setOutputStream(os)
    records.foreach(rec => codec.encode(rec.asSam))
    os.close()

    // Read back with libdeflate
    val is        = new BlockCompressedInputStream(output.toFile, new LibdeflateInflaterFactory)
    val readCodec = new BAMRecordCodec(header)
    readCodec.setInputStream(is, output.toFile.toString)
    val readBack = Iterator.continually(readCodec.decode()).takeWhile(_ != null).toList
    is.close()

    readBack.size shouldBe records.size
    readBack.zip(records).foreach { case (read, original) =>
      read.getReadName shouldBe original.asSam.getReadName
      read.getReadBases shouldBe original.asSam.getReadBases
      read.getBaseQualities shouldBe original.asSam.getBaseQualities
    }
  }

  // ---------------------------------------------------------------------------
  // Benchmark tests
  // ---------------------------------------------------------------------------

  "LibdeflateSupport" should "be available" in {
    if (!supported) cancel("jlibdeflate is not available on this platform")
  }

  levels.foreach { level =>
    it should s"deflate faster than the JDK Deflater at level $level" taggedAs Retryable in {
      if (!supported) cancel("jlibdeflate is not available on this platform")
      val testBam = buildTestBam()
      val source  = SamSource(testBam)
      val records = source.toList
      val header  = source.header
      source.close()

      /** Writes all records to a BAM using the given deflater factory and returns elapsed time in ms. */
      def run(factory: DeflaterFactory): Long = {
        val output    = makeTempFile("test", ".bam")
        val startTime = System.currentTimeMillis()
        val os        = new BlockCompressedOutputStream(output, level, factory)
        val codec     = new BAMRecordCodec(header)
        codec.setOutputStream(os)
        val repetitions = levels.max - level
        Range.inclusive(0, repetitions).foreach { _ =>
          records.foreach(rec => codec.encode(rec.asSam))
        }
        os.close()
        System.currentTimeMillis() - startTime
      }

      val jdkTime        = run(new DeflaterFactory)
      val libdeflateTime = run(new LibdeflateDeflaterFactory)

      info(f"Libdeflate: ${libdeflateTime}ms JDK: ${jdkTime}ms speedup: ${jdkTime / libdeflateTime.toFloat}%.2fx")
      libdeflateTime.toDouble should be <= (jdkTime * 1.05)
    }
  }

}
