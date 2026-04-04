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

import com.fulcrumgenomics.jlibdeflate.{LibdeflateCompressor, LibdeflateDecompressor}
import htsjdk.samtools.util.zip.{DeflaterFactory, InflaterFactory}

import java.util.zip.{Deflater, Inflater}
import scala.util.Using


/**
  * A `java.util.zip.Deflater` implementation backed by jlibdeflate for hardware-accelerated DEFLATE compression.
  *
  * Only the methods called by HTSJDK's `BlockCompressedOutputStream.deflateBlock()` are implemented:
  * `reset()`, `setInput()`, `finish()`, `deflate()`, and `finished()`.
  *
  * @param level the compression level (0-12; jlibdeflate supports higher levels than zlib's 0-9)
  * @param nowrap if true, use raw DEFLATE (no zlib/gzip wrapper). Must be true for BGZF.
  */
class LibdeflateDeflater(level: Int, nowrap: Boolean) extends Deflater(level, nowrap) {
  // Release the JDK's native zlib allocation immediately; we manage our own native handle
  super.end()

  private val compressor = new LibdeflateCompressor(level)
  private var inputBuf: Array[Byte] = _
  private var inputOff: Int = 0
  private var inputLen: Int = 0
  private var _finished: Boolean = false

  /** Resets the deflater, clearing any buffered input and marking it as not finished. */
  override def reset(): Unit = {
    inputBuf = null
    inputOff = 0
    inputLen = 0
    _finished = false
  }

  /** Stores the input data to be compressed on the next call to [[deflate]]. */
  override def setInput(input: Array[Byte], off: Int, len: Int): Unit = {
    inputBuf = input
    inputOff = off
    inputLen = len
  }

  /** No-op: jlibdeflate compresses the entire buffer in one call in [[deflate]]. */
  override def finish(): Unit = { }

  /** Compresses the input into the output buffer in a single call. Returns the number of compressed
    * bytes written, or 0 if the output buffer was too small (in which case [[finished]] returns false).
    */
  override def deflate(output: Array[Byte], off: Int, len: Int): Int = {
    val result = compressor.deflateCompress(inputBuf, inputOff, inputLen, output, off, len)
    if (result == -1) {
      _finished = false
      0
    } else {
      _finished = true
      result
    }
  }

  /** Returns true if the last call to [[deflate]] successfully compressed all the input. */
  override def finished(): Boolean = _finished

  /** Releases the native jlibdeflate compressor resources. */
  override def end(): Unit = compressor.close()
}


/**
  * A `java.util.zip.Inflater` implementation backed by jlibdeflate for hardware-accelerated DEFLATE decompression.
  *
  * Only the methods called by HTSJDK's `BlockGunzipper.unzipBlock()` are implemented:
  * `reset()`, `setInput()`, and `inflate()`.
  *
  * HTSJDK reads the uncompressed size from the BGZF block footer and passes it as the `len`
  * parameter to `inflate()`, so we always know the expected output size.
  *
  * @param nowrap if true, use raw DEFLATE (no zlib/gzip wrapper). Must be true for BGZF.
  */
class LibdeflateInflater(nowrap: Boolean) extends Inflater(nowrap) {
  // Release the JDK's native zlib allocation immediately; we manage our own native handle
  super.end()

  private val decompressor = new LibdeflateDecompressor()
  private var inputBuf: Array[Byte] = _
  private var inputOff: Int = 0
  private var inputLen: Int = 0

  /** Resets the inflater, clearing any buffered input. */
  override def reset(): Unit = {
    inputBuf = null
    inputOff = 0
    inputLen = 0
  }

  /** Stores the compressed input data to be decompressed on the next call to [[inflate]]. */
  override def setInput(input: Array[Byte], off: Int, len: Int): Unit = {
    inputBuf = input
    inputOff = off
    inputLen = len
  }

  /** Decompresses the input into the output buffer in a single call. The caller must provide
    * the exact expected uncompressed size as `len` (as HTSJDK does from the BGZF block footer).
    * Returns `len` on success.
    */
  override def inflate(output: Array[Byte], off: Int, len: Int): Int = {
    decompressor.deflateDecompress(inputBuf, inputOff, inputLen, output, off, len)
    len
  }

  /** Releases the native jlibdeflate decompressor resources. */
  override def end(): Unit = decompressor.close()
}


/** Factory that creates [[LibdeflateDeflater]] instances for use with HTSJDK. */
class LibdeflateDeflaterFactory extends DeflaterFactory {
  /** Creates a new [[LibdeflateDeflater]] at the given compression level.
    *
    * @param compressionLevel the compression level (0-12)
    * @param gzipCompatible htsjdk passes true for BGZF; the JDK Deflater(level, nowrap) constructor
    *                       treats this as nowrap=true (raw DEFLATE, no zlib header), which is what BGZF requires
    */
  override def makeDeflater(compressionLevel: Int, gzipCompatible: Boolean): Deflater = {
    new LibdeflateDeflater(compressionLevel, gzipCompatible)
  }
}


/** Factory that creates [[LibdeflateInflater]] instances for use with HTSJDK. */
class LibdeflateInflaterFactory extends InflaterFactory {
  /** Creates a new [[LibdeflateInflater]].
    *
    * @param gzipCompatible htsjdk passes true for BGZF; the JDK Inflater(nowrap) constructor
    *                       treats this as nowrap=true (raw DEFLATE, no zlib header), which is what BGZF requires
    */
  override def makeInflater(gzipCompatible: Boolean): Inflater = {
    new LibdeflateInflater(gzipCompatible)
  }
}


/** Utility object to detect whether jlibdeflate is available on the current platform. */
object LibdeflateSupport {
  /** True if jlibdeflate's native library loads and can successfully compress and decompress data. */
  lazy val isSupported: Boolean = try {
    val input = new Array[Byte](100)
    java.util.Arrays.fill(input, 65.toByte)
    val compressed   = Using.resource(new LibdeflateCompressor())(_.deflateCompress(input))
    val decompressed = Using.resource(new LibdeflateDecompressor())(_.deflateDecompress(compressed, input.length))
    java.util.Arrays.equals(input, decompressed)
  } catch {
    case _: Throwable => false
  }
}
