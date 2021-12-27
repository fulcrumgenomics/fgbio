/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.cmdline

import java.nio.file.Paths

import com.fulcrumgenomics.FgBioDef.SafelyClosable
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.commons.util.SystemUtil.IntelCompressionLibrarySupported
import com.fulcrumgenomics.testing.UnitSpec
import com.intel.gkl.compression.{IntelDeflaterFactory, IntelInflaterFactory}
import htsjdk.samtools.BAMRecordCodec
import htsjdk.samtools.util._
import htsjdk.samtools.util.zip.{DeflaterFactory, InflaterFactory}


class IntelCompressionTest extends UnitSpec {
  private val intelSupported = IntelCompressionLibrarySupported
  private val testBam        = Paths.get("src/test/resources/com/fulcrumgenomics/bam/estimate_pooling_fractions/HG01583.bam")
  private val levels         = Seq(2, 5, 9)

  "IntelDeflater" should "be available" in {
    if (!intelSupported) cancel("IntelDeflater is not available on this platform")
  }

  levels.foreach { level =>
    it should s"deflate faster than the JDK Deflater on level $level" in {
      if (!intelSupported) cancel("IntelDeflater is not available on this platform")
      val source = SamSource(testBam)
      val records = source.toList
      val header  = source.header
      source.safelyClose()

      // a little method to deflate given an deflater factory
      def run(factory: DeflaterFactory): Long = {
        val output    = makeTempFile("test", ".txt")
        val startTime = System.currentTimeMillis()
        val os        = new BlockCompressedOutputStream(output, level, factory)
        val codec     = new BAMRecordCodec(header)
        codec.setOutputStream(os)
        val repetitions = levels.max - level
        Range.inclusive(0, repetitions).foreach { _ =>
          records.foreach { rec =>
            codec.encode(rec.asSam)
          }
        }
        os.close()
        val endTime = System.currentTimeMillis()
        endTime - startTime
      }

      val jdkTime   = run(new DeflaterFactory)
      val intelTime = run(new IntelDeflaterFactory)

      info(f"Intel: ${intelTime}ms JDK: ${jdkTime}ms speedup: ${jdkTime/intelTime.toFloat}%.2fx")

      intelTime should be <= jdkTime
    }
  }

  "IntelInflater" should "be available" in {
    if (!intelSupported) cancel("IntelDeflater is not available on this platform")
  }

  levels.foreach { level =>
    it should s"inflate faster than the JDK Inflater on level $level" in {
      if (!intelSupported) cancel("IntelInflater is not available on this platform")

      val source = SamSource(testBam)
      val records = source.toList
      val header  = source.header
      source.safelyClose()

      // create a new compressed file
      val output = makeTempFile("test", ".txt")
      val os     = new BlockCompressedOutputStream(output, level, new IntelDeflaterFactory)
      val codec  = new BAMRecordCodec(header)
      codec.setOutputStream(os)
      records.foreach { rec =>
        codec.encode(rec.asSam)
      }
      os.close()

      // a little method to inflate given an inflater factory
      def run(factory: InflaterFactory): Long = {
        val startTime = System.currentTimeMillis()
        Range.inclusive(1, 25).foreach { i =>
          val is    = new BlockCompressedInputStream(output.toFile, factory)
          val codec = new BAMRecordCodec(header)
          codec.setInputStream(is, testBam.toFile.toString)
          while (null != codec.decode()) { () }
          is.close()
        }
        val endTime = System.currentTimeMillis()
        endTime - startTime
      }

      val jdkTime   = run(new InflaterFactory)
      val intelTime = run(new IntelInflaterFactory)

      info(f"Intel: ${intelTime}ms JDK: ${jdkTime}ms speedup: ${jdkTime/intelTime.toFloat}%.2fx")

      intelTime should be <= jdkTime
    }
  }
}
