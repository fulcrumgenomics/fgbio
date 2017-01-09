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

import java.io._
import java.nio.channels.{Channels, ReadableByteChannel, WritableByteChannel}
import java.nio.file.Path

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import dagr.sopt.cmdline.ValidationException

@clp(description =
  """
    |Performs parallelized block compression or decompression (BGZF).
    |
    |This tool compresses and decompresses data into BGZF compression format (see
    |https://samtools.github.io/hts-specs/SAMv1.pdf).
    |
    |Threads are create for reading and writing data respectively, and the number of compression and decompression
    |threads are specified via an option.   This tool attempts to allocate as much memory as possible to store the
    |data while it is being compressed or decompressed.  Use the scaling factor to reduce the amount of memory that
    |should be used.
    |
    |The tool will take a few seconds to being processing data while it parses the command line options and starts all
    |the relevant threads.
  """,
  group = ClpGroups.Utilities)
class ZipBlocks
(@arg(flag="i", doc = "Input file.") var input: Path = Io.StdIn,
 @arg(flag="o", doc = "Output file.") var output: Path = Io.StdOut,
 @arg(flag="d", doc = "Decompress input if true, compress otherwise.") var decompress: Boolean = false,
 @arg(flag="n", doc = "Number of compression/decompression threads.") var numThreads: Int = Zipper.systemCores,
 @arg(flag="b", doc = "The scaling factor for the maximum memory used for blocks storing compressed and uncompressed data.") var blockScalingFactor: Double = Zipper.DefaultBlockScalingFactor
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  if (blockScalingFactor <= 0) throw new ValidationException("block scaling factor must be greater than zero")

  override def execute(): Unit = {
    val numBlocks = Zipper.numBlocks(blockScalingFactor = blockScalingFactor)

    val readableByteChannel: ReadableByteChannel = {
      if (input == Io.StdIn) Channels.newChannel(System.in)
      else new FileInputStream(input.toFile).getChannel
    }

    val writableByteChannel: WritableByteChannel = {
      if (output == Io.StdOut) Channels.newChannel(System.out)
      else new FileOutputStream(output.toFile).getChannel
    }

    val zipper = new Zipper(
      decompress = decompress,
      in         = readableByteChannel,
      out        = writableByteChannel,
      numThreads = numThreads,
      numBlocks  = numBlocks
    )

    val zipperThread = new Thread(zipper)

    // go for it
    zipperThread.start()
    zipperThread.join()

    // shut it down
    zipper.close()
  }
}