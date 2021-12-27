/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

import java.io.{InputStream, OutputStream}
import java.nio.file.{Files, Path, Paths}
import java.util.zip.{GZIPInputStream, GZIPOutputStream}

import com.fulcrumgenomics.commons.CommonsDef.DirPath
import com.fulcrumgenomics.commons.io.{IoUtil, PathUtil}
import htsjdk.samtools.util.BlockCompressedOutputStream

/**
  * Provides common IO utility methods.  Can be instantiated to create a custom factory, or
  * the companion object can be used as a singleton version.
  */
class Io(var compressionLevel: Int = 5,
         override val bufferSize: Int = 128*1024,
         var tmpDir: DirPath = Paths.get(System.getProperty("java.io.tmpdir"))) extends IoUtil {

  /** Adds the automatic handling of gzipped files when opening files for writing. */
  override def toOutputStream(path: Path): OutputStream = {
    PathUtil.extensionOf(path) match {
      case Some(".bgz") | Some(".bgzip") => new BlockCompressedOutputStream(super.toOutputStream(path), path.toFile, compressionLevel)
      case _           => super.toOutputStream(path)
    }
  }

  /** Overridden to ensure tmp directories are created within the given tmpDir. */
  override def makeTempDir(name: String): DirPath = Files.createTempDirectory(tmpDir, name)

  /** Overridden to ensure that tmp files are created within the correct tmpDir. */
  override def makeTempFile(prefix: String, suffix: String, dir: Option[DirPath] = Some(tmpDir)): DirPath = super.makeTempFile(prefix, suffix, dir)
}

/** Singleton object that can be used when the default buffer size and compression are desired. */
object Io extends Io(compressionLevel=5, bufferSize=128*1024, Paths.get(System.getProperty("java.io.tmpdir")))
