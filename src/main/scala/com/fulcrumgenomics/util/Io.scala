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

import java.io.{BufferedWriter, InputStream, OutputStream, OutputStreamWriter}
import java.nio.file.{Files, Path, Paths}
import java.util.zip.{GZIPInputStream, GZIPOutputStream}

import com.fulcrumgenomics.commons.io.{IoUtil, PathUtil}
import com.fulcrumgenomics.commons.CommonsDef.DirPath

/**
  * Provides common IO utility methods.  Can be instantiated to create a custom factory, or
  * the companion object can be used as a singleton version.
  */
class Io(var compressionLevel: Int = 5, override val bufferSize: Int = 128*1024) extends IoUtil {
  /** Adds the automatic handling of gzipped files when opening files for reading. */
  override def toInputStream(path: Path): InputStream = {
    PathUtil.extensionOf(path) match {
      case Some(".gz") => new GZIPInputStream(Files.newInputStream(path), bufferSize)
      case _           => super.toInputStream(path)
    }
  }

  /** Adds the automatic handling of gzipped files when opening files for writing. */
  override def toOutputStream(path: Path): OutputStream = {
    PathUtil.extensionOf(path) match {
      case Some(".gz") => new GZIPOutputStream(Files.newOutputStream(path), bufferSize) { this.`def`.setLevel(compressionLevel) }
      case _           => super.toOutputStream(path)
    }
  }

  /** Returns the system default temporary directory path. */
  def defaultTempDir(): DirPath = Paths.get(System.getProperty("java.io.tmpdir"))
}

/** Singleton object that can be used when the default buffer size and compression are desired. */
object Io extends Io(compressionLevel=5, bufferSize=128*1024)
