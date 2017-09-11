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

package com.fulcrumgenomics.testing

import java.nio.file.{Files, Path}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.reference.{FastaSequenceIndexCreator, ReferenceSequenceFileFactory}
import htsjdk.samtools.{SAMFileHeader, SAMSequenceDictionary, SAMSequenceRecord, SAMTextHeaderCodec}

import scala.collection.mutable.ListBuffer

/** Class to programatically build up a reference sequence. */
class ReferenceSetBuilder {
  private val LineLength = 80

  // Class to build up a single reference sequence
  class ReferenceBuilder private[ReferenceSetBuilder](val name: String) {
    private[ReferenceSetBuilder] val bases = new StringBuilder

    /** Adds bases to the reference. */
    def add(s: String, times: Int=1): this.type = {
      forloop (from=0, until=times) { _ => bases.append(s) }
      this
    }
  }

  /** The sequences in order. */
  private val refs = new ListBuffer[ReferenceBuilder]()

  /** Generates a new ReferenceBuilder object in order and returns it. */
  def add(name: String): ReferenceBuilder = {
    this.refs += new ReferenceBuilder(name)
    this.refs.last
  }

  /** Writes the fasta out to a temp file and creates a sequence dictionary alongside it. */
  def toTempFile(deleteOnExit: Boolean = true): PathToFasta = {
    val path = Files.createTempFile("SamRecordSet.", ".fa")
    if (deleteOnExit) path.toFile.deleteOnExit()
    toFile(path, deleteOnExit=deleteOnExit)
    path
  }

  /** Writes the fasta out to a file, and creates a sequence dictionary and fasta index (FAI) file. */
  def toFile(path: Path, deleteOnExit: Boolean = true): Unit = {
    val out = Io.toWriter(path)
    val dict = new SAMSequenceDictionary()
    val header = new SAMFileHeader(dict)

    this.refs.foreach { ref =>
      out.write('>')
      out.write(ref.name)
      out.newLine()
      val bases = ref.bases.toString()
      ref.bases.toString().sliding(LineLength, LineLength).foreach { line =>
        out.write(line)
        out.newLine()
      }

      dict.addSequence(new SAMSequenceRecord(ref.name, bases.length))
    }

    out.close()

    // Create the sequence dictionary
    val dictOut    = PathUtil.replaceExtension(path, ".dict")
    val dictWriter = Io.toWriter(PathUtil.replaceExtension(path, ".dict"))
    new SAMTextHeaderCodec().encode(dictWriter, header)
    dictWriter.close()
    if (deleteOnExit) dictOut.toFile.deleteOnExit()

    // Create the FAI file
    val faiOut = ReferenceSequenceFileFactory.getFastaIndexFileName(path)
    val fai    = FastaSequenceIndexCreator.buildFromFasta(path)
    fai.write(faiOut)
    if (deleteOnExit) faiOut.toFile.deleteOnExit()

  }
}
