/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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
import java.nio.file.Path

import com.fulcrumgenomics.commons.io.AsyncStreamSink
import com.fulcrumgenomics.commons.util.LazyLogging
import htsjdk.samtools.Defaults
import htsjdk.samtools.util.CloserUtil

/** Stores the prediction result from ViennaRNA's RNAfold
  *
  * @param sequence the DNA sequence
  * @param structure the structure in dot-bracket notation
  * @param deltaG the delta-G
  */
case class DnaFoldPrediction(sequence: String, structure: String, deltaG: Double)

object DnaFoldPredictor {
  def apply(viennaRnaInstallDir: File, tm: Double, redirectStderr: Boolean = true): DnaFoldPredictor = {
    new DnaFoldPredictor(viennaRnaInstallDir.toPath, tm, redirectStderr)
  }
}


/**
  * Class that wraps ViennaRNA's RNAFold utility which can be used to estimate the minimum free energy
  * secondary structure of DNA and RNA molecules. When constructed a background process is started
  * running RNAFold, calls to predict() then pipe input and output through RNAFold.
  *
  * @param viennaRnaInstallDir the path to the viennaRNA installation directory
  * @param tm the melting temperature to use
  * @param redirectStderr true to redirect standard error to the calling processes standard error, otherwise it will be
  *                       logged as DEBUG.
  */
class DnaFoldPredictor(viennaRnaInstallDir: Path, tm: Double, redirectStderr: Boolean = true) extends LazyLogging {


  // Build the process, and set up the input and outputs
  private val process: Process = {
    val binary = viennaRnaInstallDir.resolve("bin").resolve("RNAfold")
    val params = viennaRnaInstallDir.resolve("share").resolve("ViennaRNA").resolve("dna_mathews2004.par")
    val args   = Seq(
      binary.toAbsolutePath,
      "--noconv",
      "--noPS",
      "--temp=", tm,
      "--paramFile=" + params.toAbsolutePath
    ).map(_.toString)
    val builder: ProcessBuilder = new ProcessBuilder(args: _*)
    if (redirectStderr) {
      builder.redirectError(ProcessBuilder.Redirect.INHERIT)
    }
    builder.start()
  }
  private val out: BufferedWriter = wrapIOException { new BufferedWriter(new OutputStreamWriter(process.getOutputStream), Defaults.BUFFER_SIZE) }
  private val in: BufferedReader = wrapIOException { new BufferedReader(new InputStreamReader(process.getInputStream), Defaults.BUFFER_SIZE) }
  private val errorPipe: Option[AsyncStreamSink] = wrapIOException {
    if (redirectStderr) None else {
      Some(Io.pipeStream(process.getErrorStream, s => logger.debug("RNAfold: ", s)))
    }
  }

  /** Wraps a block of code, catching an [[IOException]] and wrapping it in a [[RuntimeException]]. */
  private def wrapIOException[T](f: => T): T = {
    try { f }
    catch { case ioe: IOException =>
      throw new RuntimeException(ioe)
    }
  }

  /** Creates the secondary structure prediction for the input sequence. */
  def predict(sequence: String): DnaFoldPrediction = wrapIOException {
    this.out.write(sequence)
    this.out.newLine()
    this.out.flush()

    // Output when using stdin/stdout looks like this:
    // TGAACTCCTCAACCCTCTTCTCATCAGGAGTGATAGTGGCACATTTGACG
    // ((.(((..(..(((((.(....)..))).)).).)))..))......... ( -3.13)

    val seq2 = readLine()
    val result = readLine()

    if (seq2 != sequence) throw new IllegalStateException(f"Return sequence '$seq2' does not match entered sequence '$sequence'")

    val structure = result.substring(0, seq2.length)
    val dg = result.substring(seq2.length()).replace("(", "").replace(")", "").toDouble

    DnaFoldPrediction(seq2, structure, dg)
  }

  /** Reads a line from the input, ignoring warning lines. */
  private def readLine(): String = {
    while (true) {
      val line: String = this.in.readLine()
      if (line == null || !line.startsWith("WARNING")) return line
    }
    null
  }

  /** Kills the underlying process and makes future calls to predict() invalid. */
  def close(): Unit = {
    CloserUtil.close(this.out)
    this.errorPipe.foreach(err => CloserUtil.close(err))
    if (this.process.isAlive) this.process.destroy()
  }
}