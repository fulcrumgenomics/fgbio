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

package com.fulcrumgenomics.util

import java.nio.file.Path

import dagr.commons.util.LazyLogging

import scala.io.Source
import scala.util.{Success, Try}

/**
  * Object that enables running of R scripts via the Rscript command line too.
  */
object Rscript extends LazyLogging {
  /** The name of the Rscript executable. */
  private val Executable = "Rscript"

  /** Exception class that holds onto an executable's exit/status code. */
  case class RscriptException(status: Int) extends RuntimeException {
    override def getMessage: String = s"Rscript failed with exit code ${status}."
  }

  /** Returns true if R is available and false otherwise. */
  lazy val Available: Boolean = {
    try {
      val process = new ProcessBuilder(Executable, "-e", "require(ggplot2)").redirectErrorStream(true).start()
      process.waitFor() == 0
    }
    catch { case e: Exception => false }
  }

  /** Executes an Rscript from the classpath if Rscript is available. */
  def execIfAvailable(scriptResource: String, args: String*): Try[Unit] =
    if (Available) exec(scriptResource, args:_*) else Success(Unit)

  /** Executes an Rscript from a script stored at a Path if Rscript is available. */
  def execIfAvailable(script: Path, args: String*): Try[Unit] =
    if (Available) exec(script, args:_*) else Success(Unit)

  /** Executes an Rscript from the classpath. */
  def exec(scriptResource: String, args: String*): Try[Unit] =
    Try { writeResourceToTempFile(scriptResource) }.map(path => exec(path, args:_*))

  /** Executes an Rscript from a script stored at a Path. */
  def exec(script: Path, args: String*): Try[Unit] = Try {
    val command = Rscript.Executable +: script.toAbsolutePath.toString +: args
    val process = new ProcessBuilder(command:_*).redirectErrorStream(false).start()
    val pipe1   = Io.pipeStream(process.getErrorStream, logger.info)
    val pipe2   = Io.pipeStream(process.getInputStream, logger.debug)
    val retval  = process.waitFor()
    pipe1.close()
    pipe2.close()

    if (retval != 0) throw RscriptException(retval)
  }

  /** Extracts a resource from the classpath and writes it to a temp file on disk. */
  private def writeResourceToTempFile(resource: String): Path = {
    val in = getClass.getClassLoader.getResourceAsStream(resource)
    if (in == null) throw new IllegalArgumentException(s"Resource '${resource} not found in classpath")
    val path = Io.makeTempFile("script.", ".R")
    path.toFile.deleteOnExit()
    val out = Io.toWriter(path)

    Source.fromInputStream(in).getLines().foreach { l => out.write(l); out.newLine() }
    in.close()
    out.close()
    path
  }
}
