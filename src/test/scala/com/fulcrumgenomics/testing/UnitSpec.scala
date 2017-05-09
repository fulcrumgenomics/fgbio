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
package com.fulcrumgenomics.testing

import java.io.PrintStream
import java.nio.file.{Files, Path}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.commons.reflect.ReflectionUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, LogLevel, Logger}
import com.fulcrumgenomics.sopt.cmdline.CommandLineProgramParser
import com.fulcrumgenomics.sopt.util.ParsingUtil
import htsjdk.samtools.{SAMRecord, SamReaderFactory}
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

import scala.reflect.ClassTag
import scala.reflect.runtime.universe._

/** Base class for unit and integration testing */
trait UnitSpec extends FlatSpec with Matchers {
  // Turn down HTSJDK logging
  htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)

  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

  /** Reads all the records from a SAM or BAM file into an indexed seq. */
  protected def readBamRecs(bam: PathToBam): IndexedSeq[SAMRecord] = {
    val in = SamReaderFactory.make().open(bam.toFile)
    yieldAndThen(in.toIndexedSeq) { in.safelyClose() }
  }

  /** Reads all the records from a VCF file into an indexed seq. */
  protected def readVcfRecs(vcf: PathToVcf): IndexedSeq[VariantContext] = {
    val in = new VCFFileReader(vcf.toFile, false)
    yieldAndThen(in.toIndexedSeq) { in.safelyClose() }
  }

  /** Generates a command line parser for a class to check that the argument annotations are valid. */
  protected def checkClpAnnotations[T <: FgBioTool](implicit ct: ClassTag[T], tt: TypeTag[T]): Unit = {
    val cl   = ReflectionUtil.typeTagToClass[T]
    val name = cl.getName

    ParsingUtil.findClpAnnotation(cl).getOrElse(fail(s"${name} is missing the clp annotation."))
    new CommandLineProgramParser(cl)
  }

  /**
    * Executes the provided tool and returns the tools logging output as a list of String.
    */
  protected def executeFgbioTool(tool: FgBioTool with LazyLogging): Seq[String] = {
    val log = makeTempFile(tool.getClass.getSimpleName + ".", ".log")
    val stream = new PrintStream(log.toFile)

    // This is a little icky, but works without having to increase the visibility of 'logger' in LazyLogging
    val loggerAccessor = classOf[LazyLogging].getMethod("logger")
    val logger         = loggerAccessor.invoke(tool).asInstanceOf[Logger]

    val previousStream = logger.out
    logger.out = Some(stream)
    try {
      tool.execute()
    }
    finally {
      stream.close()
      logger.out = previousStream
    }

    Io.readLines(log).toSeq
  }
}

/** Base class that turns up logging to [[LogLevel.Error]] before all the tests and restores
  * the log level after all the tests.
  */
trait ErrorLogLevel extends UnitSpec with BeforeAndAfterAll {
  private var logLevel = Logger.level

  override protected def beforeAll(): Unit = {
    this.logLevel = Logger.level
    Logger.level  = LogLevel.Error
  }

  override protected def afterAll(): Unit = {
    Logger.level = LogLevel.Info
    Logger.level = this.logLevel
  }
}
