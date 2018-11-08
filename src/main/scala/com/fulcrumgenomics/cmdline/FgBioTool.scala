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
package com.fulcrumgenomics.cmdline

import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.{SAMFileHeader, SAMProgramRecord}

import scala.annotation.tailrec

/** Stores meta information about the command line use to invoke a tool.
  *
  * @param name the name of the tool.
  * @param args the list of arguments as given on the command line; this will contain any arguments also given to
  *             [[FgBioMain]] in [[FgBioCommonArgs]].
  * @param commandLineWithDefaults the command line as given to the tool, along with the defaults for all other
  *              arguments.  This will include arguments given to [[FgBioMain]] in [[FgBioCommonArgs]].
  * @param description the description of the tool
  * @param version the version of the tool.
  */
case class FgBioToolInfo(name: String, args: Seq[String], commandLineWithDefaults: String, description: String, version: String) extends Metric {

  /** The command line as given to the tool, without any defaults added.  This will include arguments given to
    * [[FgBioMain]] in [[FgBioCommonArgs]].*/
  def commandLineWithoutDefaults: String = args.mkString(" ")

  /** Adds a program group to the SAMFileHeader returning the ID. */
  def addProgramGroupTo(header: SAMFileHeader, id: Option[String]=None): SAMProgramRecord = {
    // Get the id
    val pgId = id.getOrElse {
      @tailrec
      def getPgId(intId: Int): String = if (header.getProgramRecord(intId.toString) == null) intId.toString else getPgId(intId + 1)
      getPgId(1)
    }
    val pg = new SAMProgramRecord(pgId)
    pg.setProgramName(name)
    pg.setCommandLine(commandLineWithDefaults)
    pg.setProgramVersion(version)
    header.addProgramRecord(pg)
    pg
  }
}


/** All fgbio tools should extend this. */
// Todo: extend Lazy Logging?
trait FgBioTool {
  /** Meta information about the command line use to invoke this tool, or [[None]] if unset. */
  private var _toolInfo: Option[FgBioToolInfo] = None

  /** All tools should implement this method. */
  def execute(): Unit

  /** Fail with just an exit code. */
  def fail(exit: Int) = throw new FailureException(exit=exit)

  /** Fail with the default exit code and a message. */
  def fail(message: String) = throw new FailureException(message=Some(message))

  /** Fail with a specific error code and message. */
  def fail(exit: Int, message: String) = throw new FailureException(exit=exit, message=Some(message))

  /** Generates a new validation exception with the given message. */
  def invalid(message: String) = throw new ValidationException(message)

  /** Generates a validation exception if the test value is false. */
  def validate(test: Boolean, message: => String) = if (!test) throw new ValidationException(message)

  /** Meta information about the command line use to invoke this tool, or [[None]] if unset. */
  def toolInfo: Option[FgBioToolInfo] = _toolInfo

  /** Sets the command line used to invoke this tool. */
  private[cmdline] def toolInfo_=(commandLine: FgBioToolInfo): Unit = this._toolInfo = Some(commandLine)
}
