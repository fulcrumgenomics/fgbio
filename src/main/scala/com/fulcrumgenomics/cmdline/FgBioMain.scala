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

import java.net.InetAddress
import java.text.DecimalFormat
import java.util.Date

import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.util.Io
import dagr.commons.util.{LazyLogging, StringUtil}
import dagr.sopt.cmdline.{CommandLineParser, CommandLineProgramParserStrings}

/**
  * Main program for fgbio that loads everything up and runs the appropriate sub-command
  */
object FgBioMain {
  /** The main method */
  def main(args: Array[String]): Unit = new FgBioMain().makeItSoAndExit(args)

  /**
    * Exception class intended to be used by [[FgBioMain]] and [[FgBioTool]] to communicate
    * non-exceptional(!) failures when running a tool.
    */
  case class FailureException private[cmdline] (exit:Int = 1, message:Option[String] = None) extends RuntimeException
}

class FgBioMain extends LazyLogging {
  /** A main method that invokes System.exit with the exit code. */
  def makeItSoAndExit(args: Array[String]): Unit = {
    System.exit(makeItSo(args))
  }

  /** A main method that returns an exit code instead of exiting. */
  def makeItSo(args: Array[String]): Int = {
    // Turn down HTSJDK logging
    htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)

    val startTime = System.currentTimeMillis()
    val parser = new CommandLineParser[FgBioTool](name)

    val exitCode = parser.parseSubCommand(args=args, packageList=packageList) match {
      case None => 1
      case Some(tool) =>
        val name = tool.getClass.getSimpleName
        try {
          printStartupLines(name, args)
          tool.execute()
          printEndingLines(startTime, name, true)
          0
        }
        catch {
          case ex: FailureException =>
            val banner = "#" * ex.message.map(_.length).getOrElse(80)
            logger.fatal(banner)
            logger.fatal("Execution failed!")
            ex.message.foreach(logger.fatal)
            logger.fatal(banner)
            printEndingLines(startTime, name, false)
            ex.exit
          case ex: Throwable =>
            printEndingLines(startTime, name, false)
            throw ex
        }
    }

    exitCode
  }

  protected def name: String = "fgbio"

  /** Prints a line of useful information when a tool starts executing. */
  protected def printStartupLines(tool: String, args: Array[String]): Unit = {
    val version    = CommandLineProgramParserStrings.version(getClass, color=false).replace("Version: ", "")
    val host       = InetAddress.getLocalHost.getHostName
    val user       = System.getProperty("user.name")
    val jreVersion = System.getProperty("java.runtime.version")
    logger.info(s"Executing $tool from $name version $version as $user@$host on JRE $jreVersion")
  }

  /** Prints a line of useful information when a tool stops executing. */
  protected def printEndingLines(startTime: Long, name: String, success: Boolean): Unit = {
    val elapsedMinutes: Double = (System.currentTimeMillis() - startTime) / (1000d * 60d)
    val elapsedString: String = new DecimalFormat("#,##0.00").format(elapsedMinutes)
    val verb = if (success) "completed" else "failed"
    logger.info(s"$name $verb. Elapsed time: $elapsedString minutes.")
  }

  /** The packages we wish to include in our command line **/
  protected def packageList: List[String] = List[String]("com.fulcrumgenomics")
}
