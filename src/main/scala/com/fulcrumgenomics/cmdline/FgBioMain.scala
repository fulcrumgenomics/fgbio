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
    val startDate: Date = new Date

    // Print out some basic info
    logger.info("fgbio " + CommandLineProgramParserStrings.version(classOf[FgBioMain], color=false))
    logger.info("Executing as " + sysProp("user.name") + "@" + InetAddress.getLocalHost.getHostName +
      " on " + sysProp("os.name") + " " + sysProp("os.version") + " " + sysProp("os.arch"))
    logger.info(sysProp("java.vm.name") + " " + sysProp("java.runtime.version"))
    StringUtil.wordWrap(args.mkString(" "), 80).split("\n").filter(_.nonEmpty).zipWithIndex.foreach { case (str, idx) =>
      if (0 == idx) logger.info("Args: " + str)
      else logger.info("        " + str)
    }

    val parser = new CommandLineParser[FgBioTool]("fgbio")

    var name = this.getClass.getName
    val exitCode = parser.parseSubCommand(args=args, packageList=packageList) match {
      case None => 1
      case Some(tool) =>
        try {
          name = tool.getClass.getName
          logger.info(name + " starting " + name)
          tool.execute()
          0
        }
        catch {
          case ex: FailureException =>
            ex.message.foreach(logger.fatal)
            ex.exit
        }
    }

    val endDate: Date = new Date
    val elapsedMinutes: Double = (endDate.getTime - startDate.getTime) / (1000d * 60d)
    val elapsedString: String = new DecimalFormat("#,##0.00").format(elapsedMinutes)
    logger.info(s"$name done. Elapsed time: $elapsedString minutes.")
    logger.info("Runtime.totalMemory()=" + Runtime.getRuntime.totalMemory)

    exitCode
  }

  /** The packages we wish to include in our command line **/
  protected def packageList: List[String] = List[String]("com.fulcrumgenomics")

  private def sysProp(key: String): String = System.getProperty(key)
}
