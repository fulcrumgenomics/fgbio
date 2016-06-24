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
import dagr.commons.util.LazyLogging
import dagr.sopt.cmdline.CommandLineParser

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
    val parser = new CommandLineParser[FgBioTool]("fgbio")

    parser.parseSubCommand(args=args, packageList=packageList) match {
      case None => 1
      case Some(tool) =>
        try {
          tool.execute()
          0
        }
        catch {
          case ex: FailureException =>
            ex.message.foreach(logger.fatal)
            ex.exit
        }
    }
  }

  /** The packages we wish to include in our command line **/
  protected def packageList: List[String] = List[String]("com.fulcrumgenomics")
}
