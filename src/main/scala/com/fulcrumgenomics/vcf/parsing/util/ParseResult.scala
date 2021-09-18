/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.vcf.parsing.util

import com.fulcrumgenomics.commons.util.{LogLevel => FgLogLevel, Logger => FgLogger}
import com.fulcrumgenomics.sopt.util.{KERROR, KNRM}
import fastparse.Parsed.{Failure, Success}
import fastparse.internal.Util
import fastparse.{P, Parsed, ParserInput, parse}


/** Methods for dealing with parse failures. */
object ParseResult {

  /** Returns a parser that succeeds with the given value. */
  def success[T](value: T)(implicit ctx: P[_]): P[T] = ctx.freshSuccess(value=value)

  /** Returns a parser that fails at the given start index in the input and with the given message. */
  def fail(startIndex: Int = 0, message: String = "fail")(implicit ctx: P[_]): P[Nothing] = {
    val res: P[Nothing] = ctx.freshFailure()
    if (ctx.verboseFailures) {
      ctx.aggregateTerminal(
        startIndex = startIndex,
        f          = () => message
      )
    }
    res
  }

  /** Parses the input with the given parser and returns the value if successful, or throwing an exception upon failure. */
  def parseOrExcept[T](input: ParserInput, parser: P[_] => P[T]): T = {
    parse(input, parser(_)) match {
      case Success(value, _) => value
      case failure: Failure  => throw new Exception(formatFailure(failure=failure))
    }
  }

  implicit class LogFailure(logger: FgLogger) {
    /** Logs the failure to the given logger.
      *
      * @param failure the failure to log
      * @param lineNumber the line number for the failure
      * @param logBuffer the log buffer from fastparse
      * @return the logged message
      */
    def fail(failure: Parsed.Failure, lineNumber: Int, logBuffer: Iterator[String] = Iterator.empty): String = {
      val tracedFailure = failure.trace()
      val message = f"On line $lineNumber%,d: " + formatFailure(tracedFailure)
      logger.error(message)
      if (FgLogLevel.Debug.compareTo(FgLogger.level) >= 0) {
        val terminalStack = Failure.formatStack(tracedFailure.input, List(tracedFailure.terminalAggregateString -> tracedFailure.index))
        logger.debug(s"terminal: $terminalStack")
        val groupStack = Failure.formatStack(tracedFailure.input, List(tracedFailure.groupAggregateString -> tracedFailure.index))
        logger.debug(s"group: $groupStack")
        logger.debug(s"log:")
        logBuffer.foreach { logLine => logger.debug(s"\t$logLine") }
      }
      message
    }
  }

  /** Returns a failure message from the traced failure */
  def formatFailure(tracedFailure: Parsed.TracedFailure): String = {
    val Parsed.Failure(label, index, extra) = tracedFailure.failure
    val input = {
      val prefix = literalize(extra.input.slice(0, index))
      val suffix = literalize(extra.input.slice(index, extra.input.length))
      prefix + KERROR.code + "[" + suffix + "]" + KNRM.code
    }

    val stack = List(tracedFailure.terminalAggregateString -> index)
    Failure.formatStack(input, stack)

    label match {
      case "" =>
        val columnNumber = index + 1
        s"column `$columnNumber`, found `$input`"
      case s =>
        val lineNumberLookup = Util.lineNumberLookup(extra.input.toString)
        val lineNumber = lineNumberLookup.indexWhere(_ > index) match{
          case -1 => lineNumberLookup.length - 1
          case n => math.max(0, n - 1)
        }
        val columnNumber = index - lineNumberLookup(lineNumber) + 1
        s"column `$columnNumber`: expected `$s`, found `$input`"
    }
  }

  /** Returns a failure message from the given failure */
  def formatFailure(failure: Parsed.Failure): String = {
    formatFailure(failure.trace())
  }

  /** Creates a literal for the given character sequence. */
  def literalize(s: IndexedSeq[Char], unicode: Boolean = false): String = {
    Util.literalize(s=s, unicode=unicode).drop(1).dropRight(1)
  }
}
