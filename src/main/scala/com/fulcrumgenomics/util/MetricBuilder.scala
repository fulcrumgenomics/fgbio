/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.commons.CommonsDef.{forloop, unreachable}
import com.fulcrumgenomics.commons.reflect.{ReflectionUtil, ReflectiveBuilder}
import com.fulcrumgenomics.commons.util.{LazyLogging, Row}

import java.io.{PrintWriter, StringWriter}
import scala.reflect.runtime.{universe => ru}
import scala.util.{Failure, Success}

case class MetricBuilderException private[util](message: Option[String] = None) extends Exception

/** Class for building metrics of type [[T]].
  *
  * This is not thread-safe.
  *
  * @param source optionally, the source of reading (e.g. file)
  * @tparam T the metric type
  */
class MetricBuilder[T <: Metric](source: Option[String] = None)(implicit tt: ru.TypeTag[T]) extends LazyLogging {
  // The main reason why a builder is necessary is to cache some expensive reflective calls.
  private val clazz: Class[T]   = ReflectionUtil.typeTagToClass[T]
  private val reflectiveBuilder = new ReflectiveBuilder(clazz)
  private val names             = Metric.names[T]
  private val namesSet          = names.toSet

  /** Builds a metric from a [[ com.fulcrumgenomics.commons.util.Row]] produced from a
    * [[ com.fulcrumgenomics.commons.util.DelimitedDataParser]].
    *
    * @param row        the row to parse.  All required fields must be given.  Can be in any order.
    * @param lineNumber optionally, the line number when building a metric from a line in a file
    * @param ignoreExtra ignore extra fields (that are not in the metric)
    * @return
    */
  def fromRow(row: Row, headers: Seq[String], lineNumber: Option[Int] = None, ignoreExtra: Boolean = true): T = {
    val argMap = headers  // NB: filter after map so the row lookup is faster
      .zipWithIndex
      .map { case (header, i) => header -> row[String](i) }
      .filter { case (key, _) => !ignoreExtra || namesSet.contains(key) }
      .toMap
    fromArgMap(argMap=argMap, lineNumber=lineNumber)
  }

  /** Builds a metric from values for the complete set of metric fields
    *
    * @param values     the values in the same order as the names defined in the class
    * @param lineNumber optionally, the line number when building a metric from a line in a file
    * @return
    */
  def fromValues(values: Iterable[String], lineNumber: Option[Int] = None): T = {
    val vals = values.toIndexedSeq
    if (names.length != vals.length) {
      fail(message = f"Failed decoding: expected '${names.length}' fields, found '${vals.length}'.", lineNumber = lineNumber)
    }
    fromArgMap(argMap = names.zip(values).toMap, lineNumber = lineNumber)
  }

  /** Builds a metric of type [[T]]
    *
    * @param argMap     map of field names to values.  All required fields must be given.  Can be in any order.
    * @param lineNumber optionally, the line number when building a metric from a line in a file
    * @return a new instance of type [[T]]
    */
  private[util] def fromArgMap(argMap: Map[String, String], lineNumber: Option[Int] = None): T = {
    reflectiveBuilder.reset() // reset the arguments to their initial values

    val names = argMap.keys.toIndexedSeq
    forloop(from = 0, until = names.length) { i =>
      reflectiveBuilder.argumentLookup.forField(names(i)) match {
        case Some(arg) =>
          val value = {
            val tmp = argMap(names(i))
            if (tmp.isEmpty && arg.argumentType == classOf[Option[_]]) ReflectionUtil.SpecialEmptyOrNoneToken else tmp
          }

          val argumentValue = ReflectionUtil.constructFromString(arg.argumentType, arg.unitType, value) match {
            case Success(v) => v
            case Failure(thr) =>
              fail(
                message    = s"Could not construct value for column '${arg.name}' of type '${arg.typeDescription}' from '$value'",
                throwable  = Some(thr),
                lineNumber = lineNumber
              )
          }
          arg.value = argumentValue
        case None =>
          fail(
            message    = s"Did not have a field with name '${names(i)}'.",
            lineNumber = lineNumber
          )
      }
    }

    // build it.  NB: if arguments are missing values, then an exception will be thrown here
    // Also, we don't use the default "build()" method since if a collection or option is empty, it will be treated as
    // missing.
    val params = reflectiveBuilder.argumentLookup.ordered.map(arg => arg.value getOrElse buildFailure(s"Arguments not set: ${arg.name}"))
    try {
      reflectiveBuilder.build(params)
    } catch {
      case e: Exception => throw buildFailure(message=s"Failure building metric: ${e.getMessage}")
    }
  }

  /** Logs the throwable, if given, and builds a [[MetricBuilderException]] with information about when reading metrics
    *  fails
    *
    * @param message    the message to include in the exception thrown
    * @param throwable  optionally, a throwable that should be logged
    * @param lineNumber optionally, the line number when building a metric from a line in a file
    */
  private[util] def buildFailure(message: String, throwable: Option[Throwable] = None, lineNumber: Option[Int] = None): MetricBuilderException = {
    throwable.foreach { thr =>
      val stringWriter = new StringWriter
      thr.printStackTrace(new PrintWriter(stringWriter))
      val banner = "#" * 80
      logger.debug(banner)
      logger.debug(stringWriter.toString)
      logger.debug(banner)
    }
    val sourceMessage = source.map("\nIn source: " + _).getOrElse("")
    val prefix = lineNumber match {
      case None => "For metric"
      case Some(n) => s"On line #$n for metric"
    }
    val fullMessage = s"$prefix '${clazz.getSimpleName}'$sourceMessage\n$message"

    MetricBuilderException(message = Some(fullMessage))
  }

  /** Logs the throwable, if given, and throws a [[MetricBuilderException]] with information about when reading metrics fails
    *
    * @param message    the message to include in the exception thrown
    * @param throwable  optionally, a throwable that should be logged
    * @param lineNumber optionally, the line number when building a metric from a line in a file
    */
  private[util] def fail(message: String, throwable: Option[Throwable] = None, lineNumber: Option[Int] = None): Unit = {
    throw buildFailure(message=message, throwable=throwable, lineNumber=lineNumber)
  }
}
