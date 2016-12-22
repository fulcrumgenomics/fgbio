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
 *
 */

package com.fulcrumgenomics.util

import java.io.{BufferedWriter, Closeable, Writer}
import java.nio.file.Path

import dagr.commons.CommonsDef.unreachable
import dagr.commons.reflect.{ReflectionUtil, ReflectiveBuilder}
import htsjdk.samtools.util.FormatUtil

import scala.reflect.runtime.{universe => ru}
import scala.util.{Failure, Success}

object Metric {
  val Delimiter: Char = '\t'
  val DelimiterAsString: String = s"$Delimiter"

  /** A class that provides streaming writing capability for metrics. */
  class MetricWriter[T <: Metric] private[Metric](val writer: Writer)(implicit tt: ru.TypeTag[T]) extends Closeable {
    this.writer.write(names[T].mkString(DelimiterAsString))
    this.writer.write("\n")

    /** Writes one or more metric values to the output. */
    def write(metrics: T*): Unit = metrics.foreach { metric =>
      writer.write(metric.values.mkString(DelimiterAsString))
      writer.write("\n")
    }

    def flush(): Unit = this.writer.flush()
    override def close(): Unit = this.writer.close()
  }

  /** Get the names of the arguments in the order they were defined for the type [T]. */
  def names[T <: Metric](implicit tt: ru.TypeTag[T]): Seq[String] = {
    val clazz             = ReflectionUtil.typeTagToClass[T]
    val reflectiveBuilder = new ReflectiveBuilder(clazz)
    reflectiveBuilder.argumentLookup.ordered.map(_.name)
  }

  /** Construct a metric of the given type.  The arguments given should be a sequence of tuples: the name of the field
    * and value as a string. */
  def build[T <: Metric](args: Seq[(String, String)])(implicit tt: ru.TypeTag[T]): T = {
    val clazz             = ReflectionUtil.typeTagToClass[T]
    val reflectiveBuilder = new ReflectiveBuilder(clazz)
    // set the arguments
    args.foreach { case (name, value) =>
      reflectiveBuilder.argumentLookup.forField(name) match {
        case Some(arg) =>
          val argumentValue = ReflectionUtil.constructFromString(arg.argumentType, arg.unitType, value) match {
            case Success(v) => v
            case Failure(thr) => throw thr
          }
          arg.value_=(argumentValue)
        case None =>
          throw new IllegalArgumentException(s"The class '${clazz.getSimpleName}' did not have a field with name '$name'.")
      }
    }
    // build it.  NB: if arguments are missing values, then an exception will be thrown here
    // Also, we don't use the default "build()" method since if a collection or option is empty, it will be treated as
    // missing.
    reflectiveBuilder.build(reflectiveBuilder.argumentLookup.ordered.map(arg => arg.value getOrElse unreachable(s"Arguments not set: ${arg.name}")))
  }

  /** Reads metrics from a set of lines.  The first line should be the header with the field names.  Each subsequent
    * line should be a single metric. */
  def read[T <: Metric](lines: Iterator[String], source: Option[String] = None)(implicit tt: ru.TypeTag[T]): Seq[T] = {
    if (lines.isEmpty) throw new IllegalArgumentException(s"No header found in metrics" + source.map(" in source: " + _).getOrElse(""))
    val parser = new DelimitedDataParser(lines=lines, delimiter=Delimiter, ignoreBlankLines=false, trimFields=true)
    val names  = parser.headers
    parser.zipWithIndex.map { case (row, lineNumber) =>
      val values = names.indices.map { i => row[String](i) }.map { value =>
        if (value.nonEmpty) value else ReflectionUtil.SpecialEmptyOrNoneToken
      }
      build[T](names.zip(values))
    }.toSeq
  }

  /** Reads metrics from the given path.  The first line should be the header with the field names.  Each subsequent
    * line should be a single metric. */
  def read[T <: Metric](path: Path)(implicit tt: ru.TypeTag[T]): Seq[T] = read[T](Io.readLines(path), Some(path.toString))

  /** Writes a metric to the given writer.  The first line will be a header with the field names.  Each subsequent
    * line is a single metric.
    */
  def write[T <: Metric](metrics: Seq[T], writer: Writer)(implicit tt: ru.TypeTag[T]): Unit = {
    val out = new MetricWriter[T](writer)
    out.write(metrics:_*)
    out.close()
  }

  /** Writes a metric to the given path.  The first line will be a header with the field names.  Each subsequent
    * line is a single metric.
    */
  def write[T <: Metric](metrics: Seq[T], path: Path)(implicit tt: ru.TypeTag[T]): Unit = write(metrics, Io.toWriter(path))

  /** Returns a MetricWriter that can be used to stream metrics out to a file. */
  def writer[T <: Metric](path: Path)(implicit tt: ru.TypeTag[T]): MetricWriter[T] = new MetricWriter[T](Io.toWriter(path))
}

/**
  * Base trait for metrics.
  *
  * All classes extending this class should be a case class.  By convention, all fields should be lower case with
  * words separated by underscores.
  */
trait Metric extends Product {
  private lazy val reflectiveBuilder = new ReflectiveBuilder(this.getClass)
  private val formatter = new FormatUtil {
    override def format(value: Any): String = value match {
      case None              => ""
      case Some(x)           => super.format(x)
      case y                 => super.format(y)
    }
  }

  /** Get the names of the arguments in the order they were defined. */
  def names: Seq[String] = {
    this.reflectiveBuilder.argumentLookup.ordered.map(_.name)
  }

  /** Get the values of the arguments in the order they were defined. */
  def values: Seq[String] = productIterator.map(formatValues).toSeq

  /** Override this method to customize how values are formatted. */
  protected def formatValues(value: Any): String = formatter.format(value)

  override def toString: String = names.zip(values).toMap.toString
}
