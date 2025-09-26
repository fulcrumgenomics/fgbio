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

import com.fulcrumgenomics.commons.CommonsDef.DirPath

import scala.reflect.runtime.{universe => ru}

/** Disk-backed metrics sorter
  *
  * @param maxObjectsInRam the maximum number of metrics to keep in memory before spilling to disk
  * @param keyfunc method to convert a metric to an ordered key
  * @param tmpDir the temporary directory in which to spill to disk
  * @param tt the type tag for `T`
  * @tparam Key the key to use for sorting metrics
  * @tparam T the metric type
  */
class MetricSorter[Key <: Ordered[Key], T <: Metric](maxObjectsInRam: Int = MetricSorter.MaxInMemory,
                                                     keyfunc: T => Key,
                                                     tmpDir: DirPath = Io.tmpDir,

                                                    )(implicit tt: ru.TypeTag[T]) extends Sorter[T, Key](
  maxObjectsInRam = maxObjectsInRam,
  codec           = new MetricSorter.MetricSorterCodec[T](),
  keyfunc         = keyfunc,
  tmpDir          = tmpDir
)

object MetricSorter {
  /** The default maximum # of records to keep and sort in memory. */
  val MaxInMemory: Int = 1e6.toInt

  /** The codec for encoding and decoding a metric */
  class MetricSorterCodec[T <: Metric]()(implicit tt: ru.TypeTag[T])
    extends Sorter.Codec[T] {
    private val builder = new MetricBuilder[T]()

    /** Encode the metric into an array of bytes. */
    def encode(metric: T): Array[Byte] = metric.values.mkString(Metric.DelimiterAsString).getBytes

    /** Decode a metric from an array of bytes. */
    def decode(bs: Array[Byte], start: Int, length: Int): T = {
      val fields = new String(bs.slice(from = start, until = start + length)).split(Metric.DelimiterAsString)
      builder.fromValues(fields)
    }
  }
}