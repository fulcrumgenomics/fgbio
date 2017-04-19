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
package com.fulcrumgenomics.util

import java.nio.file.Paths

import com.fulcrumgenomics.testing.UnitSpec

import scala.collection.mutable

/**
  * Tests for the AsyncStreamSink class
  */
class AsyncStreamSinkTest extends UnitSpec {
  "AsyncStreamSink" should "capture all the output" in {
    val file = Paths.get("src/test/resources/com/fulcrumgenomics/util/async-stream-sink-test.txt")
    val expected = Io.toSource(file).getLines().toList

    val proc = new ProcessBuilder("cat", file.toString).start()
    val actual = mutable.ListBuffer[String]()
    val sink = new AsyncStreamSink(proc.getInputStream, s => actual.append(s))
    val exit = proc.waitFor()
    sink.close(1000)

    exit shouldBe 0
    actual.toList shouldBe expected
  }
}
