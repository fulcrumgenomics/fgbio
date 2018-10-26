/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.testing.UnitSpec
import org.scalatest.concurrent.ScalaFutures

class AsyncAlignerTest extends UnitSpec with ScalaFutures {
  import scala.concurrent.ExecutionContext.Implicits.global

  private val testTasks = Seq(("ACGTAACC", "ACGTAACC"), ("CCGG", "AACCGGTT"), ("AAAATTTT", "AAAATTTTGGGG"))

  private def newAligner(mode: Mode) = Aligner(1, -4, -6, -1, mode=mode)

  Seq(Mode.Local, Mode.Glocal, Mode.Global).foreach { mode =>
    "AsyncAligner" should s"align in $mode mode" in {
      val asyncAligner = new AsyncScalaAligner(newAligner(mode))
      val aligner      = newAligner(mode)
      val futures      = testTasks.map { case (query, target) =>
          asyncAligner.align(query, target)
      }
      futures.length shouldBe testTasks.length
      testTasks.zip(futures).map { case ((query, target), future) =>
        whenReady(future) { actual =>
          actual shouldBe aligner.align(query, target)
        }
      }
    }
  }
}
