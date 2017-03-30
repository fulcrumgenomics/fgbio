/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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
import com.fulcrumgenomics.FgBioDef._

class MatrixTest extends UnitSpec {
  "Matrix.apply" should "build a matrix of appropriate size" in {
    val m = Matrix[Int](3, 2)
    m.dimensions shouldBe (3,2)
    m.size       shouldBe 6
  }

  "Matrix" should "remember values placed in the matrix" in {
    val m = Matrix[Int](2,10)
    forloop(from=0, until=2) { i =>
      forloop(from=0, until=10) { j =>
        m(i, j) = i * j
      }
    }

    forloop(from=0, until=2) { i =>
      forloop(from=0, until=10) { j =>
        m(i, j) shouldBe i * j
      }
    }
  }

  it should "throw exceptions for invalid indices" in {
    val m = Matrix[Int](4,5)
    an[Exception] shouldBe thrownBy { m(-1, -1) }
    an[Exception] shouldBe thrownBy { m( 1,  7) }
    an[Exception] shouldBe thrownBy { m( 6,  1) }
    an[Exception] shouldBe thrownBy { m( 5,  5) }

    an[Exception] shouldBe thrownBy { m(-1, -1) = 1}
    an[Exception] shouldBe thrownBy { m( 3, -1) = 1}
    an[Exception] shouldBe thrownBy { m(-1,  3) = 1}
    an[Exception] shouldBe thrownBy { m( 1,  7) = 1}
    an[Exception] shouldBe thrownBy { m( 6,  1) = 1}
    an[Exception] shouldBe thrownBy { m( 5,  5) = 1}
  }

  it should "display a matrix in a reasonable format for debugging" in {
    val m = Matrix[Int](3,2)
    forloop(from=0, until=3) { i =>
      forloop(from=0, until=2) { j =>
        m(i, j) = i * j
      }
    }

    val table = m.toString
    val expected =
      """
        |0  1
        |0  2
        |0  3
      """.stripMargin.trim.replaceAll(" +", "\t")
  }
}
