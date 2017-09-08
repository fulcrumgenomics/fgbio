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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.testing.UnitSpec

class SimpleConsensusCallerTest extends UnitSpec {

  "SimpleConsensusCaller" should "fail if no sequences are given" in {
    val caller = new SimpleConsensusCaller()
    an[Exception] should be thrownBy caller.callConsensus(Seq.empty)
  }

  it should "fail if the sequences have different lengths" in {
    val caller = new SimpleConsensusCaller()
    an[Exception] should be thrownBy caller.callConsensus(Seq("A", "AC"))
  }

  it should "fail if not all strings had the same non-ACGTN character at the same position" in {
    val caller = new SimpleConsensusCaller()
    an[Exception] should be thrownBy caller.callConsensus(Seq("GATT-ACA", "GATT-ACA", "GATTAACA"))
  }

  it should "create a consensus from the sequences that all agree" in {
    val caller = new SimpleConsensusCaller()
    caller.callConsensus(Seq("A", "A")) shouldBe "A"
    caller.callConsensus(Seq("GATTACA", "GATTACA")) shouldBe "GATTACA"
  }

  it should "create a consensus from sequences that differ" in {
    val caller = new SimpleConsensusCaller()
    caller.callConsensus(Seq("A", "C", "G", "T")) shouldBe "N"
    caller.callConsensus(Seq("A", "C", "C", "C")) shouldBe "C"
    caller.callConsensus(Seq("C", "C", "C", "A")) shouldBe "C"
    caller.callConsensus(Seq("GATTACA", "GATTACA", "GATTACA", "NNNNNNN")) shouldBe "GATTACA"
  }

  it should "gracefully handle non-ACGTN bases" in {
    val caller = new SimpleConsensusCaller()
    caller.callConsensus(Seq("GATT-ACA", "GATT-ACA", "GATT-ACA")) shouldBe "GATT-ACA"
    caller.callConsensus(Seq("XGAT", "XGAT")) shouldBe "XGAT"
    caller.callConsensus(Seq("GATY", "GATY")) shouldBe "GATY"
  }
}
