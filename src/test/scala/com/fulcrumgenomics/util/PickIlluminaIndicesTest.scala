/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

import java.nio.file.Path

import com.fulcrumgenomics.testing.UnitSpec

/**
  * Tests for PickIlluminaIndices.
  */
class PickIlluminaIndicesTest extends UnitSpec {

  def newOutput = makeTempFile("pick_illumina_indices_test.", ".txt")
  
  def parseLine(line: String): (String, Int) = {
    val data = line.split("\t").take(2)
    (data.head, data.last.toInt)
  }

  def checkOutput(output: Path, expectedNumIndices: Int, expectedEditDistance: Int): Unit = {
    val lines = Io.readLines(output).drop(1).toSeq
    lines.length shouldBe expectedNumIndices
    for (i <- lines.indices) {
      for (j <- i+1 until lines.length) {
        val (leftIndex, leftEditDistance) = parseLine(lines(i))
        val (rightIndex, rightEditDistance) = parseLine(lines(j))
        leftEditDistance should be >= expectedEditDistance
        rightEditDistance should be >= expectedEditDistance
        leftIndex.zip(rightIndex).count { case (l, r) => l != r } should be >= expectedEditDistance
      }
    }
  }

  "PickIlluminaIndices" should "pick indices with length four and edit distance one" in {
    Range(1, 5, 1).foreach { numIndices =>
      val output = newOutput
      new PickIlluminaIndices(
        length = 4,
        indices = numIndices,
        editDistance=1,
        output=output,
        allowReverses=true,
        allowReverseComplements=true,
        allowPalindromes=true,
        maxHomopolymer=8,
        threads=1
      ).execute()
      checkOutput(output=output, expectedNumIndices=numIndices, expectedEditDistance=1)
    }
  }

  it should "ask to pick four indices with length four and edit distance four but return only two" in {
    val output = newOutput
    new PickIlluminaIndices(
      length = 4,
      indices = 4,
      editDistance=4,
      output=output,
      allowReverses=true,
      allowReverseComplements=true,
      allowPalindromes=true,
      maxHomopolymer=8,
      threads=1
    ).execute()
    checkOutput(output=output, expectedNumIndices=2, expectedEditDistance=4)
  }

  it should "ask to pick three indices with length four and edit distance four from provided candidates" in {
    val candidates = {
      val path = makeTempFile("candidates", ".txt")
      Io.writeLines(path=path, lines=Seq("ATATATAT", "ATATATAA", "AAAATTTT", "GGGGGGGG"))
      path
    }
    val output = newOutput
    new PickIlluminaIndices(
      length = 8,
      indices = 3,
      editDistance=2,
      output=output,
      allowReverses=true,
      allowReverseComplements=true,
      allowPalindromes=true,
      minGc=0,
      maxGc=1,
      maxHomopolymer=8,
      threads=1,
      candidates=Some(candidates)
    ).execute()
    checkOutput(output=output, expectedNumIndices=3, expectedEditDistance=2)
    val indices = Io.readLines(output).drop(1).map(_.split('\t').head).toSeq
    indices should contain theSameElementsInOrderAs Seq("AAAATTTT", "ATATATAT", "GGGGGGGG")
  }

  // NB: this test fails due to a System.exit(1)
  /*
  it should "ask to pick four indices with length one and edit distance one but fail" in {
    val output = newOutput
    new PickIlluminaIndices(
      length = 1,
      indices = 4,
      editDistance=1,
      output=output,
      allowReverses=true,
      allowReverseComplements=true,
      allowPalindromes=true,
      maxHomopolymer=8,
      threads=1
    ).execute()
    Io.readLines(output).length shouldBe 3
  }
  */
}
