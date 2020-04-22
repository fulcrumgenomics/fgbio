/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.testing.UnitSpec

class UpdateBedContigNamesTest extends UnitSpec {

  private val mappingLines = Seq(
    "NC_000001.10\tchr1",
    "NC_000002.11\tchr2",
    "NC_000003.11\tchr3",
    "NC_000004.11\tchr4"
  )

  private def mappingInput(skipLast: Boolean = false) = {
    val path = makeTempFile("test.", ".txt")
    if (skipLast) Io.writeLines(path, mappingLines.dropRight(1)) else Io.writeLines(path, mappingLines)
    path
  }

  private val bedLines =
    Seq(
      "NC_000001.10\t123\t500",
      "NC_000002.11\t444\t888\tother\tfields",
      "NC_000003.11\t8284\t10000\tmore",
      "NC_000004.11\t112345\t223456",
    )

  private def bed = {
    val out = makeTempFile("test.", ".bed")
    Io.writeLines(out, bedLines)
    out
  }

  "UpdateBedContigNames" should "update a BED" in {
    val output = makeTempFile("test.", ".bed")
    val tool = new UpdateBedContigNames(
      input   = bed,
      output  = output,
      mapping = mappingInput()
    )

    executeFgbioTool(tool)

    val lines = Io.readLines(output).toSeq
    lines should contain theSameElementsInOrderAs Seq(
      "chr1\t123\t500",
      "chr2\t444\t888\tother\tfields",
      "chr3\t8284\t10000\tmore",
      "chr4\t112345\t223456",
    )
  }

  it should "throw an exception if there are missing source contigs" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateBedContigNames(
      input   = bed,
      output  = output,
      mapping = mappingInput(skipLast = true)
    )

    val ex = intercept[Exception] {executeFgbioTool(tool) }
    ex.getMessage should include ("Did not find contig")
  }

  it should "skip missing source contigs when using --skip-missing" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateBedContigNames(
      input       = bed,
      output      = output,
      mapping     = mappingInput(skipLast = true),
      skipMissing = true
    )

    executeFgbioTool(tool)

    val lines = Io.readLines(output).toSeq
    lines should contain theSameElementsInOrderAs Seq(
      "chr1\t123\t500",
      "chr2\t444\t888\tother\tfields",
      "chr3\t8284\t10000\tmore",
    )
  }
}