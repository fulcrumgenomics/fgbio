/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.FgBioDef.PathToFasta

class UpdateFastaContigNamesTest extends UnitSpec {

  private val dir   = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/fasta")
  private val fasta = dir.resolve("soft-masked.fa")


  private def mappingInput(skipLast: Boolean = false) = {
    val lines = Seq(
      "gi|7|emb|X51700.1|\t1",
      "gi|20|emb|X52703.1|\t2",
      "gi|595576|gb|U13080.1|BPU13080\t3"
    )
    val path = makeTempFile("test.", ".txt")
    if (skipLast) Io.writeLines(path, lines.dropRight(1)) else Io.writeLines(path, lines)
    path
  }

  private def getNames(path: PathToFasta): Seq[String] = {
    ReferenceSequenceIterator(path, stripComments = true).map(_.getName).toSeq
  }

  "UpdateFastaContigNames" should "update a simple FASTA" in {
    val output = makeTempFile("test.", ".fasta")

    val tool = new UpdateFastaContigNames(
      input   = fasta,
      output  = output,
      mapping = mappingInput()
    )

    executeFgbioTool(tool)

    getNames(output) should contain theSameElementsInOrderAs Seq("1", "2", "3")
  }

  it should "throw an exception if there are missing source contigs" in {
    val output = makeTempFile("test.", ".fasta")

    val tool = new UpdateFastaContigNames(
      input   = fasta,
      output  = output,
      mapping = mappingInput(skipLast = true)
    )
    val ex = intercept[Exception] {executeFgbioTool(tool) }
    ex.getMessage should include ("Did not find contig")
  }

  it should "skip missing source contigs when using --skip-missing" in {
    val output = makeTempFile("test.", ".fasta")

    val tool = new UpdateFastaContigNames(
      input       = fasta,
      output      = output,
      mapping     = mappingInput(skipLast = true),
      skipMissing = true
    )
    executeFgbioTool(tool)

    getNames(output) should contain theSameElementsInOrderAs Seq("1", "2")
  }
}
