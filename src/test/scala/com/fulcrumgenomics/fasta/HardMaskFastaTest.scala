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
package com.fulcrumgenomics.fasta

import java.nio.file.Files

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import dagr.commons.CommonsDef._
import dagr.commons.io.PathUtil
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

/**
  * Tests for HardMaskFasta
  */
class HardMaskFastaTest extends UnitSpec {
  val dir   = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/fasta")
  val fasta = dir.resolve("soft-masked.fa")

  "HardMaskFasta" should "convert soft masked bases to Ns" in {
    val out = makeTempFile("hard-masked.", ".fa")
    val LineLength: Int = 60
    val masker = new HardMaskFasta(input=fasta, output=out, lineLength = LineLength)
    masker.execute()

    val soft = ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta, false, false)
    val hard = ReferenceSequenceFileFactory.getReferenceSequenceFile(out, false, false)

    // Check over the sequence names and bases
    var (sOpt, hOpt) = (Option(soft.nextSequence()), Option(hard.nextSequence()))
    while (sOpt.isDefined || hOpt.isDefined) {
      (sOpt, hOpt) match {
        case (None, None) => unreachable("While loop checks that at least one is defined.")
        case (Some(x), None) => fail(s"Soft masked contained ${x.getName} after hard masked ran out of sequences.")
        case (None, Some(x)) => fail(s"Hard masked contained ${x.getName} after soft masked ran out of sequences.")
        case (Some(s), Some(h)) =>
          h.getName shouldBe s.getName
          h.getBases.length shouldBe s.getBases.length
          s.getBases.zip(h.getBases).foreach(pair => {
             if (pair._1.toChar.isLower) assert(pair._2 == 'N')
             else pair._2 shouldBe pair._1
          })
      }

      sOpt = Option(soft.nextSequence())
      hOpt = Option(hard.nextSequence())
    }

    // Check that the line length limit worked
    Io.readLines(out).foreach(line =>
      if (!line.startsWith(">")) line.length should be <= LineLength
    )

    Files.delete(out)
  }
}
