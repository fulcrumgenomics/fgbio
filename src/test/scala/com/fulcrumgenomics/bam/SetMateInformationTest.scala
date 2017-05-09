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
package com.fulcrumgenomics.bam

import java.nio.file.Paths

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import htsjdk.samtools.SamReaderFactory

/**
  * Tests running SetMateInformation on a few different input files. Does not exhaustively
  * verify the mate information that was set because that is all done using HTSJDK
  * functionality that is tested there.
  */
class SetMateInformationTest extends UnitSpec {
  val dir = Paths.get("src/test/resources/com/fulcrumgenomics/bam")
  val querySortedSam      = dir.resolve("set_mate_querysorted.sam")
  val queryGroupedSam     = dir.resolve("set_mate_querygrouped.sam")
  val coordinateSortedSam = dir.resolve("set_mate_coordinatesorted.sam")

  "SetMateInformation" should "correctly set mate information in a query sorted file" in {
    val out = makeTempFile("mated.", ".bam")
    val fixer = new SetMateInformation(input=querySortedSam, output=out)
    fixer.execute()
    val in = SamReaderFactory.make().open(out.toFile)
    in.iterator().filterNot(r => r.getReadUnmappedFlag || r.getMateUnmappedFlag).foreach(rec => {
      rec.getAttribute("MC") should not be null
      rec.getAttribute("MQ") should not be null
    })
    in.close()
  }

  it should "correctly set mate information in a query GROUPED file" in {
    val out = makeTempFile("mated.", ".bam")
    val fixer = new SetMateInformation(input=queryGroupedSam, output=out)
    fixer.execute()
    val in = SamReaderFactory.make().open(out.toFile)
    in.iterator().filterNot(r => r.getReadUnmappedFlag || r.getMateUnmappedFlag).foreach(rec => {
      rec.getAttribute("MC") should not be null
      rec.getAttribute("MQ") should not be null
    })
    in.close()
  }

  it should "throw an exception on a coordinate sorted file" in {
    val out = makeTempFile("mated.", ".bam")
    an[ValidationException] should be thrownBy new SetMateInformation(input=coordinateSortedSam, output=out)
  }
}
