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
package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.commons.CommonsDef.PathToFastq

object TrimFastqTest {
  private def fq(name:String, r1: String, r2:String): (FastqRecord, FastqRecord) = {
    (FastqRecord(name=name, bases=r1, quals=r1.replaceAll(".", "C"), comment=None, readNumber=None),
      FastqRecord(name=name, bases=r2, quals=r2.replaceAll(".", "C"), comment=None, readNumber=None))
  }

  // A set of test fastq records
  val FastqRecords : Seq[(FastqRecord, FastqRecord)] = Seq(
    fq("10x10", "ACGTACGTGT", "ACTGATCGAT"),
    fq("10x20", "ACGTACGTGT", "ACTGATCGATACTGATCGAT"),
    fq("20x20", "ACGTACGTGTACGTACGTGT", "ACTGATCGATACTGATCGAT")
  )
}

/**
  * Runs TrimFastq to test that it works!
  */
class TrimFastqTest extends UnitSpec {
  /** Writes the test records to a pair of files and returns them. */
  def fqFiles: (PathToFastq, PathToFastq) = {
    val (r1, r2) = (makeTempFile("r1.", ".fq"), makeTempFile("r2.", ".fq"))
    val (w1, w2) = (FastqWriter(r1), FastqWriter(r2))
    TrimFastqTest.FastqRecords.foreach(fqs => {
      w1.write(fqs._1)
      w2.write(fqs._2)
    })
    w1.close()
    w2.close()

    (r1, r2)
  }

  "TrimFastq" should "trim a single file and not discard any records" in {
    val (r1, r2) = fqFiles
    val out = makeTempFile("trimmed.", ".fq")
    new TrimFastq(input=Seq(r1), output=Seq(out), length=15, exclude=false).execute()
    val r1Map = FastqSource(out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 3
    r1Map("10x10").length shouldBe 10
    r1Map("10x20").length shouldBe 10
    r1Map("20x20").length shouldBe 15
  }

  it should "trim a single file and discard 2 records" in {
    val (r1, r2) = fqFiles
    val out = makeTempFile("trimmed.", ".fq")
    new TrimFastq(input=Seq(r1), output=Seq(out), length=15, exclude=true).execute()
    val r1Map = FastqSource(out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 1
    r1Map("20x20").length shouldBe 15
  }

  it should "trim a single file and discard 0 records because they are all long enough" in {
    val (r1, r2) = fqFiles
    val out = makeTempFile("trimmed.", ".fq")
    new TrimFastq(input=Seq(r1), output=Seq(out), length=5, exclude=true).execute()
    val r1Map = FastqSource(out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 3
    r1Map("10x10").length shouldBe 5
    r1Map("10x20").length shouldBe 5
    r1Map("20x20").length shouldBe 5
  }

  it should "not trim or discard any reads" in {
    val (r1, r2) = fqFiles
    val (r1Out, r2Out) = (makeTempFile("r1out.", ".fq"), makeTempFile("r2out.", ".fq"))
    new TrimFastq(input=Seq(r1, r2), output=Seq(r1Out, r2Out), length=25, exclude=false).execute()
    val r1Map = FastqSource(r1Out).map(r => r.name -> r).toMap
    val r2Map = FastqSource(r2Out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 3
    r2Map.size shouldBe r1Map.size
    r1Map("10x10").length shouldBe 10
    r1Map("10x20").length shouldBe 10
    r1Map("20x20").length shouldBe 20
    r2Map("10x10").length shouldBe 10
    r2Map("10x20").length shouldBe 20
    r2Map("20x20").length shouldBe 20
  }

  it should "trim but not discard some reads" in {
    val (r1, r2) = fqFiles
    val (r1Out, r2Out) = (makeTempFile("r1out.", ".fq"), makeTempFile("r2out.", ".fq"))
    new TrimFastq(input=Seq(r1, r2), output=Seq(r1Out, r2Out), length=15, exclude=false).execute()
    val r1Map = FastqSource(r1Out).map(r => r.name -> r).toMap
    val r2Map = FastqSource(r2Out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 3
    r2Map.size shouldBe r1Map.size
    r1Map("10x10").length shouldBe 10
    r1Map("10x20").length shouldBe 10
    r1Map("20x20").length shouldBe 15
    r2Map("10x10").length shouldBe 10
    r2Map("10x20").length shouldBe 15
    r2Map("20x20").length shouldBe 15
  }

  it should "trim some reads and discard others by pair in" in {
    val (r1, r2) = fqFiles
    val (r1Out, r2Out) = (makeTempFile("r1out.", ".fq"), makeTempFile("r2out.", ".fq"))
    new TrimFastq(input=Seq(r1, r2), output=Seq(r1Out, r2Out), length=15, exclude=true).execute()
    val r1Map = FastqSource(r1Out).map(r => r.name -> r).toMap
    val r2Map = FastqSource(r2Out).map(r => r.name -> r).toMap
    r1Map.size shouldBe 1
    r2Map.size shouldBe r1Map.size
    r1Map("20x20").length shouldBe 15
    r2Map("20x20").length shouldBe 15
  }
}
