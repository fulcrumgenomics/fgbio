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
import com.fulcrumgenomics.commons.io.Io

object FastqIoTest {
  val someFastq =
    """
      |@Foo:279:000000000-ABCDE:1:1101:15033:1749 1:N:0:18
      |ATGACTGCGTATATATGACGTCGTGCTA
      |+
      |-,86,,;:C,<;-,86,,;:C,<;-,86
      |@Foo:279:000000000-ABCDE:1:1101:12430:1790 1:N:0:18
      |ATAGCTATGAGCTGCTGCTGA
      |+
      |8,,8,++++8788,,8,++++
      |@Foo:279:000000000-ABCDE:1:1101:9471:1291/1 1:N:0:18
      |CCCCCCCTCCTTGGGGGGGAGAGGG
      |+Foo:279:000000000-ABCDE:1:1101:9471:1291/1 1:N:0:18
      |AAA,,++78,66AAA,,++78,666
    """.stripMargin.trim().linesIterator.toSeq.dropWhile(_.isEmpty)

  val someFastqRecords = Seq(
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:15033:1749", "ATGACTGCGTATATATGACGTCGTGCTA", "-,86,,;:C,<;-,86,,;:C,<;-,86", Some("1:N:0:18"), None),
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:12430:1790", "ATAGCTATGAGCTGCTGCTGA",        "8,,8,++++8788,,8,++++",        Some("1:N:0:18"), None),
    FastqRecord("Foo:279:000000000-ABCDE:1:1101:9471:1291",  "CCCCCCCTCCTTGGGGGGGAGAGGG",    "AAA,,++78,66AAA,,++78,666",    Some("1:N:0:18"), Some(1))
  )
}

/**
  * Tests for the FastqSource
  */
class FastqIoTest extends UnitSpec {
  "FastqSource" should "read valid fastq from various kinds of input resources" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, FastqIoTest.someFastq)

    Seq(FastqSource(fq), FastqSource(fq.toFile), FastqSource(Io.toInputStream(fq)), FastqSource(FastqIoTest.someFastq)).foreach(source => {
      source.hasNext shouldBe true
      source.next() shouldBe FastqIoTest.someFastqRecords(0)
      source.hasNext shouldBe true
      source.next() shouldBe FastqIoTest.someFastqRecords(1)
      source.hasNext shouldBe true
      source.next() shouldBe FastqIoTest.someFastqRecords(2)
      source.hasNext shouldBe false
      an[NoSuchElementException] shouldBe thrownBy { source.next() }
    })
  }

  it should "throw an exception if a header line doesn't start with an @" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("foo", "ACGT", "+", "7777"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if a quality header line doesn't start with an +" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "@", "7777"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if there are not n*4 lines" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "+foo"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "throw an exception if seqs and quals are not the same length for a record" in {
    val fq = makeTempFile("some", ".fq")
    Io.writeLines(fq, Seq("@foo", "ACGT", "+foo", "12345"))
    an[IllegalStateException] shouldBe thrownBy { FastqSource(fq) }
  }

  it should "handle fastq that has read numbers greater than 2" in {
    val fq = makeTempFile("some.", ".fq")
    val out = FastqWriter(fq)
    Range.inclusive(1, 9).foreach { n =>
      out += FastqRecord(name=s"q$n", bases="ACGT", quals="####", comment=Some("Y"), readNumber=Some(n))
    }
    out.close()

    FastqSource(fq).zipWithIndex.foreach { case (rec, index) =>
      rec.readNumber shouldBe Some(index+1)
      rec.comment shouldBe Some("Y")
      rec.name.contains("/") shouldBe false
    }
  }

  "FastqWriter" should "write valid records to a file" in {
    val files = Seq(makeTempFile("some", ".fq"), makeTempFile("some", ".fq.gz"))
    files.foreach(fq => {
      val writer = FastqWriter(fq)
      FastqIoTest.someFastqRecords.foreach(writer.write)
      writer.close()
      FastqSource(fq).toSeq shouldBe FastqIoTest.someFastqRecords
    })
  }
}
