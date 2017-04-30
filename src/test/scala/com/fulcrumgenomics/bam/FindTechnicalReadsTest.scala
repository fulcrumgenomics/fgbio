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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.IlluminaAdapters
import htsjdk.samtools.util.SequenceUtil

class FindTechnicalReadsTest extends UnitSpec {
  "FindTechnicalReads" should "run and output the expected reads to a file." in {
    val builder  = new SamBuilder(readLength=36)
    IlluminaAdapters.PairedEnd.both.flatMap(a => Seq(a, SequenceUtil.reverseComplement(a))).flatMap(_.sliding(36)).foreach { bases =>
      builder.addPair(start1=0, start2=0, unmapped1=true, unmapped2=true, bases1=bases, bases2="ACGA" * 9)
    }

    val expected = builder.map(_.id).toSeq

    builder.addPair(start1=300, start2=500, bases1="A" * 36, bases2="A" * 36)
    builder.addPair(start1=300, start2=500, bases1="C" * 36, bases2="C" * 36)
    builder.addPair(start1=300, start2=500, bases1="G" * 36, bases2="G" * 36)
    builder.addPair(start1=300, start2=500, bases1="T" * 36, bases2="T" * 36)

    val in  = builder.toTempFile()
    val out = makeTempFile("technical_sequences.", ".bam")
    new FindTechnicalReads(input=in, output=out).execute()
    val actual = readBamRecs(out).map(_.id)

    actual should contain theSameElementsAs expected
  }

  "FindTechnicalReads" should "output all reads, tagging only those that are technical" in {
    val tseqs = IlluminaAdapters.DualIndexed.both ++ IlluminaAdapters.NexteraV2.both
    val builder  = new SamBuilder(readLength=30)
    builder.addFrag(name="ok1", start=100,     bases="AACCGGTTAACCGGTTAACCGGTTAACCGG")
    builder.addFrag(name="di1", unmapped=true, bases="CTCTTTCCCTACACGACGCTCTTCCGATCT")
    builder.addFrag(name="di2", unmapped=true, bases="AATGATACGGCGACCACCGAGATCTNNNNN")
    builder.addFrag(name="nx1", unmapped=true, bases="TCGTCGGCAGCGTCAGATGTGTATAAGAGA")

    val in  = builder.toTempFile()
    val out = makeTempFile("technical_sequences.", ".bam")
    new FindTechnicalReads(input=in, output=out, sequences=tseqs, allReads=true, tag=Some("XS")).execute()
    val actual = readBamRecs(out).map(r => r.name -> r).toMap

    actual("di1")[Int]("XS") shouldBe 0
    actual("di2")[Int]("XS") shouldBe 0
    actual("nx1")[Int]("XS") shouldBe 2
    actual should contain key "ok1"
    actual("ok1").get[Int]("XS") shouldBe None
  }

  "Matcher" should "match reads that come from the adapters and don't contain Ns" in {
    val matcher = new Matcher(matchLength=20, maxErrors=1, sequences=IlluminaAdapters.DualIndexed.both)
    IlluminaAdapters.DualIndexed.both.flatMap(_.sliding(36)).foreach(s => {
      val matched = matcher.matches(s.getBytes)
      if (s.substring(0, 20).contains('N'))
        assert(!matched, s"${s} should not have matched.")
      else
        assert(matched, s"${s} should have matched.")
    })
  }

  it should "match reads with errors in them" in {
    val matcher = new Matcher(matchLength=20, maxErrors=2, sequences=IlluminaAdapters.PairedEnd.both)
    IlluminaAdapters.PairedEnd.both.flatMap(_.sliding(36)).foreach(s => {
      val bs = s.getBytes
      assert(matcher.matches(bs), s"Should have matched unedited sequence ${s}")
      bs(3) = 'N'
      assert(matcher.matches(bs), s"Should have matched sequence with one errors ${new String(bs)}")
      bs(7) = 'N'
      assert(matcher.matches(bs), s"Should have matched sequence with two errors ${new String(bs)}")
      bs(13) = 'N'
      assert(!matcher.matches(bs), s"Should not have matched sequence with three errors ${new String(bs)}")
    })
  }

  it should "find the earliest sequence in the set that matches" in {
    val matcher = new Matcher(matchLength=4, maxErrors=0, sequences=Seq("ACGTACGT", "ACACGTAC", "ACACTGAT", "ACGTAAAA"))
    matcher.findMatch("ACGT".getBytes()) shouldBe Some("ACGTACGT")
    matcher.findMatch("ACAC".getBytes()) shouldBe Some("ACACGTAC")
    matcher.findMatch("AAAA".getBytes()) shouldBe Some("ACGTAAAA")
    matcher.findMatch("ACTG".getBytes()) shouldBe Some("ACACTGAT")
    matcher.findMatch("GGGG".getBytes()) shouldBe None
  }

}
