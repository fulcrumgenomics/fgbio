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

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.IlluminaAdapters
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.util.SequenceUtil
import scala.collection.JavaConversions.asScalaIterator

class FindTechnicalReadsTest extends UnitSpec {
  "FindTechnicalReads" should "run and output the expected reads to a file." in {
    val builder  = new SamRecordSetBuilder(readLength=36)
    IlluminaAdapters.PairedEnd.both.flatMap(a => Seq(a, SequenceUtil.reverseComplement(a))).flatMap(_.sliding(36)).foreach(bases => {
      val Seq(r1, r2) = builder.addPair(start1=0, start2=0, record1Unmapped=true, record2Unmapped=true)
      r1.setReadString(bases)
      r2.setReadString("ACGA" * 9)
    })

    val expected = builder.map(_.getSAMString).toSeq

    builder.addPair(start1=300, start2=500).foreach(_.setReadString("A" * 36))
    builder.addPair(start1=300, start2=500).foreach(_.setReadString("C" * 36))
    builder.addPair(start1=300, start2=500).foreach(_.setReadString("G" * 36))
    builder.addPair(start1=300, start2=500).foreach(_.setReadString("T" * 36))

    val in  = builder.toTempFile()
    val out = makeTempFile("technical_sequences.", ".bam")
    new FindTechnicalReads(input=in, output=out).execute()
    val actual = SamReaderFactory.make().open(out.toFile).iterator().map(_.getSAMString).toSeq

    actual should contain theSameElementsAs expected
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
}
