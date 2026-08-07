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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

import java.nio.file.Paths

class CallDuplexConsensusReadsTest extends UnitSpec {
  private val MI = ConsensusTags.MolecularId
  private val RX = ConsensusTags.UmiBases

  "CallDuplexConsensusReads" should "throw an exception if the input file doesn't exist" in {
    an[Throwable] should be thrownBy {
      new CallDuplexConsensusReads(input=Paths.get("/tmp/path/to/no/where/foo.bam"), output=Paths.get("/tmp")).execute()
    }
  }

  it should "throw an exception if the output file isn't writable" in {
    an[Throwable] should be thrownBy {
      val in = makeTempFile("in.", ".bam")
      val out = Paths.get("/tmp/path/to/no/where.bam")
      new CallDuplexConsensusReads(input=in, output=out).execute()
    }
  }

  it should "throw an exception if either error rate is set too low" in {
      val in  = makeTempFile("in.", ".bam")
      val out = makeTempFile("out.", ".bam")
    an[Exception] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, errorRatePreUmi=0.toByte).execute() }
    an[Exception] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, errorRatePostUmi=0.toByte).execute() }
  }

  it should "throw a validation exception if --max-reads-per-strand is less than one" in {
    val in  = makeTempFile("in.", ".bam")
    val out = makeTempFile("out.", ".bam")
    // Validation happens at construction, so this fails only when the validation exists - not because the
    // empty input file is unreadable.
    an[ValidationException] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, maxReadsPerStrand=Some(0)) }
    an[ValidationException] should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, maxReadsPerStrand=Some(-1)) }
    noException should be thrownBy { new CallDuplexConsensusReads(input=in, output=out, maxReadsPerStrand=Some(1)) }
  }

  it should "have working CLP and arg annotations" in {
    checkClpAnnotations[CallDuplexConsensusReads]
  }

  it should "not generate a consensus if AB-R1s are not on the same strand ads BA-R2s" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="ab1", start1=100, start2=200, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Plus)
    builder.addPair(name="ba1", start1=200, start2=100, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Minus)

    val in  = builder.toTempFile()
    val out = makeTempFile("duplex.", ".bam")
    new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq

    reader.header.getReadGroups should have size 1
    reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
    recs should have size 0
  }

  it should "not generate a consensus if AB-R2s are not on the same strand ads BA-R1s" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="ab1", start1=100, start2=200, attrs=Map(MI -> "1/A"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Minus)
    builder.addPair(name="ba1", start1=200, start2=100, attrs=Map(MI -> "1/B"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA", strand1=Plus, strand2=Plus)

    val in  = builder.toTempFile()
    val out = makeTempFile("duplex.", ".bam")
    new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ").execute()
    val reader = SamSource(out)
    val recs = reader.toSeq

    reader.header.getReadGroups should have size 1
    reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
    recs should have size 0
  }

  Seq(1, 2, 4).foreach { threads =>
    it should s"run successfully and create consensus reads with $threads threads" in {
      val specialCellTag = "XX"
      val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
      builder.addPair(name="ab1", start1=100, start2=100, attrs=Map(MI -> "1/A", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
      builder.addPair(name="ab2", start1=100, start2=100, attrs=Map(MI -> "1/A", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
      builder.addPair(name="ab3", start1=100, start2=100, attrs=Map(MI -> "1/A", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
      builder.addPair(name="ba1", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
      builder.addPair(name="ba2", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")
      builder.addPair(name="ba3", start1=100, start2=100, strand1=Minus, strand2=Plus, attrs=Map(MI -> "1/B", "XX" -> "AB"), bases1="AAAAAAAAAA", bases2="AAAAAAAAAA")

      // Add the original UMI bases to each read
      builder.foreach { rec =>
        val mi = rec[String](MI)
        // The UMIs for ABR1s+BAR2s are called together, and the UMIs for ABR2s and BAR1s are called together.  But
        // before they are called together, they BAR1s and ABR2s have their UMI swapped, so we have to swap them here
        // too.
        (rec.firstOfPair, mi.endsWith("/A")) match {
          case (true,  true)  => rec(RX) = "AAT-CCG" // ABR1
          case (false, false) => rec(RX) = "CCG-AAT" // BAR2
          case (false,  true) => rec(RX) = "CCG-AAT" // ABR2
          case (true,  false) => rec(RX) = "AAT-CCG" // BAR1
        }
      }

      val in  = builder.toTempFile()
      val out = makeTempFile("duplex.", ".bam")
      new CallDuplexConsensusReads(input=in, output=out, readGroupId="ZZ", cellTag = Some(specialCellTag), threads=threads).execute()
      val reader = SamSource(out)
      val recs = reader.toSeq

      reader.header.getReadGroups should have size 1
      reader.header.getReadGroups.iterator().next().getId shouldBe "ZZ"
      recs should have size 2
      recs.foreach { rec =>
        rec[String](MI) shouldBe "1"
        rec[String](RX) shouldBe (if (rec.firstOfPair) "AAT-CCG" else "CCG-AAT")
        rec[String](specialCellTag) shouldBe "AB"
      }
    }
  }

  it should "produce identical output with one thread and with many when downsampling strands" in {
    val readLength     = 10
    val readsPerStrand = 6
    val maxPerStrand   = 3

    // Build many molecules, each with more reads per strand than we'll allow, and with a mismatching base at a
    // different offset in each read so that which reads are retained is visible in the consensus.  The molecule
    // count exceeds ConsensusCallingIterator's chunk size (threads * 16) so the input spans several chunks rather
    // than relying on the fork-join pool to split a single chunk across threads.
    val builder = new SamBuilder(readLength=readLength, sort=Some(SamOrder.TemplateCoordinate))
    Range.inclusive(1, 300).foreach { mi =>
      Range.inclusive(1, readsPerStrand).foreach { i =>
        val bases = ("A" * readLength).updated(i, 'C')
        builder.addPair(name=f"ab$mi%03d:$i", start1=100, start2=100, attrs=Map(MI -> s"$mi/A"),
          bases1=bases, bases2=bases)
        builder.addPair(name=f"ba$mi%03d:$i", start1=100, start2=100, strand1=Minus, strand2=Plus,
          attrs=Map(MI -> s"$mi/B"), bases1=bases, bases2=bases)
      }
    }
    val in = builder.toTempFile()

    /** Calls duplex consensus reads with the given number of threads and returns the output records. */
    def callWithThreads(threads: Int): Seq[SamRecord] = {
      val out = makeTempFile("duplex.", ".bam")
      new CallDuplexConsensusReads(input=in, output=out, maxReadsPerStrand=Some(maxPerStrand), threads=threads).execute()
      val source = SamSource(out)
      yieldAndThen(source.toIndexedSeq) { source.safelyClose() }
    }

    val singleThreaded = callWithThreads(threads=1)
    val multiThreaded  = callWithThreads(threads=8)

    // Sanity check that the strands really were downsampled, otherwise there's nothing to be non-deterministic about
    singleThreaded should have size 600
    singleThreaded.foreach { rec =>
      rec[Int](ConsensusTags.PerRead.AbRawReadCount) shouldBe maxPerStrand
      rec[Int](ConsensusTags.PerRead.BaRawReadCount) shouldBe maxPerStrand
    }

    multiThreaded.map(_.asSam.getSAMString) should contain theSameElementsInOrderAs singleThreaded.map(_.asSam.getSAMString)
  }
}
