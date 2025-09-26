/*
 * The MIT License
 *
 * Copyright (c) 2023 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamOrder, SamSource}
import com.fulcrumgenomics.bam.pileup.PileupBuilder
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

class DownsampleAndNormalizeBamTest extends UnitSpec {
  def newBuilder(readLength: Int = 50): SamBuilder = {
    /** Makes a new builder with two short chromosomes. */
    new SamBuilder(
      readLength = readLength,
      sort       = Some(SamOrder.Coordinate),
      sd         = Some(SequenceDictionary(
        SequenceMetadata(name="chr1", length=1000, index=0),
        SequenceMetadata(name="chr2", length=1000, index=0),
      ))
    )
  }

  /** Generates a temp file and returns its path. */
  def newBam(): PathToBam = makeTempFile("downsampled.", ".bam")


  "DownsampleAndNormalizeBam" should "run ok on an empty bam" in {
    val empty = newBuilder().toTempFile()
    val out   = newBam()
    new DownsampleAndNormalizeBam(input=empty, output=out, coverage=100).execute()
    val recs  = SamSource(out).toIndexedSeq
    recs.isEmpty shouldBe true
  }

  it should "downsample a BAM with excess coverage" in {
    val builder  = newBuilder()
    val out      = newBam()
    val coverage = 30

    // Add lots of reads - 25 read pairs starting at each position until we run out of room
    builder.dict.foreach { chrom =>
      Range.inclusive(1, chrom.length - 125).foreach { start =>
        forloop (from=0, until=5) { _ =>
          val _ = builder.addPair(contig=chrom.index, start1=start , start2=start+75)
        }
      }
    }

    new DownsampleAndNormalizeBam(input = builder.toTempFile(), output = out, coverage = coverage).execute()

    val pileIn  = PileupBuilder(source=builder.toSource, accessPattern=PileupBuilder.BamAccessPattern.Streaming)
    val pileOut = PileupBuilder(source=SamSource(out), accessPattern=PileupBuilder.BamAccessPattern.Streaming)
    var totalIn  = 0L
    var totalOut = 0L

    builder.dict.foreach { chrom =>
      Range.inclusive(1, chrom.length).foreach { pos =>
        val covIn  = pileIn.pileup(chrom.name, pos).withoutOverlaps.depth
        val covOut = pileOut.pileup(chrom.name, pos).withoutOverlaps.depth
        totalIn += covIn
        totalOut += covOut


        covOut <= covIn shouldBe true
        if (covIn < coverage) {
          covOut shouldBe covIn
        }
        else {
          covOut should be >= coverage
        }
      }
    }

    val meanIn  = totalIn / builder.dict.sumBy(_.length).toDouble
    val meanOut = totalOut / builder.dict.sumBy(_.length).toDouble

    meanOut should be < meanIn
    meanOut should be < (coverage * 2).toDouble
    meanOut should be > coverage.toDouble
  }
}
