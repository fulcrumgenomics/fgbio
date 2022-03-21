/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric

class CallOverlappingConsensusBasesTest extends UnitSpec {

  "CallOverlappingConsensusBases" should "run end to end" in {
    val builder = new SamBuilder(readLength=10)

    // reads that are not consensus called
    builder.addFrag(start=1) // fragment read
    builder.addPair(start1=1, unmapped2=true) // paired end read with R2 unmapped
    builder.addPair(start2=1, unmapped1=true) // paired end read with R1 unmapped
    builder.addPair(start1=1, start2=11) // no overlap (they abut)

    // reads that are consensus called
    builder.addPair(bases1="A"*10, bases2="A"*10, start1=1, start2=1) // fully overlapping, complete agreement
    builder.addPair(bases1="A"*10, bases2="C"*10, start1=1, start2=10) // 1bp overlap, disagreement
    builder.addPair(bases1="A"*10, bases2="C"*10, start1=9, start2=1, strand2=SamBuilder.Plus) // odd pair orientation, extend past mate, overlaps 2bp, disagreement

    // run the tool
    val output  = makeTempFile("test.", ".bam")
    val metrics = makeTempFile("test.", ".bam")
    val tool    = new CallOverlappingConsensusBases(
      input   = builder.toTempFile(),
      output  = output,
      metrics = metrics
    )
    tool.execute()

    // metrics
    val Seq(templateStats, basesStats) = Metric.read[CallOverlappingConsensusBasesMetric](metrics).toIndexedSeq
    templateStats shouldBe CallOverlappingConsensusBasesMetric(
      tpe         = CountType.Templates,
      total       = 7,
      overlapping = 3,
      corrected   = 2
    )
    basesStats shouldBe CallOverlappingConsensusBasesMetric(
      tpe         = CountType.Bases,
      total       = 1*10 + 6*10*2, // 1 fragment, 6 paired end, so 13 reads of 10bp each
      overlapping = 26, // 1 PE with 10bp overlap, 1 PE with 1bp overlap, 1PE with 2bp overlap; 2*13bp
      corrected   = 6   // 1 PE with 1bp corrected, 1PE with 2bp corrected; 2*3
    )
  }
}
