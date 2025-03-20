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
 *
 */

package com.fulcrumgenomics.basecalling

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.illumina.RunInfo
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.{Metric, ReadStructure}
import htsjdk.samtools.util.Iso8601Date

class ExtractIlluminaRunInfoTest extends UnitSpec {
  import com.fulcrumgenomics.illumina.RunInfoTest.runInfo

  "ExtractRunInfo" should "run end-to-end" in {
    val in = runInfo(date="20170204", readStructure=ReadStructure("8B150T"))
    val out = makeTempFile("out", ".txt")
    new ExtractIlluminaRunInfo(input=in, output=out).execute()
    val metrics = Metric.read[RunInfo](out)
    metrics should have size 1
    val metric  = metrics.head

    metric.run_barcode             shouldBe "BCDEFGHIJ_NS123456"
    metric.flowcell_barcode        shouldBe "NS123456"
    metric.instrument_name         shouldBe "BCDEFGHIJ"
    metric.run_date                shouldBe new Iso8601Date("2017-02-04")
    metric.read_structure.toString shouldBe "8B150T"
    metric.num_lanes               shouldBe 4
  }
}
