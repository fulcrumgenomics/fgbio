/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics LLC
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
package com.fulcrumgenomics.testing

import com.fulcrumgenomics.bam.api.SamRecord

/** Test fixtures for dovetail FR read pairs whose 5' ends coincide. */
object DovetailPairs {

  /** The cigar of the forward-strand read (R1) of [[pair]]. */
  val ForwardCigar: String = "68S53M8S"

  /** The cigar of the reverse-strand read (R2) of [[pair]]. */
  val ReverseCigar: String = "28S48M53S"

  /** Returns a dovetail FR pair as (R1, R2) whose 5' ends coincide at position 96 on the first contig: R1 is on the
    * forward strand and aligned to 96-148 (`68S53M8S`), R2 is on the reverse strand and aligned to 49-96
    * (`28S48M53S`).
    *
    * When `withMateCigar` is false, the MC tags are removed and TLEN is set to 0, as bwa writes it for this pair. The
    * per-record `SamRecord.isFrPair` then reports false for R1 and true for R2.
    */
  def pair(withMateCigar: Boolean, name: String = "q1"): (SamRecord, SamRecord) = {
    val builder = new SamBuilder(readLength=129)
    val Seq(r1, r2) = builder.addPair(name=name, contig=0, start1=96, start2=49, cigar1=ForwardCigar, cigar2=ReverseCigar).toSeq
    if (!withMateCigar) Seq(r1, r2).foreach { r => r.remove("MC"); r.insertSize = 0 }
    (r1, r2)
  }

  /** Returns two copies of [[pair]]: the first with MC tags, the second without. */
  def withAndWithoutMateCigar(): ((SamRecord, SamRecord), (SamRecord, SamRecord)) = {
    (pair(withMateCigar=true, name="with_mc"), pair(withMateCigar=false, name="without_mc"))
  }
}
