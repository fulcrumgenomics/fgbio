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

/**
  * Object that encapsulates the various consensus related tags that are added to consensus reads
  * at both the per-read and per-base level.
  *
  * Currently only contains tags for single-strand consensus reads, but with a view to using the following
  * names for consistency if/when we add duplex calling:
  *     Value                 AB  BA  Final
  *     ===================== ==  ==  =====
  *     bases                 ac  bc  bases
  *     quals                 aq  bq  quals
  *     per-read-depth        aD  bD  cD
  *     per-read-min-depth    aM  bM  cM
  *     per-read-error-rate   aE  bE  cE
  *     per-base-depth        ad  bd  cd
  *     per-base-error-count  ae  be  ce
  */
object ConsensusTags {
  /** The default field in which to look for UMI sequences. */
  val UmiBases    = "RX"

  /** Post-grouping ID that is file-wide unique per source molecule. */
  val MolecularId = "MI"

  object PerBase {
    /** The per-base number of raw-reads contributing to the consensus (stored as a short[]). */
    val RawReadCount    = "cd" // consensus depth
    /** The number of bases at each position that disagreed with the final consensus call (stored as a short[]). */
    val RawReadErrors   = "ce" // consensus errors

    // Duplex versions of the above tags for the two single strand consensus reads
    val AbRawReadCount  = "ad"
    val BaRawReadCount  = "bd"
    val AbRawReadErrors = "ae"
    val BaRawReadErrors = "be"

    // NOTE: Important that this is updated if any new tags are added!
    val AllPerBaseTags = Seq(RawReadCount, RawReadErrors, AbRawReadCount, AbRawReadErrors, BaRawReadCount, BaRawReadErrors)
  }

  object PerRead {
    /** The number of raw reads that contributed to the consensus. I.e. the maximum of the per-base raw-read counts. */
    val RawReadCount    = "cD" // consensus Depth
    /** The minimum number of raw reads contributing to a consensus call anywhere in the read. */
    val MinRawReadCount = "cM" // consensus Min-depth
    /** The number of bases in the raw reads that contributed to the consensus but disagreed with the consensus call. */
    val RawReadErrorRate   = "cE" // consensus Error rate

    // Duplex versions of the above tags for the two single strand consensus reads
    val AbRawReadCount     = "aD"
    val BaRawReadCount     = "bD"
    val AbMinRawReadCount  = "aM"
    val BaMinRawReadCount  = "bM"
    val AbRawReadErrorRate = "aE"
    val BaRawReadErrorRate = "bE"
  }
}
