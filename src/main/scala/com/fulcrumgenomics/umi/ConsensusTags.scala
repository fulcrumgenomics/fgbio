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
  *     per-read-depth        aD  bD  cD
  *     per-read-min-depth    aM  bM  cM
  *     per-read-error-rate   aE  bE  cE
  *     per-base-depth        ad  bd  cd
  *     per-base-error-count  ae  be  ce
  *     per-base-bases        ac  bc  bases
  *     per-base-quals        aq  bq  quals
  * The second letter in the tag is lower case if it is per-base, upper case if it is per-read.
  */
object ConsensusTags {
  /** The default field in which to look for UMI sequences. */
  val UmiBases    = "RX"

  /** The field in which the original UMI bases are stored post-correction. */
  val OriginalUmiBases = "OX"

  /** Post-grouping ID that is file-wide unique per source molecule. */
  val MolecularId = "MI"

  object PerBase {
    /** The per-base number of raw-reads contributing to the consensus (stored as a short[]). */
    val RawReadCount    = "cd" // consensus depth
    /** The number of bases at each position that disagreed with the final consensus call (stored as a short[]). If the
      * final consensus call is a no call (N), then we use the most likely consensus base instead of the final call.  */
    val RawReadErrors   = "ce" // consensus errors

    // Duplex versions of the above tags for the two single strand consensus reads
    val AbRawReadCount  = "ad"
    val BaRawReadCount  = "bd"
    val AbRawReadErrors = "ae"
    val BaRawReadErrors = "be"

    // Duplex-specific tags
    /** The single-stranded consensus from the AB raw reads. */
    val AbConsensusBases = "ac"
    /** The single-stranded consensus from the BA raw reads. */
    val BaConsensusBases = "bc"
    /** The phred-scaled qualities, as phred-33 ascii values, of the single-stranded consensus from the AB raw reads. */
    val AbConsensusQuals = "aq"
    /** The phred-scaled qualities, as phred-33 ascii values, of the single-stranded consensus from the BA raw reads. */
    val BaConsensusQuals = "bq"

    /** The set of the per-base tags produced by the consensus caller that need to be reversed after alignment. */
    val TagsToReverse = Seq(RawReadCount, RawReadErrors, AbRawReadCount, AbRawReadErrors, BaRawReadCount, BaRawReadErrors,
      AbConsensusQuals, BaConsensusQuals)

    /** The set of the per-base tags produced by the consensus caller that need to be reverse complemented after alignment. */
    val TagsToReverseComplement = Seq(AbConsensusBases, BaConsensusBases)

    // NOTE: Important that this is updated if any new tags are added!
    /** The set of all per-base tags. */
    val AllPerBaseTags = TagsToReverse ++ TagsToReverseComplement
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

    // NOTE: Important that this is updated if any new tags are added!
    /** The set of all per-read tags. */
    val AllPerReadTags = Seq(
      RawReadCount,
      MinRawReadCount,
      RawReadErrorRate,
      AbRawReadCount,
      BaRawReadCount,
      AbMinRawReadCount,
      BaMinRawReadCount,
      AbRawReadErrorRate,
      BaRawReadErrorRate
    )
  }
}
