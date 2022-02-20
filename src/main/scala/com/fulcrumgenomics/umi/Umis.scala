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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.SamRecord
import htsjdk.samtools.util.SequenceUtil

object Umis {

  /** Copies the UMI sequence from the read name.
    *
    * The read name is split by the given name delimiter, and the last field is assumed to be the UMI sequence.  The UMI
    * will be copied to the `RX` tag as per the SAM specification.
    *
    * @param rec the record to modify
    * @param removeUmi true to remove the UMI from the read name, otherwise only copy the UMI to the tag
    * @param umiDelimiter if not None, replaces any occurrences of this delimiter found in the UMI with a dash ('-')
    *                     as per the SAM specification
    * @return the modified record
    */
  def copyUmiFromReadName(rec: SamRecord,
                          removeUmi: Boolean = false,
                          nameDelimiter: Char = ':',
                          umiDelimiter: Option[Char] = Some('+')): SamRecord = {
    // extract the UMI
    val idx = rec.name.lastIndexOf(nameDelimiter)
    require(idx != -1, s"Read did not have multiple '$nameDelimiter'-separated fields: ${rec.name}")
    val rawUmi = rec.name.substring(idx + 1, rec.name.length)

    // re-delimit the UMI if desired, validating as necessary
    val umi = umiDelimiter match {
      case Some(delim) => rawUmi.replace(delim, '-')
      case None        => rawUmi
    }

    // validate the UMI
    require(umi.nonEmpty && umi.forall { char => SequenceUtil.isValidBase(char.toByte) || SequenceUtil.isNoCall(char.toByte) || char == '-' },
      f"UMI '$umi' is not valid (not `[ACGTNacgtn-]+`) for read: ${rec.name}"
    )

    // update the record
    rec(ConsensusTags.UmiBases) = umi
    if (removeUmi) {
      rec.name = rec.name.substring(0, idx)
    }
    rec
  }
}
