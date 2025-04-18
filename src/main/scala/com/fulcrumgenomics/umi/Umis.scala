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
import com.fulcrumgenomics.umi.ConsensusTags.PerRead.{AbRawReadCount, BaRawReadCount, RawReadCount}
import com.fulcrumgenomics.util.Sequences

object Umis {
  val RevcompPrefix: String = "r"

  /** Copies the UMI sequence from the read name.
    *
    * The read name is split by the given name delimiter, and the last field is assumed to be the UMI sequence.  The UMI
    * will be copied to the `RX` tag as per the SAM specification.
    *
    * An exception will be thrown if the read name does not contain a valid UMI in the last delimited segment.
    *
    * @param rec the record to modify
    * @param removeUmi true to remove the UMI from the read name, otherwise only copy the UMI to the tag
    * @param fieldDelimiter the delimiter of fields within the read name
    * @param umiDelimiter the delimiter between sequences in the UMI string
    * @param reverseComplementPrefixedUmis whether to reverse-compliment UMIs prefixed with 'r'
    * @return the modified record
    */
  def copyUmiFromReadName(rec: SamRecord,
                          removeUmi: Boolean = false,
                          fieldDelimiter: Char = ':',
                          umiDelimiter: Char = '+',
                          reverseComplementPrefixedUmis: Boolean = true): SamRecord = {
    // Extract and set the UMI
    val umi = extractUmisFromReadName(
      name                          = rec.name, 
      fieldDelimiter                = fieldDelimiter,
      strict                        = false,
      umiDelimiter                  = umiDelimiter,
      reverseComplementPrefixedUmis = reverseComplementPrefixedUmis,
    )
    require(umi.nonEmpty, f"No valid UMI found in: ${rec.name}")
    umi.foreach(u => rec(ConsensusTags.UmiBases) = u)

    // Remove the UMI from the read name if requested
    if (removeUmi) rec.name = rec.name.substring(0, rec.name.lastIndexOf(fieldDelimiter.toInt))

    rec
  }

  /**
    * Extracts the UMI from an Illumina fastq style read name.  Illumina documents their FASTQ read names as:
    *   `@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>`
    *
    *  See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    *  The UMI field is optional, so read names may or may not contain it.  Illumina also specifies that the UMI
    *  field may contain multiple UMIs, in which case they will delimit them with `umiDelimiter` characters, which
    *  will be translated to hyphens before returning.
    *
    *  If `strict` is true the name _must_ contain either 7 or 8 colon-separated segments, 
    with the UMI being the last in the case of 8 and `None` in the case of 7.
    * 
    * If `strict` is false the last segment is returned so long as it appears to be a valid UMI.
    */
  def extractUmisFromReadName(name: String,
                              fieldDelimiter: Char = ':',
                              strict: Boolean,
                              umiDelimiter: Char = '+',
                              reverseComplementPrefixedUmis: Boolean = true): Option[String] = {
    // If strict, check that the read name actually has eight parts, which is expected
    val rawUmi = if (strict) {
      val colons = name.count(_ == fieldDelimiter)
      if (colons == 6) None
      else if (colons == 7) Some(name.substring(name.lastIndexOf(fieldDelimiter.toInt) + 1, name.length))
      else throw new IllegalArgumentException(s"Trying to extract UMI from read with ${colons + 1} parts (7-8 expected): ${name}")
    } else {
      val idx = name.lastIndexOf(fieldDelimiter.toInt)
      require(idx != -1, s"Read did not have multiple '${fieldDelimiter}'-separated fields: ${name}")
      Some(name.substring(idx + 1, name.length))
    }

    // Remove 'r' prefixes, optionally reverse-complementing the prefixed UMIs if 
    // reverseComplementPrefixedUmis = true, replace the delimiter (if any) with '-',
    // and make sure the sequence is upper-case.
    val umi = rawUmi.map(raw =>
      (raw.indexOf(RevcompPrefix) >= 0, raw.indexOf(umiDelimiter.toInt) > 0) match {
        case (true, true) if reverseComplementPrefixedUmis => 
          raw.split(umiDelimiter).map(seq => 
            if (seq.startsWith(RevcompPrefix)) Sequences.revcomp(seq.stripPrefix(RevcompPrefix))
            else seq.stripPrefix(RevcompPrefix)
          ).mkString("-").toUpperCase
        case (true, false) if reverseComplementPrefixedUmis =>
          Sequences.revcomp(raw.stripPrefix(RevcompPrefix)).toUpperCase
        case (true, true) => raw.replace(RevcompPrefix, "").replace(umiDelimiter, '-').toUpperCase
        case (true, false) => raw.replace(RevcompPrefix, "").toUpperCase
        case (false, true) => raw.replace(umiDelimiter, '-').toUpperCase
        case (false, false) => raw.toUpperCase
      }
    )

    val valid  = umi.forall(u => u.forall(isValidUmiCharacter))

    if (strict && !valid) throw new IllegalArgumentException(s"Invalid UMI '${umi.get}' extracted from name '${name}")
    else if (!valid) None
    else umi
  }

  @inline private def isValidUmiCharacter(ch: Char): Boolean = {
    ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N' || ch == '-'
  }

  /** Returns True if the record appears to be a consensus sequence typically produced by fgbio's
    * CallMolecularConsensusReads or CallDuplexConsensusReads.
    *
    * @param rec the record to check
    * @return boolean indicating if the record is a consensus or not
    */
  def isFgbioStyleConsensus(rec: SamRecord): Boolean = isFgbioSimplexConsensus(rec) || isFgbioDuplexConsensus(rec)

  /** Returns true if the record appears to be a simplex consensus sequence. */
  def isFgbioSimplexConsensus(rec: SamRecord): Boolean = rec.contains(RawReadCount) && !isFgbioDuplexConsensus(rec)

  /** Returns true if the record appears to be a duplex consensus sequence. */
  def isFgbioDuplexConsensus(rec: SamRecord): Boolean = rec.contains(AbRawReadCount) && rec.contains(BaRawReadCount)
}
