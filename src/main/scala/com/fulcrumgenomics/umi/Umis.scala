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

object Umis {

  /** Copies the UMI sequence from the read name.
    *
    * The read name is split by the given name delimiter, and the last field is assumed to be the UMI sequence.  The UMI
    * will be copied to the `RX` tag as per the SAM specification.
    *
    * An exception will be thrown if the read name does not contain a valid UMI in the last delimited segment.
    *
    * @param rec the record to modify
    * @param removeUmi true to remove the UMI from the read name, otherwise only copy the UMI to the tag
    * @param delimiter the delimiter of fields within the read name
    * @return the modified record
    */
  def copyUmiFromReadName(rec: SamRecord, removeUmi: Boolean = false, delimiter: Char = ':'): SamRecord = {
    // Extract and set the UMI
    val umi = extractUmisFromReadName(rec.name, delimiter, strict=false)
    require(umi.nonEmpty, f"No valid UMI found in: ${rec.name}")
    umi.foreach(u => rec(ConsensusTags.UmiBases) = u)

    // Remove the UMI from the read name if requested
    if (removeUmi) rec.name = rec.name.substring(0, rec.name.lastIndexOf(delimiter))

    rec
  }

  /**
    * Extracts the UMI from an Illumina fastq style read name.  Illumina documents their FASTQ read names as:
    *   @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
    *
    *  See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    *  The UMI field is optional, so read names may or may not contain it.  Illumina also specifies that the UMI
    *  field may contain multiple UMIs, in which case they will delimit them with `+` characters.  Pluses will be
    *  translated to hyphens before returning.
    *
    *  If `strict` is true the name _must_ contain either 7 or 8 colon-separated segments, 
    with the UMI being the last in the case of 8 and `None` in the case of 7.
    * 
    * If `strict` is false the last segment is returned so long as it appears to be a valid UMI.
    */
  def extractUmisFromReadName(name: String, delimiter: Char = ':', strict: Boolean): Option[String] = {
    // If strict, check that the read name actually has eight parts, which is expected
    val rawUmi = if (strict) {
      val colons = name.count(_ == delimiter)
      if (colons == 6) None
      else if (colons == 7) Some(name.substring(name.lastIndexOf(delimiter) + 1, name.length))
      else throw new IllegalArgumentException(s"Trying to extract UMI from read with ${colons + 1} parts (7-8 expected): ${name}")
    } else {
      val idx = name.lastIndexOf(delimiter)
      require(idx != -1, s"Read did not have multiple '${delimiter}'-separated fields: ${name}")
      Some(name.substring(idx + 1, name.length))
    }

    val umi   = rawUmi.map(raw => (if (raw.indexOf('+') > 0) raw.replace('+', '-') else raw).toUpperCase)
    val valid = umi.forall(u => u.forall(isValidUmiCharacter))

    if (strict && !valid) throw new IllegalArgumentException(s"Invalid UMI '${umi.get}' extracted from name '${name}")
    else if (!valid) None
    else umi
  }

  @inline private def isValidUmiCharacter(ch: Char): Boolean = {
    ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N' || ch == '-'
  }
}
