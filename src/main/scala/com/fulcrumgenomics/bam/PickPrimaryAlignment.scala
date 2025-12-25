/*
 * The MIT License
 *
 * Copyright (c) 2024 Fulcrum Genomics LLC
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
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger

@clp(group=ClpGroups.SamOrBam, description="""
Resets which record is marked as the primary alignment per read.
""")
class PickPrimaryAlignment(
  @arg(flag='i', doc="") val input: PathToBam,
  @arg(flag='o', doc="") val output: PathToBam
) extends FgBioTool {

  override def execute(): Unit = {
    val in       = SamSource(input)
    val out      = SamWriter(output, in.header)
    val progress = ProgressLogger(this.logger, noun="templates")

    Bams.templateIterator(in).foreach { template =>
      val fixed = pickPrimaries(template)
      if (!(fixed eq template)) {
        fixed.fixMateInfo()
      }

      out ++= fixed.allReads
      progress.record()
    }

    in.safelyClose()
    out.close()
  }

  /** Re-picks the primary mappings for any reads that have supplementary reads. */
  private[bam] def pickPrimaries(t: Template): Template = {
    if (t.r1Supplementals.isEmpty && t.r2Supplementals.isEmpty) {
      t
    }
    else {
      val r1Primary = if (t.r1Supplementals.isEmpty) t.r1 else {
        Some(pickPrimary(t.r1 ++ t.r1Supplementals))
      }
      val r2Primary = if (t.r2Supplementals.isEmpty) t.r2 else {
        Some(pickPrimary(t.r2 ++ t.r2Supplementals))
      }

      if ((r1Primary eq t.r1) && (r2Primary eq t.r2)) t else {
        r1Primary.foreach(_.supplementary = false)
        r2Primary.foreach(_.supplementary = false)

        Template(
          r1 = r1Primary,
          r2 = r2Primary,
          r1Supplementals = (t.r1 ++ t.r1Supplementals).filter(r => !r1Primary.exists(_ eq r)).tapEach(_.supplementary=true).toSeq,
          r2Supplementals = (t.r2 ++ t.r2Supplementals).filter(r => !r2Primary.exists(_ eq r)).tapEach(_.supplementary=true).toSeq,
          r1Secondaries   = t.r1Secondaries,
          r2Secondaries   = t.r2Secondaries
        )
      }
    }
  }

  /** Selects the record, among mappings of the same read, with the earliest query base aligned. */
  private [bam] def pickPrimary(recs: Iterable[SamRecord]): SamRecord = {
    recs.minBy(firstMappedBaseInQueryOrder)
  }

  /** Returns the 1-based position of the first query base that is aligned to the reference. */
  private [bam] def firstMappedBaseInQueryOrder(rec: SamRecord): Int = {
    require(rec.mapped, s"Record is unmapped: ${rec}")
    val iter = (if (rec.negativeStrand) rec.cigar.elems.reverseIterator else rec.cigar.elems.iterator).bufferBetter
    iter.takeWhile(!_.operator.isAlignment).sumBy(_.lengthOnQuery)
  }
}
