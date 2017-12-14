/*
 * The MIT License
 *
 * Copyright (c) $year Fulcrum Genomics
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

import java.nio.file.Path
import java.text.DecimalFormat

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.PathToBam
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import htsjdk.samtools.util.IntervalList

import scala.math.abs

@clp(description = 
  """
     |Makes reads unmapped if based on their insert size.
     |
     |Fragment reads will always be made unmapped. Secondary and supplementary reads whose become unmapped will instead
     |be removed.
  """,
  group = ClpGroups.SamOrBam)
class MakeUnmappedByInsertSize
( @arg(flag='i', doc="Input BAM file.")                                 val input: PathToBam,
  @arg(flag='o', doc="Output BAM file.")                                val output: PathToBam,
  @arg(flag='m', doc="Remove all reads with insert size < this value.") val minInsertSize: Option[Int] = None,
  @arg(flag='M', doc="Remove all reads with insert size > this value.") val maxInsertSize: Option[Int] = None
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  validate(minInsertSize.isDefined || maxInsertSize.isDefined, "Either --min-insert-size or --max-insert-size must be given.")

  override def execute(): Unit = {
    val in       = SamSource(input)
    val out      = SamWriter(output, in.header)
    var secondary: Long = 0
    var supplementary: Long = 0
    val unmapped: Long = in
      .count { rec =>
        val insertSize = abs(rec.insertSize)
        val makeUnmapped = !rec.paired ||
          minInsertSize.exists(isize => insertSize < isize) ||
          maxInsertSize.exists(isize => insertSize > isize)
        if (makeUnmapped && rec.mapped) {
          if (rec.secondary) {
            secondary += 1
            false
          }
          else if (rec.supplementary) {
            supplementary += 1
            false
          }
          else {
            rec.unmapped = true
            if (rec.paired) rec.unmapped = true
            rec.mapq = 0
            out += rec
            true
          }
        }
        else {
          out += rec
          false
        }
    }

    logger.info(f"Made $unmapped%,d records unmapped.")
    logger.info(f"Removed $secondary%,d secondary and $supplementary%,d supplementary records.")

    out.close()
    in.safelyClose()
  }
}
