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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.SafelyClosable
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io

@clp(
  description = "Removes SAM tags from a SAM or BAM file.  If no tags to remove are given, the original file is produced.",
  group = ClpGroups.SamOrBam
)
class RemoveSamTags
( @arg(flag='i', doc = "Input SAM or BAM.") val input: PathToBam,
  @arg(flag='o', doc = "Output SAM or BAM.") val output: PathToBam,
  @arg(flag='t', doc = "The tags to remove.", minElements = 0) val tagsToRemove: Seq[String] = Seq.empty
) extends FgBioTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  validate(tagsToRemove.forall(_.length == 2), "Tags to remove must be of length two: " + tagsToRemove.filter(_.length != 2).mkString(", "))

  override def execute(): Unit = {
    val in  = SamSource(input)
    val out = SamWriter(output, in.header)
    in.foreach { rec =>
      tagsToRemove.foreach { tag => rec(tag) = null }
      out += rec
    }
    in.safelyClose()
    out.close()
  }
}
