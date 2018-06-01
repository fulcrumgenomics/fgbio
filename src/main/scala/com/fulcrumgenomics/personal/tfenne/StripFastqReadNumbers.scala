/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.fastq.{FastqSource, FastqWriter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io

@clp(group=ClpGroups.Personal, description=
  """
    |Removes trailing /# from read names in fastq.  Read names that end in a slash and
    |a single digit will have the slash and digit removed.  Read names that do not end
    |with a slash and a number will be untouched.
  """)
class StripFastqReadNumbers
( @arg(flag='i', doc="Input fastq file.") val input: PathToFastq,
  @arg(flag='o', doc="Output fastq file.") val output: PathToFastq
) extends FgBioTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  override def execute(): Unit = {
    val in = FastqSource(input)
    val out = FastqWriter(output)

    in.foreach { rec =>
      out += rec.copy(readNumber=None)
    }

    in.safelyClose()
    out.close()
  }
}
