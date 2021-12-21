/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics LLC
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
import com.fulcrumgenomics.alignment.{Aligner, Mode}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.fastq.{FastqSource, FastqWriter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Sequences}

@clp(group=ClpGroups.Personal, description=
  """
    |Align one query sequence to one target sequence.
  """)
class PairwiseAlign
( @arg(flag='q', doc="Query sequence.")  val query: String,
  @arg(flag='t', doc="Target sequence.") val target: String,
  @arg(flag='m', doc="Alignment mode.") val mode: Mode = Mode.Global,
) extends FgBioTool {

  override def execute(): Unit = {
    val aligner = Aligner(1, -4, -6, -1, mode)

    val left  = aligner.align(query, target)
    val right = aligner.align(Sequences.revcomp(query), Sequences.revcomp(target)).revcomped

    println("Left Aligned:")
    left.paddedString().foreach(println)

    println()
    println("Right Aligned:")
    right.paddedString().foreach(println)
  }
}
