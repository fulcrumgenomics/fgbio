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

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.testing.UnitSpec

class InteractiveAlignmentScorerTest extends UnitSpec {

  private case class Bases(seq: String) extends Alignable {
    val bases: Array[Byte] = seq.getBytes
  }

  private def f(query: String, target: String): AlignmentTask[Bases, Bases] = AlignmentTask(Bases(seq = query), Bases(seq = target))

  private val testTasks =  Seq(f("ACGTAACC", "ACGTAACC"), f("CCGG", "AACCGGTT"), f("AAAATTTT", "AAAATTTTGGGG"))

  Seq(Mode.Local, Mode.Glocal, Mode.Global).foreach { mode =>
    "InteractiveAlignmentScorer" should s"align in $mode mode" in {
      val iAligner = InteractiveAlignmentScorer[Bases, Bases](1, -4, -6, -1, mode=mode)
      val aligner  = Aligner(1, -4, -6, -1, mode=mode)

      // Add the tasks to the interactive aligner
      testTasks.foreach { task => iAligner.append(task) }

      // Compare the results
      testTasks.zip(iAligner.iterator.toSeq).foreach { case (task, actualPartialAlignment) =>
        val expectedAlignment        = aligner.align(task.query.bases, task.target.bases)
        val expectedPartialAlignment = PartialAlignment(expectedAlignment)

        actualPartialAlignment shouldBe expectedPartialAlignment
      }
      iAligner.numAdded shouldBe testTasks.length
      iAligner.numAligned shouldBe testTasks.length
    }
  }
}
