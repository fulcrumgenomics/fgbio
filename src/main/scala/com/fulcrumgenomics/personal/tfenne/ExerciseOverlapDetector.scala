/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

import java.nio.file.Path

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.{Interval, OverlapDetector}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

import scala.util.Random


@clp(group=ClpGroups.Personal, description="Test out overlap detector.")
class ExerciseOverlapDetector
( @arg(flag='l', doc="Interval size") val length: Int = 20,
  @arg(flag='n', doc="Number of things.") val nThings: Int = 1000000,
  @arg(flag='N', doc="Number of hits per thing") val hitsPerThing: Int = 50,
  @arg(flag='d', doc="Sequence dictionary.") val sequenceDictionary: Path
) extends FgBioTool with LazyLogging {
  override def execute(): Unit = {
    val od     = new OverlapDetector[Interval](0, 0)
    val sd     = SAMSequenceDictionaryExtractor.extractDictionary(sequenceDictionary)
    val random = new Random(42)

    logger.info("Assembling things.")
    val things = Range(0, nThings).map { _ =>
      Range(0, hitsPerThing).map { _ =>
        val chrIdx = random.nextInt(sd.size())
        val chrom  = sd.getSequence(chrIdx)
        val start  = random.nextInt(chrom.getSequenceLength - length)
        new Interval(chrom.getContig, start, start + length - 1)
      }
    }

    logger.info("Querying things.")
    val progress = new ProgressLogger(logger, "things", "queried", 10000)
    things.foreach { thing =>
      if (thing.exists(r => !od.overlapsAny(r))) {
        thing.foreach(r => od.addLhs(r, r))
      }

      progress.record()
    }
  }
}
