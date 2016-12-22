/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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
import htsjdk.samtools.util.SamLocusIterator
import htsjdk.samtools.util.SamLocusIterator.LocusInfo
import scala.collection.JavaConversions.iterableAsScalaIterable

object BaseCounts {
  /**
    * Generates a BaseCounts object from a collection of LocusInfo object
    * representing a pileup.  Ensures that each template is only counted once.
    */
  def apply(locus: LocusInfo): BaseCounts = apply(locus.getRecordAndOffsets)

  /**
    * Generates a BaseCounts object from a collection of RecordAndOffest objects
    * representing a pileup.  Ensures that each template is only counted once.
    */
  def apply(recs: Traversable[SamLocusIterator.RecordAndOffset]): BaseCounts = {
    val countsByBase = recs.groupBy(r => Character.toUpperCase(r.getReadBase.toChar))
      .map { case (ch, rs) => (ch, rs.map(_.getRecord.getReadName).toSeq.distinct.size) }

    BaseCounts(
      a=countsByBase.getOrElse('A', 0),
      c=countsByBase.getOrElse('C', 0),
      g=countsByBase.getOrElse('G', 0),
      t=countsByBase.getOrElse('T', 0),
      n=countsByBase.getOrElse('N', 0))
  }
}

/** Utility class to count the number of observation of each kind of base. */
case class BaseCounts(a:Int , c: Int, g: Int, t: Int, n: Int) {
  /** Returns the number of bases counted of the given base. Base must be in [acgtnACGTN]. */
  def apply(ch: Char): Int = ch.toUpper match {
        case 'A' => a
        case 'C' => c
        case 'G' => g
        case 'T' => t
        case 'N' => n
        case x   => throw new IllegalArgumentException(s"'$x' is not a valid base.")
  }
}
