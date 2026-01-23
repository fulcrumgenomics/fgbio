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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._

import htsjdk.samtools.util.{CoordMath, Locatable}

/**
  * Trait that can be added to any class that defines chrom/start/end to implement Locatable
  */
trait LocatableAdapter[A <: LocatableAdapter[A]] extends Locatable with Ordered[A] {
  def chrom: String
  def start: Int
  def end: Int
  def name: String
  def negativeStrand: Boolean
  def positiveStrand: Boolean = !negativeStrand

  // Implementations of Locatable interface methods
  final override def getContig: String = chrom
  final override def getStart: Int     = start
  final override def getEnd: Int       = end

  final def length: Int = getLengthOnReference

  /** Calculates the overlap between two locatables assuming that they do in fact overlap. */
  def overlap(other: Locatable): Int = CoordMath.getOverlap(start, end, other.getStart, other.getEnd)

  /** Provides a basic comparison that orders on chrom (lexicographically) and then start and end position. */
  def compare(that: A): Int = {
    var result = this.chrom.compareTo(that.chrom)
    if (result == 0) result = this.start - that.start
    if (result == 0) result = this.end - that.end
    result
  }
}
