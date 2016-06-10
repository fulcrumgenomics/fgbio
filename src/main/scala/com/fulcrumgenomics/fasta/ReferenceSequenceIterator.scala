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
package com.fulcrumgenomics.fasta

import htsjdk.samtools.reference.{ReferenceSequence, ReferenceSequenceFile, ReferenceSequenceFileFactory}
import dagr.commons.CommonsDef._

object ReferenceSequenceIterator {
  /** Constructs an iterator over a reference sequence from a Path to the FASTA file. */
  def apply(path: PathToFasta, stripComments: Boolean = false) = {
    new ReferenceSequenceIterator(ReferenceSequenceFileFactory.getReferenceSequenceFile(path, stripComments, false))
  }

  /** Constructs an iterator over a reference sequence from a ReferenceSequenceFile. */
  def apply(ref: ReferenceSequenceFile) = {
    new ReferenceSequenceIterator(ref)
  }
}

/**
  * Simple Iterator that wraps a ReferenceSequence and provides a standard Iterator over
  * all the reference sequences in order.
  */
class ReferenceSequenceIterator private(private val ref: ReferenceSequenceFile) extends Iterator[ReferenceSequence] {
  // The next reference sequence; will be null if there is no next :(
  private var nextSequence: ReferenceSequence = ref.nextSequence()

  override def hasNext: Boolean = this.nextSequence != null

  override def next(): ReferenceSequence = {
    assert(hasNext, "next() called on empty iterator")
    yieldAndThen(this.nextSequence) { this.nextSequence = ref.nextSequence() }
  }
}
