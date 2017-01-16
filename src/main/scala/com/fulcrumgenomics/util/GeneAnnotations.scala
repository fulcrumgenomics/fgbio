/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import htsjdk.samtools.util.{CoordMath, Interval}

/** Stores classes useful for storing annotation information for genes and their transcripts and exons. */
object GeneAnnotations {

  /** A gene with zero or more transcripts, that all are on the same transcription strand. */
  case class Gene(contig: String, start: Int, end: Int, negativeStrand: Boolean, name: String, transcripts: Seq[Transcript])
    extends Interval(contig, start, end, negativeStrand, name) with Iterable[Transcript] {
    def iterator: Iterator[Transcript] = transcripts.toIterator
  }

  /** A transcript associated with a given gene (which stores the transcription strand).  Contains zero or more exons.
    * The exons should be given in transcript order. */
  case class Transcript(name: String, start: Int, end: Int, cdsStart: Int, cdsEnd: Int, exons: Seq[Exon]) {
    if (exons.length > 1) {
      val exonsOverlap = genomicOrder.sliding(2).map { e => (e.head, e.last) }.exists { case (e1, e2) => CoordMath.overlaps(e1.start, e1.end, e2.start, e2.end) }
      require(!exonsOverlap, s"exons overlap for transcript: $name")
    }
    /** The order in which exons appear in the transcripts */
    def transcriptOrder: Iterator[Exon] = exons.toIterator
    /** The order in which exons appear in the genome */
    def genomicOrder: Iterator[Exon] = exons.sortBy { exon => (exon.start, exon.end) }.toIterator
  }

  /** Defines an exonic sequence within a transcript. */
  case class Exon(start: Int, end: Int) {
    require(start <= end, "start is greater than end when creating an exon")
  }
}