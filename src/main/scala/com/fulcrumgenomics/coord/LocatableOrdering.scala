package com.fulcrumgenomics.coord

import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util.Locatable

/** Methods for building orderings of [[Locatable]] instances. */
object LocatableOrdering {

  /** Build a coordinate-based ordering of [[Locatable]] instances. */
  def apply(dict: SequenceDictionary): Ordering[Locatable] = (x: Locatable, y: Locatable) => {
    var compare = (dict.get(x.getContig), dict.get(y.getContig)) match {
      case (Some(meta1), Some(meta2)) => meta1.index.compare(meta2.index)
      case (None, _) => throw new NoSuchElementException(s"Sequence dictionary does not contain contig: ${x.getContig}")
      case (_, None) => throw new NoSuchElementException(s"Sequence dictionary does not contain contig: ${y.getContig}")
    }
    if (compare == 0) compare = x.getStart.compare(y.getStart)
    if (compare == 0) compare = x.getEnd.compare(y.getEnd)
    compare
  }
}
