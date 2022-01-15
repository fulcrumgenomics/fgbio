package com.fulcrumgenomics.coord

import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util.Locatable

/** Methods for working with classes that implement [[Locatable]] */
object LocatableUtil {

  /** Coordinate ordering for [[Locatable]] instances. */
  case class GenomicOrdering(dict: SequenceDictionary) extends Ordering[Locatable] {

    /** Compare two locatables on a genomic coordinate system that is contig-aware. */
    override def compare(x: Locatable, y: Locatable): Int = {
      require(dict.contains(x.getContig), s"Sequence dictionary does not contain contig: $x")
      require(dict.contains(y.getContig), s"Sequence dictionary does not contain contig: $y")
      var compare = dict(x.getContig).index.compare(dict(y.getContig).index)
      if (compare == 0) compare = x.getStart.compare(y.getStart)
      if (compare == 0) compare = x.getEnd.compare(y.getEnd)
      compare
    }
  }
}
