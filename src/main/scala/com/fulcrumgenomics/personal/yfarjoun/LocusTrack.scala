package com.fulcrumgenomics.personal.yfarjoun

import htsjdk.samtools.util.{Interval, Locatable}
import scala.collection.mutable

object LocusTrack {
  val Empty = new LocusTrack(None)


  /**
    * A small class which implements a "connected slice", i.e. a slice of a mutable.Seq which, when modified
    * simply modifies its "parent" Seq. Unlike Seq.slice and Array.slice, copies of data are not made.
    * @param mySeq original Seq[A] of which this is a slice
    * @param from index from which slice starts (0-based, inclusive.) Must be >=0
    * @param to index to which slice continues (0-based, exclusive). Must be >= from
    * @tparam A Type of Seq
    */
  class ConnectedSlice[A](val mySeq: mutable.Seq[A], val from: Int, val to: Int) extends mutable.Seq[A] {

    assert(from >= 0)
    assert(to >= from)
    assert(to <= mySeq.size)

    /**
      * Updates element idx in slice to value elem. Will (only) update appropriate value in origSeq
      * @param idx index of element in the slice to update
      * @param elem new value of element
      */
    @throws[IndexOutOfBoundsException]
    override def update(idx: Int, elem: A): Unit = {
      if (idx < 0 || idx >= length) {
        throw new IndexOutOfBoundsException(idx)
      }
      mySeq.update(idx + from, elem)
    }

    /** Get the element at the specified index. */
    @throws[IndexOutOfBoundsException]
    override def apply(i: Int): A = {
      if (i < 0 || i >= length) {
        throw new IndexOutOfBoundsException(i)
      }
      mySeq.apply(i + from)
    }

    override val length: Int = to - from

    /** Provides an iterator over the elements between from and to in mySeq.
      * */
    override def iterator: Iterator[A] = mySeq.slice(from, to).iterator
  }
}

class LocusTrack(val locus: Option[Interval], val track: mutable.Seq[Short]) {
  def this(l: Option[Interval]) = {
    this(l, l match {
      case None => mutable.Seq.empty
      case Some(loc) => new Array[Short](loc.length())
    })
  }

  def this(l: Interval) = {
    this(Some(l))
  }

  // return a slice of the track which matches a region overlapping a
  // given locus
  def sliceToLocus(l: Locatable): LocusTrack = {
    locus match {
      case None =>
        LocusTrack.Empty
      case Some(loc) if !loc.contigsMatch(l) =>
        LocusTrack.Empty
      case Some(loc) =>
        val start: Int = Math.max(loc.getStart, l.getStart)
        val end: Int = Math.min(loc.getEnd, l.getEnd)
        if (start > end) {
          LocusTrack.Empty
        } else {
          val overlapLocus = new Interval(loc.getContig, start, end)
          val overlapTrack = new LocusTrack.ConnectedSlice(track, start - loc.getStart, end + 1 - loc.getStart)
          new LocusTrack(Some(overlapLocus), overlapTrack)
        }
    }
  }
}
