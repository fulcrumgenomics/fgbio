package com.fulcrumgenomics.personal.yfarjoun

import htsjdk.samtools.util.{Interval, Locatable}

import scala.collection.mutable

object LocusTrack {
  val Empty = new LocusTrack(None)

  class ConnectedSlice[A](val mySeq: mutable.Seq[A], val from: Int, val to: Int) extends mutable.Seq[A] {

    assert(from >= 0)
    assert(to >= from)
    assert(to <= mySeq.size)

    override def update(idx: Int, elem: A): Unit = mySeq.update(idx + from, elem)

    override def apply(i: Int): A = mySeq.apply(i + from)

    override def length: Int = to - from

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
      case None => LocusTrack.Empty
      case Some(loc) if !loc.contigsMatch(l) => LocusTrack.Empty
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
