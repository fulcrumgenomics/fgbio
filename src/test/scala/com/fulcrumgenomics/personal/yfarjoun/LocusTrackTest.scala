package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.personal.yfarjoun.LocusTrack.Empty
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.Interval

class LocusTrackTest extends UnitSpec {
  val i1: Interval = new Interval("chr1", 100, 200)
  val i2: Interval = new Interval("chr1", 150, 300)
  val i3: Interval = new Interval("chr1", 300, 400)

  val i4: Interval = new Interval("chr2", 300, 400)

  val lt1 = new LocusTrack(i1)
  val lt2 = new LocusTrack(i2)
  val lt3 = new LocusTrack(i3)
  val lt4 = new LocusTrack(i4)
  val ltEmpty = new LocusTrack(None)

  "LocusTrack" should "Start off with a seq[Short] of the correct length" in {
    lt1.track should have length i1.length()
    lt2.track should have length i2.length()
    lt3.track should have length i3.length()
    lt4.track should have length i4.length()
  }

  it should "start off containing only zeros" in {
    lt1.track.count(_ == 0) shouldEqual i1.length()
    lt2.track.count(_ == 0) shouldEqual i2.length()
    lt3.track.count(_ == 0) shouldEqual i3.length()
    lt4.track.count(_ == 0) shouldEqual i4.length()
  }

  it should "return the emtpy locusTrack when its empty" in{
    ltEmpty.sliceToLocus(i1) shouldBe Empty
    lt3.sliceToLocus(i1) shouldBe Empty
    lt4.sliceToLocus(i1) shouldBe Empty
  }

  it should "get incremented when incrementing a slice" in {
    val lt2_slice_i1 = lt2.sliceToLocus(i1)
    // this increments the overlap of i1 and i2, so chr1:50-100
    lt2_slice_i1.track.forall(_==0) shouldBe true
    lt2_slice_i1.track.indices.foreach(i => lt2_slice_i1.track(i) = (lt2_slice_i1.track(i) + 1).toShort)
    lt2_slice_i1.track.forall(_==1) shouldBe true

    val map: Map[Short, Int] = lt2.track
      .indices
      .groupBy(i => lt2.track(i))
      .map { case (k: Short, v: Seq[Short]) => (k, v.length) }
    // 150->200 should be incremented, and 201->300 not.
    map shouldEqual Map(0 -> 100, 1 -> 51)
  }


  it should "not fall over when there's no overlap" in {
    val lt3_slice_i1 = lt3.sliceToLocus(i1)
    // this increments the overlap of i1 and i2, so chr1:50-100
    lt3_slice_i1 shouldBe LocusTrack.Empty
    lt3_slice_i1.track.indices
      .foreach(i => lt3_slice_i1.track(i) = (lt3_slice_i1.track(i) + 1).toShort)
    lt3_slice_i1.track.forall(_==1) shouldBe true

    val map: Map[Short, Int] = lt3.track.indices.groupBy(i => lt3.track(i)).map { case (k: Short, v: Seq[Short]) => (k, v.length) }
    // 150->200 should be incremented, and 201->300 not.
    map shouldEqual Map(0 -> 101)
  }
}
