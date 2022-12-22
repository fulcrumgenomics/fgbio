package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.{Interval, IntervalList}

import scala.jdk.CollectionConverters.SeqHasAsJava

class CoverageManagerTest extends UnitSpec {

  private val dict = SequenceDictionary(
    SequenceMetadata(name = "chr1", length = 10000),
    SequenceMetadata(name = "chr2", length = 50000)
  ).asSam

  private val intervals = new IntervalList(this.dict)

  val i1: Interval = new Interval("chr1", 100, 200)
  val i2: Interval = new Interval("chr1", 100, 300)
  val i3: Interval = new Interval("chr1", 300, 400)

  val i4: Interval = new Interval("chr2", 300, 400)

  val il: IntervalList = new IntervalList(dict)
  il.addall(List(i1, i2, i3).asJava)
  val covMan: CoverageManager = new CoverageManager(il)

  "CoverageManager" should "start its life with coverage = zero" in {
    covMan.getMinCoverage(i1) shouldEqual 0
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
  }

  "CoverageManager.getMinCoverage" should "return Short.MaxValue when out of range" in {
    covMan.getMinCoverage(i4) shouldEqual Short.MaxValue
  }

  "CoverageManager" should "remain coverage zero when only partially covered" in {
    covMan.resetCoverage()
    covMan.incrementCoverage(i1)
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
    covMan.getMinCoverage(i1) shouldEqual 1
  }

  "CoverageManager.getMinCoverage" should "remain Short.MaxValue when out of range" in {
    covMan.getMinCoverage(i4) shouldEqual Short.MaxValue
  }


  "CoverageManager" should "get get to 2 when adding coverage twice" in {
    covMan.resetCoverage()
    covMan.incrementCoverage(i1)
    covMan.incrementCoverage(i1)
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
    covMan.getMinCoverage(i1) shouldEqual 2
  }


}
