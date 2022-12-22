package com.fulcrumgenomics.personal.yfarjoun

//import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.util.Locatable
//import com.fulcrumgenomics.commons.CommonsDef.PathToBam
import htsjdk.samtools.util.Interval

import org.scalatest._
import matchers._

trait LocatableMatchers {

  class SameLocusAs(right: Locatable) extends Matcher[Locatable] {

    def apply(left: Locatable): MatchResult = {
      MatchResult(
        left.getContig.equals(right.getContig) &&
          left.getStart.equals(right.getStart) &&
          left.getEnd.equals(right.getEnd),
        s"""Locatable $left is not the same as $right"""",
        s"""Locatable $left is the same as $right""""
      )
    }
  }

  def haveSameLocusAs(right: Locatable) = new SameLocusAs(right)
}

object LocatableMatchers extends LocatableMatchers

class NormalizeCoverageOptionsTest extends UnitSpec with LocatableMatchers {

//  val nc = new NormalizeCoverage(PathToBam("/dev/null"), Path(""), Path(""))

  val i1: Interval = new Interval("chr1", 100, 200)
  val i2: Interval = new Interval("chr1", 100, 300)
  val i3: Interval = new Interval("chr1", 300, 400)


//
//  var subset_1_2: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i1), Seq(i2))
//  var subset_2_1: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i2), Seq(i1))
//  var subset_1_3: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i1), Seq(i3))
//  var subset_3_1: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i3), Seq(i1))
//  var subset_2_3: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i2), Seq(i3))
//  var subset_3_2: Seq[Locatable] = CoverageManager.subsetReadToLocus(Seq(i3), Seq(i2))
//
//  "subsetReadToLocus" should "correctly subset one locatable to another" in {
//    subset_1_2 should have size 1
//    subset_1_2.head should haveSameLocusAs(i1)
//
//    subset_2_1 should have size 1
//    subset_2_1.head should haveSameLocusAs(i1)
//
//    subset_1_3 shouldBe empty
//    subset_3_1 shouldBe empty
//
//    subset_2_3 should have size 1
//    subset_2_3.head should haveSameLocusAs(new Interval("chr1",300,300))
//
//    subset_3_2 should have size 1
//    subset_3_2.head should haveSameLocusAs(new Interval("chr1",300,300))
//  }
}
