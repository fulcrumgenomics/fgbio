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

    def apply(left: Locatable) = {
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

  var union_1_2: Locatable = NormalizeCoverage.union(i1, i2)
  var union_2_1: Locatable = NormalizeCoverage.union(i2, i1)
  var union_1_3: Locatable = NormalizeCoverage.union(i1, i3)
  var union_3_1: Locatable = NormalizeCoverage.union(i3, i1)
  var union_2_3: Locatable = NormalizeCoverage.union(i2, i3)
  var union_3_2: Locatable = NormalizeCoverage.union(i3, i2)

  "Union" should " correctly find the union of two locatables" in {

    union_1_2 should haveSameLocusAs(new Interval("chr1",100,300))
    union_1_3 should haveSameLocusAs(new Interval("chr1",100,400))
    union_2_3 should haveSameLocusAs(new Interval("chr1",100,400))
  }

  it should "be symmetric" in {
    union_1_2 should haveSameLocusAs(union_2_1)
    union_1_3 should haveSameLocusAs(union_3_1)
    union_3_2 should haveSameLocusAs(union_2_3)
  }

  var subset_1_2: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i1, Set(i2))
  var subset_2_1: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i2, Set(i1))
  var subset_1_3: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i1, Set(i3))
  var subset_3_1: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i3, Set(i1))
  var subset_2_3: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i2, Set(i3))
  var subset_3_2: Option[Locatable] = NormalizeCoverage.subsetReadToLocus(i3, Set(i2))

  "subsetReadToLocus" should "correctly subset one locatable to another" in {
    subset_1_2 should not be empty
    subset_1_2.get should haveSameLocusAs(i1)

    subset_2_1 should not be empty
    subset_2_1.get should haveSameLocusAs(i1)

    subset_1_3 shouldBe empty
    subset_3_1 shouldBe empty

    subset_2_3 should not be empty
    subset_2_3.get should haveSameLocusAs(new Interval("chr1",300,300))

    subset_3_2 should not be empty
    subset_3_2.get should haveSameLocusAs(new Interval("chr1",300,300))
  }
}
