/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.umi

import java.nio.file.Files

import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.umi.GroupReadsByUmi.{AdjacencyUmiAssigner, EarlierReadComparator}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.{SAMRecord, SamPairUtil, SamReaderFactory}
import com.fulcrumgenomics.FgBioDef._
import SamRecordSetBuilder._

import scala.collection.JavaConversions._
import scala.collection.mutable.ListBuffer


/**
  * Tests for the tool that groups reads by position and UMI to attempt to identify
  * read pairs that arose from the same original molecule.
  */
class GroupReadsByUmiTest extends UnitSpec {
  // Returns a List of the element 't' repeated 'n' times
  private def n[T](t: T, n:Int =1): List[T] = List.tabulate(n)(x => t)

  {
    "IdentityUmiAssigner" should "group UMIs together that have the exact same sequence" in {
      val assigner = new GroupReadsByUmi.IdentityUmiAssigner
      val umis = Seq("AAAAAA", "AAAAAT", "ACGTAC", "ACGTAC", "AAAAAA")
      assigner.assign(umis) should contain theSameElementsAs Map(
        "AAAAAA" -> "AAAAAA",
        "AAAAAT" -> "AAAAAT",
        "ACGTAC" -> "ACGTAC",
        "AAAAAA" -> "AAAAAA"
      )
    }
  }

  {
    "SimpleErrorUmiAssigner" should "group UMIs together by mismatches" in {
      val assigner = new GroupReadsByUmi.SimpleErrorUmiAssigner(1)
      val umis = Seq("AAAAAA", "AAAATT", "AAAATA",
                     "TGCACC", "TGCACG",
                     "GGCGGC", "GGCGGC", "GGCGGC")

      assigner.assign(umis) should contain theSameElementsAs Map(
        "AAAAAA" -> "AAAAAA",
        "AAAATT" -> "AAAAAA",
        "AAAATA" -> "AAAAAA",
        "GGCGGC" -> "GGCGGC",
        "TGCACC" -> "TGCACC",
        "TGCACG" -> "TGCACC"
      )
    }

    it should "(stupidly) assign everything to the same tag" in {
      val assigner = new GroupReadsByUmi.SimpleErrorUmiAssigner(6)
      val umis = Seq("AAAAAA", "AAAATT", "AAAATA", "TGCACC", "TGCACG", "GGCGGC", "GGCGGC", "GGCGGC")

      assigner.assign(umis) should contain theSameElementsAs Map(
        "AAAAAA" -> "AAAAAA",
        "AAAATA" -> "AAAAAA",
        "AAAATT" -> "AAAAAA",
        "GGCGGC" -> "AAAAAA",
        "TGCACC" -> "AAAAAA",
        "TGCACG" -> "AAAAAA"
      )
    }

    {
      "AdjacencyUmiAssigner" should "assign each UMI to separate groups" in {
        val umis = Seq("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "AAATTT", "TTTAAA", "AGAGAG")
        val map  = new AdjacencyUmiAssigner(maxMismatches=2).assign(umis)
        map shouldEqual umis.map(u => u -> u).toMap
      }

      it should "assign everything into one group when all counts=1 and within mismatch threshold" in {
        val umis = Seq("AAAAAA", "AAAAAc", "AAAAAg").map(_.toUpperCase)
        val map  = new AdjacencyUmiAssigner(maxMismatches=1).assign(umis)
        map should have size 3
        map.values.toSet should have size 1
      }

      it should "assign everything into one group pointing at AAAAAA" in {
        val umis = Seq("AAAAAA", "AAAAAA", "AAAAAA", "AAAAAc", "AAAAAc", "AAAAAg", "AAAtAA").map(_.toUpperCase)
        val map  = new AdjacencyUmiAssigner(maxMismatches=1).assign(umis)
        map shouldEqual umis.map(umi => umi -> "AAAAAA").toMap
      }

      it should "make three groups" in {
        val umis: Seq[String] = n("AAAAAA", 4) ++ n("AAAAAT", 2) ++ n("AATAAT", 1) ++ n("AATAAA", 2) ++
                                n("GACGAC", 9) ++ n("GACGAT", 1) ++ n("GACGCC", 4) ++
                                n("TACGAC", 7)

        val map  = new AdjacencyUmiAssigner(maxMismatches=2).assign(umis)
        map shouldEqual umis.map(umi => {
          if (umi.startsWith("A"))      umi -> "AAAAAA"
          else if (umi.startsWith("G")) umi -> "GACGAC"
          else                          umi -> "TACGAC"
        }).toMap
      }

      // Unit test for something that failed when running on real data
      it should "correctly assign the following UMIs" in {
        val umis = Seq("CGGGGG", "GTGGGG", "GGGGGG", "CTCACA", "TGCAGT", "CTCACA", "CGGGGG")
        val map  = new AdjacencyUmiAssigner(maxMismatches=1).assign(umis)
        val exp  = Map("CGGGGG" -> "CGGGGG", // 2 obs
                       "GGGGGG" -> "CGGGGG", // 1 obs
                       "GTGGGG" -> "CGGGGG", // 1 obs
                       "CTCACA" -> "CTCACA", // 2 obs
                       "TGCAGT" -> "TGCAGT") // 1 obs
        map should contain theSameElementsAs exp
      }

      it should "handle a deep tree of UMIs" in {
        val umis = n("AAAAAA", 256) ++ n("TAAAAA", 128) ++ n("TTAAAA", 64) ++ n("TTTAAA", 32) ++ n("TTTTAA", 16)
        val map  = new AdjacencyUmiAssigner(maxMismatches=1).assign(umis)
        umis.foreach(umi => map should contain (umi -> "AAAAAA"))
      }
    }
  }

  // Tests for the comparator which is used internally
  {
    val comparator = new GroupReadsByUmi.EarlierReadComparator()

    "EarlierReadComparator" should "reject unmapped reads" in {
      val recs = new SamRecordSetBuilder(readLength=100).addPair("q1", start1=100, start2=100, record1Unmapped = true)
      an[AssertionError] should be thrownBy recs.sort(comparator)
    }

    it should "reject unpaired reads" in {
      val builder = new SamRecordSetBuilder(readLength=100)
      builder.addFrag("q1", start=100)
      builder.addFrag("q2", start=200)
      an[AssertionError] should be thrownBy builder.toList.sort(comparator)
    }

    it should "reject secondary reads" in {
      val recs = new SamRecordSetBuilder(readLength=100).addPair("q1", start1=100, start2=300)
      recs.foreach(rec => rec.setNotPrimaryAlignmentFlag(true))
      an[AssertionError] should be thrownBy recs.sort(comparator)
    }

    it should "reject supplementary reads" in {
      val recs = new SamRecordSetBuilder(readLength=100).addPair("q1", start1=100, start2=300)
      recs.foreach(rec => rec.setSupplementaryAlignmentFlag(true))
      an[AssertionError] should be thrownBy recs.sort(comparator)
    }

    it should "reject pairs with reads mapped to different chromosomes" in {
      val recs = new SamRecordSetBuilder(readLength=100).addPair("q1", start1=100, start2=300)
      recs match {
        case Seq(r1, r2) =>
          r2.setReferenceIndex(3)
          SamPairUtil.setMateInfo(r1, r2, true)
        case _           => unreachable("This should never happen!!")
      }

      an[AssertionError] should be thrownBy recs.sort(comparator)
    }

    it should "sort pairs by the 'lower' 5' position of the pair" in {
      val builder = new SamRecordSetBuilder(readLength=100, sortOrder=SortOrder.coordinate)
      val exp = ListBuffer[SAMRecord]()
      // Records are added to the builder in the order that we expect them to be sorted, but the builder
      // will coordinate sort them for us, so we can re-sort them and test the results
      exp ++= builder.addPair("q1", 0, 100, 300)
      exp ++= builder.addPair("q2", 0, 106, 300, cigar1="5S95M") // effective=101
      exp ++= builder.addPair("q3", 0, 102, 299)
      exp ++= builder.addPair("q4", 0, 300, 110, strand1=Minus, strand2=Plus)
      exp ++= builder.addPair("q5", 0, 120, 320)
      exp ++= builder.addPair("q6", 1, 1, 200)

      val f = (r:SAMRecord) => r.getReadName + "/" + (if (r.getFirstOfPairFlag) 1 else 2 )
      val expected = exp.map(f)
      val actual   = builder.toList.sorted(Ordering.comparatorToOrdering(comparator)).map(f)

      actual should contain theSameElementsInOrderAs expected
    }
  }


  // Test for running the GroupReadsByUmi command line program with some sample input
  "GroupReadsByUmi" should "group reads correctly" in {
    val builder = new SamRecordSetBuilder(readLength=100, sortOrder=SortOrder.coordinate)
    builder.addPair(name="a01", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addPair(name="a02", start1=100, start2=300, attrs=Map("RX" -> "AAAAgAAA"))
    builder.addPair(name="a03", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addPair(name="a04", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAt"))
    builder.addPair(name="a05", start1=100, start2=300, record2Unmapped=true, attrs=Map("RX" -> "AAAAAAAt"))
    builder.addPair(name="a06", start1=100, start2=300, mapq1=5)

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", strategy="edit", edits=1).execute()
    println(out.toString)

    val reader = SamReaderFactory.make().open(out.toFile)
    val groups = reader.iterator().toSeq.groupBy(_.getReadName.charAt(0))

    // Group A: Read 5 out for unmapped mate, 6 out for low mapq, 1-4 all passed through into one umi group
    groups('a') should have size 4*2
    groups('a').map(_.getReadName).toSet shouldEqual Set("a01", "a02", "a03", "a04")
    groups('a').map(_.getAttribute("MI")).toSet should have size 1

    reader.close()
  }
}
