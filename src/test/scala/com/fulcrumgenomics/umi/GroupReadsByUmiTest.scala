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

import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.GroupReadsByUmi._

import scala.collection.mutable

/**
  * Tests for the tool that groups reads by position and UMI to attempt to identify
  * read pairs that arose from the same original molecule.
  */
class GroupReadsByUmiTest extends UnitSpec {
  // Returns a List of the element 't' repeated 'n' times
  private def n[T](t: T, n:Int =1): List[T] = List.tabulate(n)(x => t)

  /**
    * Converts a mapping of umi->id to a Set of Sets of Umis that are assigned the same ID.
    * E.g. { "AAA" -> 1, "AAC" -> 1, "CCC" -> 2 } => (("AAA", "AAC"), ("CCC"))
    */
  private def group(assignments: Map[Umi,MoleculeId], stripSuffix: Boolean = false): Set[Set[Umi]] = {
    def strip(s:String) = { if (stripSuffix) s.substring(0, s.indexOf('/')) else s }

    val groups: Map[MoleculeId, mutable.Set[Umi]] = assignments.map(kv => (strip(kv._2), mutable.Set[Umi]()))
    assignments.foreach { case (umi, id) => groups(strip(id)).add(umi) }
    groups.values.map(_.toSet).toSet
  }


  {
    "IdentityUmiAssigner" should "group UMIs together that have the exact same sequence" in {
      val assigner = new GroupReadsByUmi.IdentityUmiAssigner
      val umis = Seq("AAAAAA", "AAAAAT", "ACGTAC", "ACGTAC", "AAAAAA")
      group(assigner.assign(umis)) should contain theSameElementsAs Set(Set("AAAAAA"), Set("AAAAAT"), Set("ACGTAC"))
    }
  }

  {
    "SimpleErrorUmiAssigner" should "group UMIs together by mismatches" in {
      val assigner = new GroupReadsByUmi.SimpleErrorUmiAssigner(1)
      val umis = Seq("AAAAAA", "AAAATT", "AAAATA",
                     "TGCACC", "TGCACG",
                     "GGCGGC", "GGCGGC", "GGCGGC")

      group(assigner.assign(umis)) should contain theSameElementsAs Set(Set("AAAAAA", "AAAATT", "AAAATA"), Set("GGCGGC"), Set("TGCACC", "TGCACG"))
    }

    it should "(stupidly) assign everything to the same tag" in {
      val assigner = new GroupReadsByUmi.SimpleErrorUmiAssigner(6)
      val umis = Seq("AAAAAA", "AAAATT", "AAAATA", "TGCACC", "TGCACG", "GGCGGC", "GGCGGC", "GGCGGC")

      group(assigner.assign(umis)) should contain theSameElementsAs Set(Set("AAAAAA", "AAAATA", "AAAATT", "GGCGGC", "TGCACC", "TGCACG"))
    }

    {
      "AdjacencyUmiAssigner" should "assign each UMI to separate groups" in {
        val umis = Seq("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "AAATTT", "TTTAAA", "AGAGAG")
        val groups  = group(new AdjacencyUmiAssigner(maxMismatches=2).assign(umis))
        groups shouldBe umis.map(Set(_)).toSet
      }

      it should "assign everything into one group when all counts=1 and within mismatch threshold" in {
        val umis = Seq("AAAAAA", "AAAAAc", "AAAAAg").map(_.toUpperCase)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1).assign(umis))
        groups shouldBe Set(umis.toSet)
      }

      it should "assign everything into one group" in {
        val umis = Seq("AAAAAA", "AAAAAA", "AAAAAA", "AAAAAc", "AAAAAc", "AAAAAg", "AAAtAA").map(_.toUpperCase)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1).assign(umis))
        groups shouldBe Set(umis.toSet)
      }

      it should "make three groups" in {
        val umis: Seq[String] = n("AAAAAA", 4) ++ n("AAAAAT", 2) ++ n("AATAAT", 1) ++ n("AATAAA", 2) ++
                                n("GACGAC", 9) ++ n("GACGAT", 1) ++ n("GACGCC", 4) ++
                                n("TACGAC", 7)

        val groups  = group(new AdjacencyUmiAssigner(maxMismatches=2).assign(umis))
        groups shouldBe Set(
          Set("AAAAAA", "AAAAAT", "AATAAT", "AATAAA"),
          Set("GACGAC", "GACGAT", "GACGCC"),
          Set("TACGAC")
        )
      }

      // Unit test for something that failed when running on real data
      it should "correctly assign the following UMIs" in {
        val umis   = Seq("CGGGGG", "GTGGGG", "GGGGGG", "CTCACA", "TGCAGT", "CTCACA", "CGGGGG")
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1).assign(umis))
        groups shouldBe Set(Set("CGGGGG", "GGGGGG", "GTGGGG"), Set("CTCACA"), Set("TGCAGT"))
      }

      it should "handle a deep tree of UMIs" in {
        val umis   = n("AAAAAA", 256) ++ n("TAAAAA", 128) ++ n("TTAAAA", 64) ++ n("TTTAAA", 32) ++ n("TTTTAA", 16)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1).assign(umis))
        groups shouldBe Set(Set("AAAAAA", "TAAAAA", "TTAAAA", "TTTAAA", "TTTTAA"))
      }
    }

    {
      "PairedUmiAssigner" should "assign A-B and B-A into groups with the same prefix but different suffix" in {
        val umis = Seq("AAAA-CCCC", "CCCC-AAAA")
        val map = new PairedUmiAssigner(maxMismatches=1).assign(umis)
        group(map, stripSuffix=false) shouldBe Set(Set("AAAA-CCCC"), Set("CCCC-AAAA"))
        group(map, stripSuffix=true)  shouldBe Set(Set("AAAA-CCCC", "CCCC-AAAA"))
      }

      it should "assign A-B and B-A into groups with the same prefix but different suffix, with errors" in {
        val umis = Seq("AAAA-CCCC", "CCCC-CAAA", "AAAA-GGGG", "GGGG-AAGA")
        val map = new PairedUmiAssigner(maxMismatches=1).assign(umis)
        group(map, stripSuffix=false) shouldBe Set(Set("AAAA-CCCC"), Set("CCCC-CAAA"), Set("AAAA-GGGG"), Set("GGGG-AAGA"))
        group(map, stripSuffix=true)  shouldBe Set(Set("AAAA-CCCC", "CCCC-CAAA"), Set("AAAA-GGGG", "GGGG-AAGA"))
      }

      it should "handle errors in the first base changing lexical ordering of AB vs. BA" in {
        val umis   = n("GTGT-ACAC", 500) ++ n("ACAC-GTGT", 460) ++ n("GTGT-TCAC", 6)  ++ n("TCAC-GTGT", 6) ++ n("GTGT-TGAC", 1)

        val map = new PairedUmiAssigner(maxMismatches=1).assign(umis)
        group(map, stripSuffix=false) shouldBe Set(Set("GTGT-ACAC", "GTGT-TCAC", "GTGT-TGAC"), Set("TCAC-GTGT", "ACAC-GTGT"))
        group(map, stripSuffix=true)  shouldBe Set(Set("GTGT-ACAC", "ACAC-GTGT", "GTGT-TCAC", "TCAC-GTGT", "GTGT-TGAC"))
      }

      it should "count A-B and B-A together when constructing adjacency graph" in {
        // Since the graph only creates nodes where count(child) <= count(parent) / 2 + 1, it should
        // group everything together in the first set, but not in the second set.
        val umis1   = n("AAAA-CCCC", 256) ++ n("AAAA-CCCG", 64)  ++ n("CCCG-AAAA", 64)
        val umis2   = n("AAAA-CCCC", 256) ++ n("AAAA-CCCG", 128) ++ n("CCCG-AAAA", 128)

        val map1 = new PairedUmiAssigner(maxMismatches=1).assign(umis1)
        val map2 = new PairedUmiAssigner(maxMismatches=1).assign(umis2)
        group(map1, stripSuffix=true) shouldBe Set(Set("AAAA-CCCC", "AAAA-CCCG", "CCCG-AAAA"))
        group(map2, stripSuffix=true) shouldBe Set(Set("AAAA-CCCC"), Set("AAAA-CCCG", "CCCG-AAAA"))
      }

      it should "fail if supplied non-paired UMIs" in {
        val umis   = Seq("AAAAAAAA", "GGGGGGGG")
        an[IllegalStateException] shouldBe thrownBy { new PairedUmiAssigner(maxMismatches=1).assign(umis) }
      }
    }
  }

  // Test for running the GroupReadsByUmi command line program with some sample input
  "GroupReadsByUmi" should "group reads correctly" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addPair(name="a02", start1=100, start2=300, attrs=Map("RX" -> "AAAAgAAA"))
    builder.addPair(name="a03", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addPair(name="a04", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAt"))
    builder.addPair(name="a05", start1=100, start2=300, unmapped2=true, attrs=Map("RX" -> "AAAAAAAt"))
    builder.addPair(name="a06", start1=100, start2=300, mapq1=5)

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), rawTag="RX", assignTag="MI", strategy=Strategy.Edit, edits=1).execute()

    val groups = readBamRecs(out).groupBy(_.name.charAt(0))

    // Group A: Read 5 out for unmapped mate, 6 out for low mapq, 1-4 all passed through into one umi group
    groups('a') should have size 4*2
    groups('a').map(_.name).toSet shouldEqual Set("a01", "a02", "a03", "a04")
    groups('a').map(r => r[String]("MI")).toSet should have size 1

    hist.toFile.exists() shouldBe true
  }

  it should "correctly group reads with the paired assigner when the two UMIs are the same" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a03", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a04", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b01", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b02", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b03", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b04", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()

    val recs = readBamRecs(out)

    val aIds = recs.filter(_.name.startsWith("a")).map(r => r[String]("MI")).distinct
    val bIds = recs.filter(_.name.startsWith("b")).map(r => r[String]("MI")).distinct

    aIds should have size 1
    bIds should have size 1

    aIds.head.takeWhile(_ != '/') shouldBe bIds.head.takeWhile(_ != '/')
    aIds.head should not equal bIds.head
  }

  it should "exclude reads that contain an N in the UMI" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a03", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ANN"))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".bam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=2).execute()

    readBamRecs(out).map(_.name).distinct shouldBe Seq("a01", "a02")
  }

  it should "fail when umis have different length" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-AC"))

    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    val tool = new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1)

    an[Exception] should be thrownBy tool.execute()
  }

  it should "fail when the raw tag is not present" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus, strand2=Minus)
    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    val exception = intercept[FailureException] {
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()
    }
    exception.getMessage should include ("was missing the raw UMI tag")
  }

  Strategy.values.filterNot(_ == Strategy.Paired).foreach { strategy =>
    it should s"reject records with UMIs that are shorter than the specified minimum length with the $strategy strategy" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))

      builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACT"))
      builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTAC"))

      val in   = builder.toTempFile()
      val out  = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Identity, edits=0, minUmiLength=Some(6)).execute()

      val recs = readBamRecs(out)
      recs should have length 2
      recs.map(_.name) should contain theSameElementsInOrderAs Seq("a01", "a01")
      recs.map(r => r[String]("MI")).distinct should have length 1
    }

    it should s"truncate to the specified minimum length with the $strategy strategy" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
      builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACT"))
      builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTAC"))

      val in   = builder.toTempFile()
      val out  = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Identity, edits=0, minUmiLength=Some(5)).execute()

      val recs = readBamRecs(out)
      recs should have length 4
      recs.map(_.name) should contain theSameElementsInOrderAs Seq("a01", "a01", "a02", "a02")
      recs.map(r => r[String]("MI")).distinct should have length 1 // all should be assigned to one molecule
    }
  }
}
