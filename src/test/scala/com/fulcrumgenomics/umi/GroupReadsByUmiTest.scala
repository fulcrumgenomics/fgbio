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

import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.bam.api.SamOrder.TemplateCoordinate
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.GroupReadsByUmi._
import com.fulcrumgenomics.util.Metric
import org.scalatest.{OptionValues, PrivateMethodTester}

import java.nio.file.Files
import scala.collection.mutable

/**
  * Tests for the tool that groups reads by position and UMI to attempt to identify
  * read pairs that arose from the same original molecule.
  */
class GroupReadsByUmiTest extends UnitSpec with OptionValues with PrivateMethodTester {
  // Returns a List of the element 't' repeated 'n' times
  private def n[T](t: T, n: Int): List[T] = List.tabulate(n)(_ => t)

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

    Seq(1, 4).foreach { threads =>
      "AdjacencyUmiAssigner" should s"assign each UMI to separate groups with $threads thread(s)" in {
        val umis = Seq("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "AAATTT", "TTTAAA", "AGAGAG")
        val groups  = group(new AdjacencyUmiAssigner(maxMismatches=2, threads=threads).assign(umis))
        groups shouldBe umis.map(Set(_)).toSet
      }

      it should f"assign everything into one group when all counts=1 and within mismatch threshold with $threads thread(s)" in {
        val umis = Seq("AAAAAA", "AAAAAc", "AAAAAg").map(_.toUpperCase)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))
        groups shouldBe Set(umis.toSet)
      }

      it should f"assign everything into one group with $threads thread(s)" in {
        val umis = Seq("AAAAAA", "AAAAAA", "AAAAAA", "AAAAAc", "AAAAAc", "AAAAAg", "AAAtAA").map(_.toUpperCase)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))
        groups shouldBe Set(umis.toSet)
      }

      it should f"make three groups with $threads thread(s)" in {
        val umis: Seq[String] = n("AAAAAA", 4) ++ n("AAAAAT", 2) ++ n("AATAAT", 1) ++ n("AATAAA", 2) ++
                                n("GACGAC", 9) ++ n("GACGAT", 1) ++ n("GACGCC", 4) ++
                                n("TACGAC", 7)

        val groups  = group(new AdjacencyUmiAssigner(maxMismatches=2, threads=threads).assign(umis))
        groups shouldBe Set(
          Set("AAAAAA", "AAAAAT", "AATAAT", "AATAAA"),
          Set("GACGAC", "GACGAT", "GACGCC"),
          Set("TACGAC")
        )
      }

      // Unit test for something that failed when running on real data
      it should f"correctly assign the following UMIs with $threads thread(s)" in {
        val umis   = Seq("CGGGGG", "GTGGGG", "GGGGGG", "CTCACA", "TGCAGT", "CTCACA", "CGGGGG")
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))
        groups shouldBe Set(Set("CGGGGG", "GGGGGG", "GTGGGG"), Set("CTCACA"), Set("TGCAGT"))
      }

      it should f"handle a deep tree of UMIs with $threads thread(s)" in {
        val umis   = n("AAAAAA", 256) ++ n("TAAAAA", 128) ++ n("TTAAAA", 64) ++ n("TTTAAA", 32) ++ n("TTTTAA", 16)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))
        groups shouldBe Set(Set("AAAAAA", "TAAAAA", "TTAAAA", "TTTAAA", "TTTTAA"))
      }

      // Regression test for deterministic tie-breaking when UMIs have equal counts.
      // Without lexicographic tie-breaking in the sort, the order of equal-count UMIs
      // depends on HashMap iteration order, causing non-deterministic groupings.
      it should f"produce deterministic groupings when UMIs have equal counts with $threads thread(s)" in {
        // Setup: AAAAAA and AAAGAC both have count=2, AAAAAC has count=1
        // Both AAAAAA and AAAGAC are within 1 edit of AAAAAC (can capture it)
        // AAAGAC is 2 edits from AAAAAA (positions 3 and 5), so they won't be grouped together
        // With deterministic sorting, AAAAAA (lexicographically first) should capture AAAAAC
        val umis = n("AAAAAA", 2) ++ n("AAAGAC", 2) ++ n("AAAAAC", 1)
        val groups = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))

        // AAAAAA should capture AAAAAC because AAAAAA < AAAGAC lexicographically
        // AAAGAC should be in its own group
        groups shouldBe Set(Set("AAAAAA", "AAAAAC"), Set("AAAGAC"))

        // Verify determinism: run again and ensure same result
        val groups2 = group(new AdjacencyUmiAssigner(maxMismatches=1, threads=threads).assign(umis))
        groups2 shouldBe groups
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

  "GroupReadsByUmi.umiForRead" should "correctly assign a/b for paired UMI prefixes" in {
    val tool = new GroupReadsByUmi(rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits = 0)
    val builder = new SamBuilder(readLength=100)
    val templates = Seq(
      // These 4 should be a::AAA-b::TTT since contig1 is lower
      builder.addPair(contig = 1, contig2 = Some(2), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 1, contig2 = Some(2), start1=300, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 1, contig2 = Some(2), start1=100, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 1, contig2 = Some(2), start1=100, start2=100, strand1=Minus,  strand2=Plus, attrs=Map("RX" -> "AAA-TTT")),
      // These 4 should be b::AAA-a::TTT since contig2 is lower
      builder.addPair(contig = 2, contig2 = Some(1), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 2, contig2 = Some(1), start1=300, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 2, contig2 = Some(1), start1=100, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      builder.addPair(contig = 2, contig2 = Some(1), start1=100, start2=100, strand1=Minus,  strand2=Plus, attrs=Map("RX" -> "AAA-TTT")),
      // Should be a::AAA-b::TTT since r1 pos < r2 pos
      builder.addPair(contig = 1, start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      // Should be b::AAA-a::TTT since r1 pos < r2 pos
      builder.addPair(contig = 1, start1=300, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      // Should be a::AAA-b::TTT since same contig/pos, r1 positive strand
      builder.addPair(contig = 1, start1=100, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAA-TTT")),
      // Should be b::AAA-a::TTT since same contig/pos, r2 positive strand
      builder.addPair(contig = 1, start1=100, start2=100, strand1=Minus,  strand2=Plus, attrs=Map("RX" -> "AAA-TTT")),
    ).map { pair => Template(r1 = pair.headOption, r2 = pair.lastOption) }
    val expected = n("a::AAA-b::TTT", 4) ++ n("b::AAA-a::TTT", 4) ++ n("a::AAA-b::TTT", 1) ++ n("b::AAA-a::TTT", 1) ++ n("a::AAA-b::TTT", 1) ++ n("b::AAA-a::TTT", 1)
    val umiForRead = PrivateMethod[String](Symbol("umiForRead"))
    templates.map(t => tool invokePrivate umiForRead(t)) should contain theSameElementsInOrderAs expected
  }

  it should "correctly assign a/b for paired UMI prefixes when the UMI for one end of the source molecule is absent" in {
    val tool = new GroupReadsByUmi(rawTag = "RX", assignTag = "MI", strategy = Strategy.Paired, edits = 0)
    val builder = new SamBuilder(readLength = 100)
    val templates = Seq(
      // This should be a::AAA-b:: since contig1 is lower
      builder.addPair(contig = 1, contig2 = Some(2), start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "AAA-")),
      // This should be b::AAA-a:: since contig2 is lower
      builder.addPair(contig = 2, contig2 = Some(1), start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "AAA-")),
      // This should be a::-b::TTT since contig1 is lower
      builder.addPair(contig = 1, contig2 = Some(2), start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "-TTT")),
      // This should be b::-a::TTT since contig2 is lower
      builder.addPair(contig = 2, contig2 = Some(1), start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "-TTT")),
    ).map { pair => Template(r1 = pair.headOption, r2 = pair.lastOption) }
    val expected = n("a::AAA-b::", 1) ++ n("b::AAA-a::", 1) ++  n("a::-b::TTT", 1) ++ n("b::-a::TTT", 1)
    val umiForRead = PrivateMethod[String](Symbol("umiForRead"))
    templates.map(t => tool invokePrivate umiForRead(t)) should contain theSameElementsInOrderAs expected
  }

  "GroupReads.ReadInfo" should "extract the same ReadEnds from a Template as from an R1 with mate cigar" in {
    val specialCellTag = "XX"
    val builder = new SamBuilder(readLength=100, sort=None)
    val Seq(r1, r2) = builder.addPair(contig=2, contig2=Some(1), start1=300, start2=400, cigar1="10S90M", cigar2="90M10S")
    val template    = Template(r1=Some(r1), r2=Some(r2))

    r1(specialCellTag) = "GTTTA"

    val tReadInfo = ReadInfo(template, cellTag = specialCellTag)
    val rReadInfo = ReadInfo(r1, cellTag = specialCellTag)
    rReadInfo shouldBe tReadInfo

    tReadInfo.refIndex1 shouldBe 1
    tReadInfo.start1    shouldBe 499
    tReadInfo.refIndex2 shouldBe 2
    tReadInfo.start2    shouldBe 290
    tReadInfo.cellBarcode.value shouldBe "GTTTA"
  }

  // Test for running the GroupReadsByUmi command line program with some sample input
  "GroupReadsByUmi" should "group reads correctly" in {
    Seq(SamOrder.Coordinate, SamOrder.Queryname, SamOrder.TemplateCoordinate).foreach { sortOrder =>
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      builder.addPair(name="a01", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
      builder.addPair(name="a02", start1=100, start2=300, attrs=Map("RX" -> "AAAAgAAA"))
      builder.addPair(name="a03", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAA"))
      builder.addPair(name="a04", start1=100, start2=300, attrs=Map("RX" -> "AAAAAAAt"))
      builder.addPair(name="b01", start1=100, start2=100, unmapped2=true, attrs=Map("RX" -> "AAAAAAAt"))
      builder.addPair(name="c01", start1=100, start2=300, mapq1=5)

      val in  = builder.toTempFile()
      val out = Files.createTempFile("umi_grouped.", ".sam")
      val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
      val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
      val tool = new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), groupingMetrics=Some(metrics), rawTag="RX", assignTag="MI", strategy=Strategy.Edit, edits=1, minMapQ=Some(30))
      val logs = executeFgbioTool(tool)

      val groups = readBamRecs(out).groupBy(_.name.charAt(0))

      // Group A: 1-4 all passed through into one umi group
      groups('a') should have size 4*2
      groups('a').map(_.name).toSet shouldEqual Set("a01", "a02", "a03", "a04")
      groups('a').map(r => r[String]("MI")).toSet should have size 1

      // 5 separated out into another group due to unmapped mate
      groups('b') should have size 2
      groups('b').map(r => r[String]("MI")).toSet should have size 1
      groups('b').map(r => r[String]("MI")).head should not be groups('a').map(r => r[String]("MI")).head

      // 6 out for low mapq,
      groups.contains('c') shouldBe false

      hist.toFile.exists() shouldBe true

      // TODO: Consider creating more unit tests that vary the following metric fields
      val expectedMetric = UmiGroupingMetric(accepted_sam_records = 10, discarded_non_pf = 0, discarded_poor_alignment = 2, discarded_ns_in_umi = 0, discarded_umis_to_short = 0)
      Metric.read[UmiGroupingMetric](metrics) shouldEqual Seq(expectedMetric)

      // Make sure that we skip sorting for TemplateCoordinate
      val sortMessage = "Sorting the input to TemplateCoordinate order"
      logs.exists(_.contains(sortMessage)) shouldBe (sortOrder != TemplateCoordinate)
    }
  }

  it should "correctly mark duplicates on duplicate reads in group, when flag is passed" in {
    val builder = new SamBuilder(readLength = 100, sort = Some(SamOrder.Coordinate))
    // Mapping Quality is a tiebreaker, so use that to our advantage here.
    builder.addPair(mapq1 = 10, mapq2 = 10, name = "a01", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 30, mapq2 = 30, name = "a02", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 100, mapq2 = 10, name = "a03", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 0, mapq2 = 0, name = "a04", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))

    val in = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    val gr = new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), strategy=Strategy.Paired, edits=1, markDuplicates=true)

    gr.markDuplicates shouldBe true
    gr.execute()

    val recs = readBamRecs(out)
    recs.length shouldBe 8
    recs.filter(_.name.equals("a01")).forall(_.duplicate == true) shouldBe true
    recs.filter(_.name.equals("a02")).forall(_.duplicate == true) shouldBe true
    recs.filter(_.name.equals("a03")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a04")).forall(_.duplicate == true) shouldBe true
  }

  it should "does not mark duplicates on reads in group, when flag is not passed" in {
    val builder = new SamBuilder(readLength = 100, sort = Some(SamOrder.Coordinate))
    // Mapping Quality is a tiebreaker, so use that to our advantage here.
    builder.addPair(mapq1 = 10, mapq2 = 10, name = "a01", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 30, mapq2 = 30, name = "a02", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 100, mapq2 = 10, name = "a03", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))
    builder.addPair(mapq1 = 0, mapq2 = 0, name = "a04", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> "ACT-ACT"))

    val in = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), strategy=Strategy.Paired, edits=1).execute()

    val recs = readBamRecs(out)
    recs.length shouldBe 6
    recs.filter(_.name.equals("a01")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a02")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a03")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a04")).forall(_.duplicate == false) shouldBe true
  }

  it should "correctly mark duplicates on duplicate single-end reads with UMIs" in {
    val builder = new SamBuilder(readLength = 100, sort = Some(SamOrder.Coordinate))
    builder.addFrag(mapq = 100, name = "a01", start = 100, attrs = Map("RX" -> "AAAAAAAA"))
    builder.addFrag(mapq = 10, name = "a02", start = 100, attrs = Map("RX" -> "AAAAAAAA"))
    builder.addFrag(mapq = 100, name = "a03", start = 100, attrs = Map("RX" -> "CACACACA"))
    builder.addFrag(mapq = 10, name = "a04", start = 100, attrs = Map("RX" -> "CACACACC"))

    val in = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input = in, output = out, familySizeHistogram = Some(hist), rawTag = "RX", assignTag = "MI", strategy = Strategy.Edit, edits = 1, markDuplicates = true).execute()

    val recs = readBamRecs(out)
    recs.length shouldBe 4
    recs.filter(_.name.equals("a01")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a02")).forall(_.duplicate == true) shouldBe true
    recs.filter(_.name.equals("a03")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a04")).forall(_.duplicate == true) shouldBe true
  }

  it should "mark duplicates and discard secondary and supplementary reads" in {
    // Mapping Quality is a tiebreaker, so use that to our advantage here.
    val builder = new SamBuilder(readLength = 100, sort = Some(SamOrder.Coordinate))
    Range.inclusive(start=1, end=3).foreach { i =>
      val rec = builder.addFrag(mapq = 100, name = "a01", start = 100, attrs = Map("RX" -> "AAAAAAAA")).value
      rec.secondary = i == 2
      rec.supplementary = i == 3
      rec
    }
    builder.addFrag(mapq = 10, name = "a02", start = 100, attrs = Map("RX" -> "AAAAAAAA")).value

    val in = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input = in, output = out, familySizeHistogram = Some(hist), rawTag = "RX", assignTag = "MI", strategy = Strategy.Edit, edits = 1, markDuplicates = true).execute()

    val recs = readBamRecs(out)
    recs.length shouldBe 2
    recs.filter(_.name.equals("a01")).forall(_.duplicate == false) shouldBe true
    recs.filter(_.name.equals("a02")).forall(_.duplicate == true) shouldBe true
    recs.forall(_.secondary) shouldBe false
    recs.forall(_.supplementary) shouldBe false
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

  it should "correctly group reads with the paired assigner when the two UMIs are the same in cross-contig read pairs" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", contig = 1, contig2 = Some(2), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", contig = 1, contig2 = Some(2), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a03", contig = 1, contig2 = Some(2), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a04", contig = 1, contig2 = Some(2), start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b01", contig = 2, contig2 = Some(1), start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b02", contig = 2, contig2 = Some(1), start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b03", contig = 2, contig2 = Some(1), start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="b04", contig = 2, contig2 = Some(1), start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACT-ACT"))

    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
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

  it should "correctly group together single-end reads with UMIs" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addFrag(name="a01", start=100, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addFrag(name="a02", start=100, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addFrag(name="a03", start=100, attrs=Map("RX" -> "CACACACA"))
    builder.addFrag(name="a04", start=100, attrs=Map("RX" -> "CACACACC"))
    builder.addFrag(name="a05", start=105, attrs=Map("RX" -> "GTAGTAGG"))
    builder.addFrag(name="a06", start=105, attrs=Map("RX" -> "GTAGTAGG"))
    builder.addFrag(name="a07", start=107, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addFrag(name="a08", start=107, attrs=Map("RX" -> "AAAAAAAA"))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    val hist = Files.createTempFile("umi_grouped.", ".histogram.txt")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=Some(hist), rawTag="RX", assignTag="MI", strategy=Strategy.Edit, edits=1).execute()

    val recs = readBamRecs(out)
    recs should have size 8

    val groups = recs.groupBy(r => r[String]("MI")).values.map(rs => rs.map(_.name).toSet)
    groups should have size 4
    groups should contain theSameElementsAs Seq(Set("a01", "a02"), Set("a03", "a04"), Set("a05", "a06"), Set("a07", "a08"))
  }

  it should "correctly group together reads with UMIs across different cells using the cell barcode" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    val specialCellTag = "XX"
    builder.addFrag(name="a01", start=100, attrs=Map("RX" -> "AAAAAAAA", specialCellTag -> "AA"))
    builder.addFrag(name="a02", start=100, attrs=Map("RX" -> "AAAAAAAA", specialCellTag -> "AA"))
    builder.addFrag(name="a03", start=100, attrs=Map("RX" -> "CACACACA", specialCellTag -> "CA"))
    builder.addFrag(name="a04", start=100, attrs=Map("RX" -> "CACACACC", specialCellTag -> "NN"))
    builder.addFrag(name="a05", start=105, attrs=Map("RX" -> "GTAGTAGG", specialCellTag -> "GT"))
    builder.addFrag(name="a06", start=105, attrs=Map("RX" -> "GTAGTAGG", specialCellTag -> "GT"))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", cellTag=specialCellTag, strategy=Strategy.Edit, edits=1).execute()

    val recs = readBamRecs(out)
    recs should have size 6

    val groups = recs.groupBy(r => r[String]("MI")).values.map(rs => rs.map(_.name).toSet)
    groups should have size 4
    groups should contain theSameElementsAs Seq(Set("a01", "a02"), Set("a03"), Set("a04"), Set("a05", "a06"))
  }

  it should "exclude reads that contain an N in the UMI" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a03", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ANN"))

    val in  = builder.toTempFile()
    val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
    val out = Files.createTempFile("umi_grouped.", ".bam")
    new GroupReadsByUmi(input=in, output=out, groupingMetrics=Some(metrics), rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=2).execute()

    val expectedMetric = UmiGroupingMetric(accepted_sam_records = 4, discarded_non_pf = 0, discarded_poor_alignment = 0, discarded_ns_in_umi = 2, discarded_umis_to_short = 0)
    Metric.read[UmiGroupingMetric](metrics) shouldEqual Seq(expectedMetric)

    readBamRecs(out).map(_.name).distinct shouldBe Seq("a01", "a02")
  }

  it should "group paired UMIs of different lengths separately with the paired strategy" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    // Molecule 1 has two 6bp UMIs, molecule 2 two 3bp UMIs; each is seen from both strands (F1R2 and F2R1)
    Seq("a" -> "ACTACT-GGGGGG", "c" -> "ACT-GGG").foreach { case (prefix, umi) =>
      Range.inclusive(1, 4).foreach(i => builder.addPair(name=f"${prefix}0$i", start1=100, start2=300, strand1=Plus, strand2=Minus, attrs=Map("RX" -> umi)))
    }
    Seq("b" -> "GGGGGG-ACTACT", "d" -> "GGG-ACT").foreach { case (prefix, umi) =>
      Range.inclusive(1, 4).foreach(i => builder.addPair(name=f"${prefix}0$i", start1=300, start2=100, strand1=Minus, strand2=Plus, attrs=Map("RX" -> umi)))
    }

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()

    val idsByPrefix = readBamRecs(out).groupBy(_.name.take(1)).map { case (prefix, recs) => prefix -> recs.map(r => r[String]("MI")).distinct }
    idsByPrefix.values.foreach(ids => ids should have size 1)
    val Seq(a, b, c, d) = Seq("a", "b", "c", "d").map(p => idsByPrefix(p).head)

    // a and b are the two strands of molecule 1, and c and d the two strands of molecule 2
    a.takeWhile(_ != '/') shouldBe b.takeWhile(_ != '/')
    c.takeWhile(_ != '/') shouldBe d.takeWhile(_ != '/')
    a should not equal b
    c should not equal d
    a.takeWhile(_ != '/') should not equal c.takeWhile(_ != '/')
  }

  it should "group the two strands of a paired UMI whose two UMIs have different lengths with the paired strategy" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    // The 6bp UMI is always on the read at the lower coordinate, so F1R2 reads are ACTACT-GGG and F2R1 reads GGG-ACTACT
    Range.inclusive(1, 4).foreach(i => builder.addPair(name=f"a0$i", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACT-GGG")))
    Range.inclusive(1, 4).foreach(i => builder.addPair(name=f"b0$i", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "GGG-ACTACT")))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()

    val recs = readBamRecs(out)
    recs.filter(_.name.startsWith("a")).map(r => r[String]("MI")).distinct shouldBe Seq("0/A")
    recs.filter(_.name.startsWith("b")).map(r => r[String]("MI")).distinct shouldBe Seq("0/B")
  }

  it should "not consider read-pairs with the same coordinates but different pair orientations as being from the same molecule" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="f1r2", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addPair(name="f2r1", start1=300, start2=100, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addPair(name="ff",   start1=100, start2=300, strand1=Plus,  strand2=Plus,  attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addPair(name="rr",   start1=  1, start2=201, strand1=Minus, strand2=Minus, attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addPair(name="r1f2", start1=100, start2=300, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addPair(name="r2f1", start1=300, start2=100, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addFrag(name="Frag", start=100, strand=Minus, attrs=Map("RX" -> "ACGT-TTGA"))
    builder.addFrag(name="fRag", start=1,   strand=Plus,  attrs=Map("RX" -> "ACGT-TTGA"))

    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Adjacency, edits=1).execute()

    val recs = readBamRecs(out)
    recs.map(r => r[String]("MI")).distinct.size shouldBe 8  // each template is a separate MI
  }

  Seq(
    Seq("right", "ACT-", "-ACT"),
    Seq("left", "-ACT", "ACT-"),
  ).foreach { case Seq(orientation, leftUmi, rightUmi) =>
    it should s"succeed when run in paired mode but the $orientation end of the source molecule does not have a UMI" in {
      val builder = new SamBuilder(readLength = 100, sort = Some(SamOrder.Coordinate))
      builder.addPair(name = "a01", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> leftUmi))
      builder.addPair(name = "a02", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> leftUmi))
      builder.addPair(name = "a03", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> leftUmi))
      builder.addPair(name = "a04", start1 = 100, start2 = 300, strand1 = Plus, strand2 = Minus, attrs = Map("RX" -> leftUmi))

      builder.addPair(name = "b01", start1 = 300, start2 = 100, strand1 = Minus, strand2 = Plus, attrs = Map("RX" -> rightUmi))
      builder.addPair(name = "b02", start1 = 300, start2 = 100, strand1 = Minus, strand2 = Plus, attrs = Map("RX" -> rightUmi))
      builder.addPair(name = "b03", start1 = 300, start2 = 100, strand1 = Minus, strand2 = Plus, attrs = Map("RX" -> rightUmi))
      builder.addPair(name = "b04", start1 = 300, start2 = 100, strand1 = Minus, strand2 = Plus, attrs = Map("RX" -> rightUmi))

      val in  = builder.toTempFile()
      val out = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input = in, output = out, rawTag = "RX", assignTag = "MI", strategy = Strategy.Paired, edits = 1).execute()

      val recs = readBamRecs(out)
      val aIds = recs.filter(_.name.startsWith("a")).map(r => r[String]("MI")).distinct
      val bIds = recs.filter(_.name.startsWith("b")).map(r => r[String]("MI")).distinct

      aIds should have size 1
      bIds should have size 1
      aIds.head shouldBe "0/A"
      bIds.head shouldBe "0/B"
    }
  }

  it should "fail when the raw tag is not present" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus, strand2=Minus)
    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    val exception = intercept[FailureException] {
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()
    }
    exception.message.value should include ("was missing the raw UMI tag")
  }

  Strategy.values.filterNot(_ == Strategy.Paired).foreach { strategy =>
    it should s"reject records with UMIs that are shorter than the specified minimum length with the $strategy strategy" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))

      builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACT"))
      builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTAC"))
      builder.addPair(name="a03", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACTA"))

      val in   = builder.toTempFile()
      val out  = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=strategy, edits=0, minUmiLength=Some(6)).execute()

      // a02 is too short, and without --truncate the 6bp and 7bp UMIs are separate molecules
      val recs = readBamRecs(out)
      recs.map(_.name) should contain theSameElementsAs Seq("a01", "a01", "a03", "a03")
      recs.groupBy(r => r[String]("MI")).values.map(_.map(_.name).toSet).toSet shouldBe Set(Set("a01"), Set("a03"))
    }

    it should s"truncate to the specified minimum length with --truncate and the $strategy strategy" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
      builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACT"))
      builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTAC"))
      builder.addPair(name="a03", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACTACTA"))

      val in   = builder.toTempFile()
      val out  = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=strategy, edits=0, minUmiLength=Some(5), truncate=true).execute()

      val recs = readBamRecs(out)
      recs.map(_.name) should contain theSameElementsAs Seq("a01", "a01", "a02", "a02", "a03", "a03")
      recs.map(r => r[String]("MI")).distinct should have length 1 // all should be assigned to one molecule
      recs.map(r => r[String]("RX")).distinct should contain theSameElementsAs Seq("ACTACT", "ACTAC", "ACTACTA") // the raw UMIs are not changed
    }

    it should s"truncate to --min-umi-length even when every UMI is longer with --truncate and the $strategy strategy" in {
      // All UMIs have seven bases and differ only in the seventh, so truncating to six bases groups them together
      val umis = Seq("a01" -> "ACTACTC", "a02" -> "ACTACTC", "a03" -> "ACTACTG", "a04" -> "ACTACTG")
      groupFragsByMi(umis, strategy=strategy, edits=0, minUmiLength=Some(6), truncate=true) shouldBe Set(Set("a01", "a02", "a03", "a04"))
    }

    it should s"not group UMIs with the same number of bases but different segment lengths with the $strategy strategy" in {
      val edits = if (strategy == Strategy.Identity) 0 else 1
      // The symmetric dual UMIs differ by length (8+8 vs 7+7), and the asymmetric ones by which segment is longer
      val umis = Seq(
        "a01" -> "AAAAAAAA-CCCCCCCC", "a02" -> "AAAAAAAA-CCCCCCCC", "a03" -> "AAAAAAAA-CCCCCCCG",
        "b01" -> "AAAAAAA-CCCCCCC",   "b02" -> "AAAAAAA-CCCCCCC",   "b03" -> "AAAAAAA-CCCCCCG",
        "c01" -> "ACTACT-GGG",        "c02" -> "ACTACT-GGG",
        "d01" -> "ACTAC-TGGG",        "d02" -> "ACTAC-TGGG",
        "e01" -> "ACGT",              "e02" -> "ACGT-",
      )
      val expected = if (strategy == Strategy.Identity) {
        Set(Set("a01", "a02"), Set("a03"), Set("b01", "b02"), Set("b03"), Set("c01", "c02"), Set("d01", "d02"), Set("e01"), Set("e02"))
      } else {
        Set(Set("a01", "a02", "a03"), Set("b01", "b02", "b03"), Set("c01", "c02"), Set("d01", "d02"), Set("e01"), Set("e02"))
      }
      groupFragsByMi(umis, strategy=strategy, edits=edits, minUmiLength=None) shouldBe expected
    }
  }

  "GroupReadsByUmi" should "group single-end reads with UMIs of different lengths separately" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addFrag(name="a01", start=100, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addFrag(name="a02", start=100, attrs=Map("RX" -> "AAAAAAAA"))
    builder.addFrag(name="a03", start=100, attrs=Map("RX" -> "ACACACA"))
    builder.addFrag(name="a04", start=100, attrs=Map("RX" -> "ACACACC"))
    builder.addFrag(name="a05", start=100, attrs=Map("RX" -> "AAAAAAA"))
    builder.addFrag(name="a06", start=105, attrs=Map("RX" -> "AGTAGG"))
    builder.addFrag(name="a07", start=105, attrs=Map("RX" -> "AGTAGG"))

    val in  = builder.toTempFile()
    val out = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Edit, edits=1).execute()

    // a05 is one edit from a01 and a02, but is shorter, so it is a separate molecule
    val groups = readBamRecs(out).groupBy(r => r[String]("MI")).values.map(_.map(_.name).toSet).toSet
    groups shouldBe Set(Set("a01", "a02"), Set("a03", "a04"), Set("a05"), Set("a06", "a07"))
  }

  it should "discard single-end reads with UMIs shorter than --min-umi-length but not group UMIs of different lengths without --truncate" in {
    val umis = Seq("a01" -> "AGTAGGCC", "a02" -> "AGTAGGCC", "a03" -> "AGTAGGC", "a04" -> "AGTAGGC", "a05" -> "AGTAGG", "a06" -> "AGTAGG")
    groupFragsByMi(umis, strategy=Strategy.Edit, edits=1, minUmiLength=Some(7)) shouldBe Set(Set("a01", "a02"), Set("a03", "a04"))
  }

  it should "group single-end reads with UMIs of different lengths with --min-umi-length and --truncate" in {
    val umis = Seq("a01" -> "AGTAGGCC", "a02" -> "AGTAGGCC", "a03" -> "AGTAGGC", "a04" -> "AGTAGGC", "a05" -> "AGTAGG", "a06" -> "AGTAGG")
    groupFragsByMi(umis, strategy=Strategy.Edit, edits=1, minUmiLength=Some(7), truncate=true) shouldBe Set(Set("a01", "a02", "a03", "a04"))
  }

  it should "reject --truncate without --min-umi-length" in {
    val ex = intercept[ValidationException] { new GroupReadsByUmi(rawTag="RX", assignTag="MI", strategy=Strategy.Edit, edits=1, truncate=true) }
    ex.messages.mkString("\n") should include ("--truncate requires --min-umi-length")
  }

  /** Runs GroupReadsByUmi on single-end reads at the same position and returns the read names grouped by molecule ID.
    * With `Coordinate` input the tool sorts the reads itself, while `TemplateCoordinate` input is used as given. */
  private def groupFragsByMi(umis: Seq[(String, String)],
                             strategy: Strategy,
                             edits: Int,
                             minUmiLength: Option[Int],
                             truncate: Boolean = false,
                             sortOrder: SamOrder = SamOrder.Coordinate): Set[Set[String]] = {
    val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
    umis.foreach { case (name, umi) => builder.addFrag(name=name, start=100, attrs=Map("RX" -> umi)) }

    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=strategy, edits=edits, minUmiLength=minUmiLength, truncate=truncate).execute()

    readBamRecs(out).groupBy(r => r[String]("MI")).values.map(_.map(_.name).toSet).toSet
  }

  it should "reject a --min-umi-length less than one" in {
    Seq(0, -1).foreach { minUmiLength =>
      an[ValidationException] shouldBe thrownBy {
        new GroupReadsByUmi(rawTag="RX", assignTag="MI", strategy=Strategy.Identity, minUmiLength=Some(minUmiLength))
      }
    }
  }

  Seq(false, true).foreach { truncate =>
    it should s"reject --min-umi-length with the paired strategy and --truncate=$truncate" in {
      an[ValidationException] shouldBe thrownBy {
        new GroupReadsByUmi(rawTag="RX", assignTag="MI", strategy=Strategy.Paired, minUmiLength=Some(4), truncate=truncate)
      }
    }
  }

  "GroupReadsByUmi.basesOnly" should "return the UMI unchanged when it contains only upper-case bases" in {
    val umi = "ACGTN"
    GroupReadsByUmi.basesOnly(umi) should be theSameInstanceAs umi
  }

  it should "remove dashes and other characters that are not upper-case bases" in {
    GroupReadsByUmi.basesOnly("ACG-TAC")  shouldBe "ACGTAC"
    GroupReadsByUmi.basesOnly("-ACGT-")   shouldBe "ACGT"
    GroupReadsByUmi.basesOnly("AC+GT.aa") shouldBe "ACGT"
    GroupReadsByUmi.basesOnly("-")        shouldBe ""
    GroupReadsByUmi.basesOnly("")         shouldBe ""
  }

  "GroupReadsByUmi" should "discard reads whose UMIs contain a lower-case n" in {
    groupFragsByMi(Seq("a01" -> "ACGTn", "a02" -> "ACGTA"), strategy=Strategy.Identity, edits=0, minUmiLength=None) shouldBe Set(Set("a02"))
  }

  Seq((None, false), (Some(8), false), (Some(8), true)).foreach { case (minUmiLength, truncate) =>
    it should s"group UMIs that differ only by case with the adjacency strategy, --min-umi-length=$minUmiLength and --truncate=$truncate" in {
      groupFragsByMi(Seq("a01" -> "aaaaaaaa", "a02" -> "AAAAAAAA"), strategy=Strategy.Adjacency, edits=1, minUmiLength=minUmiLength, truncate=truncate) shouldBe Set(Set("a01", "a02"))
    }
  }

  Strategy.values.filterNot(_ == Strategy.Paired).foreach { strategy =>
    it should s"count only bases when rejecting UMIs containing dashes that are too short with the $strategy strategy" in {
      // "ACG-TA" has six characters but only five bases
      groupFragsByMi(Seq("a01" -> "ACG-TA", "a02" -> "ACG-TAC"), strategy=strategy, edits=0, minUmiLength=Some(6)) shouldBe Set(Set("a02"))
    }
  }

  Seq(Strategy.Edit, Strategy.Adjacency).foreach { strategy =>
    it should s"compare only the bases of UMIs containing dashes when truncating with the $strategy strategy" in {
      // Comparing bases gives "ACGTT" vs "ACGTA" (one mismatch); comparing characters would align the dashes against
      // bases, e.g. "AC-GTT" vs "ACGTA-" (four mismatches).
      groupFragsByMi(Seq("a01" -> "AC-GTT", "a02" -> "ACGTA-T"), strategy=strategy, edits=1, minUmiLength=Some(5), truncate=true) shouldBe Set(Set("a01", "a02"))
      groupFragsByMi(Seq("a01" -> "ACGTAC-", "a02" -> "TCGTACG-TT"), strategy=strategy, edits=1, minUmiLength=Some(6), truncate=true) shouldBe Set(Set("a01", "a02"))
      // Both UMIs have six bases, more than --min-umi-length, and "ACGTTA" vs "ACGTAC" has two mismatches.  UMIs are
      // truncated to --min-umi-length, not to the shortest UMI, so "ACGTT" vs "ACGTA" (one mismatch) are grouped.
      groupFragsByMi(Seq("a01" -> "AC-GTTA", "a02" -> "ACGTA-C"), strategy=strategy, edits=1, minUmiLength=Some(5), truncate=true) shouldBe Set(Set("a01", "a02"))
    }

    it should s"not group UMIs whose segments have different lengths without --truncate with the $strategy strategy" in {
      // One mismatch apart if the dashes are ignored, but the segment lengths differ (2+4 vs 5+1)
      groupFragsByMi(Seq("a01" -> "AC-GTTA", "a02" -> "ACGTA-A"), strategy=strategy, edits=1, minUmiLength=Some(5)) shouldBe Set(Set("a01"), Set("a02"))
    }
  }

  // When only identical UMIs are grouped, the tool groups reads by UMI while sorting.  The grouping must be the same
  // whether the tool sorts the input or the input is already TemplateCoordinate sorted.
  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate); (strategy, edits) <- Seq(Strategy.Identity -> 0, Strategy.Edit -> 0, Strategy.Adjacency -> 0)) {
    val desc = s"the $strategy strategy, edits=$edits and $sortOrder input"
    // The first eight bases are ACGTACGT for a01, a02 and a04, and ACGTACGA for a03
    val dashedUmis = Seq("a01" -> "ACG-TACGTA", "a02" -> "ACG-TACGTC", "a03" -> "ACG-TACGA", "a04" -> "ACGT-ACGTA")

    it should s"compare only the first --min-umi-length bases of UMIs containing dashes with --truncate and $desc" in {
      groupFragsByMi(dashedUmis, strategy=strategy, edits=edits, minUmiLength=Some(8), truncate=true, sortOrder=sortOrder) shouldBe Set(Set("a01", "a02", "a04"), Set("a03"))
    }

    it should s"not truncate UMIs containing dashes without --truncate with $desc" in {
      groupFragsByMi(dashedUmis, strategy=strategy, edits=edits, minUmiLength=Some(8), sortOrder=sortOrder) shouldBe Set(Set("a01"), Set("a02"), Set("a03"), Set("a04"))
    }

    it should s"truncate UMIs to --min-umi-length with --truncate and $desc" in {
      val umis = Seq("a01" -> "ACGTA", "a02" -> "ACGTC", "a03" -> "TTTT", "a04" -> "ACGT")
      groupFragsByMi(umis, strategy=strategy, edits=edits, minUmiLength=Some(4), truncate=true, sortOrder=sortOrder) shouldBe Set(Set("a01", "a02", "a04"), Set("a03"))
    }

    it should s"group only identical UMIs of the same length without --min-umi-length with $desc" in {
      val umis = Seq("a01" -> "ACGTA", "a02" -> "ACGTA", "a03" -> "ACGT", "a04" -> "ACGT-", "a05" -> "ACGTAC")
      groupFragsByMi(umis, strategy=strategy, edits=edits, minUmiLength=None, sortOrder=sortOrder) shouldBe Set(Set("a01", "a02"), Set("a03"), Set("a04"), Set("a05"))
    }
  }
}
