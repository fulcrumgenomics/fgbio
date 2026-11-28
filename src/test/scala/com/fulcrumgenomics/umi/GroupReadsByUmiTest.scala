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
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.umi.GroupReadsByUmi._
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.SAMReadGroupRecord
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

  for (cellTag <- Seq("CB", "XX"); fromTemplate <- Seq(false, true); r1First <- Seq(false, true)) {
    it should s"normalize absent and empty $cellTag values with fromTemplate=$fromTemplate and r1First=$r1First" in {
      Seq(None, Some(""), Some("AA")).foreach { barcode =>
        val builder = new SamBuilder(readLength=100, sort=None)
        val Seq(r1, r2) = builder.addPair(start1=if (r1First) 100 else 300, start2=if (r1First) 300 else 100)
        barcode.foreach(value => r1(cellTag) = value)
        val info = if (fromTemplate) ReadInfo(Template(r1=Some(r1), r2=Some(r2)), cellTag)
                   else ReadInfo(r1, cellTag)
        info.cellBarcode shouldBe Some(barcode.getOrElse(""))
      }
    }
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

  it should "fail when umis have different length" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.addPair(name="a01", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-ACT"))
    builder.addPair(name="a02", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "ACT-AC"))

    val in   = builder.toTempFile()
    val out  = Files.createTempFile("umi_grouped.", ".sam")
    val tool = new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1)

    an[Exception] should be thrownBy tool.execute()
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

      val in   = builder.toTempFile()
      val out  = Files.createTempFile("umi_grouped.", ".sam")
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=strategy, edits=0, minUmiLength=Some(6)).execute()

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
      new GroupReadsByUmi(input=in, output=out, familySizeHistogram=None, rawTag="RX", assignTag="MI", strategy=strategy, edits=0, minUmiLength=Some(5)).execute()

      val recs = readBamRecs(out)
      recs should have length 4
      recs.map(_.name) should contain theSameElementsInOrderAs Seq("a01", "a01", "a02", "a02")
      recs.map(r => r[String]("MI")).distinct should have length 1 // all should be assigned to one molecule
    }
  }

  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate); (strategy, edits) <- Seq(Strategy.Identity -> 0, Strategy.Edit -> 1, Strategy.Adjacency -> 1, Strategy.Paired -> 1)) {
    it should s"group reads by library when the names of two libraries are interleaved with the $strategy strategy and $sortOrder input" in {
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      builder.header.getReadGroups.get(0).setLibrary("libA")
      val rgB = new SAMReadGroupRecord("B")
      rgB.setLibrary("libB")
      rgB.setSample("Sample")
      builder.header.addReadGroup(rgB)

      Seq("q1", "q3").foreach(n => builder.addPair(name=n, start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT")))
      Seq("q2", "q4").foreach(n => builder.addPair(name=n, start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT", "RG" -> "B")))

      val out = Files.createTempFile("umi_grouped.", ".bam")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, rawTag="RX", assignTag="MI", strategy=strategy, edits=edits).execute()

      readBamRecs(out).groupBy(r => r[String]("MI")).values.map(_.map(_.name).toSet).toSet shouldBe Set(Set("q1", "q3"), Set("q2", "q4"))
    }
  }

  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate); strategy <- Strategy.values) {
    it should s"group absent and empty cell barcodes together with $strategy and $sortOrder input" in {
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT"))
      builder.addPair(name="q2", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT", "CB" -> ""))
      builder.addPair(name="q3", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT"))
      builder.addPair(name="q4", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT", "CB" -> "AA"))
      val out = Files.createTempFile("umi_grouped.", ".bam")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=strategy, edits=if (strategy == Strategy.Identity) 0 else 1).execute()

      val recs = readBamRecs(out)
      recs.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
        Seq("q1", "q2", "q3", "q4").flatMap(name => Seq(name -> true, name -> false))
      val families = recs.groupBy(r => r[String]("MI")).values.map(_.map(r => (r.name, r.firstOfPair)).toSet).toSeq
      families should contain theSameElementsAs Seq(Seq("q1", "q2", "q3"), Seq("q4"))
        .map(_.flatMap(name => Seq(name -> true, name -> false)).toSet)
      recs.filter(_.name == "q2").map(_.get[String]("CB")) shouldBe Seq(Some(""), Some(""))
      recs.filter(r => r.name == "q1" || r.name == "q3").map(_.get[String]("CB")).toSet shouldBe Set(None)
    }
  }

  it should "explain custom cell tag recovery and group correctly after CB is populated and input re-sorted" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
    Seq("q1" -> "AA", "q2" -> "BB", "q3" -> "AA").foreach { case (name, barcode) =>
      builder.addPair(name=name, start1=100, start2=300, attrs=Map("RX" -> "ACGT", "XX" -> barcode))
    }
    val out = Files.createTempFile("umi_grouped.", ".bam")
    val exception = intercept[FailureException] {
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, cellTag="XX", strategy=Strategy.Adjacency, edits=1).execute()
    }
    val message = exception.message.value
    message should include ("--cell-tag XX")
    message should include ("sort by CB")
    message should include ("populate CB consistently from XX")
    message should include ("re-sorting alone may not fix")

    builder.foreach(rec => rec("CB") = rec[String]("XX"))
    new GroupReadsByUmi(input=builder.toTempFile(), output=out, cellTag="XX", strategy=Strategy.Adjacency, edits=1).execute()
    val recs = readBamRecs(out)
    recs.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
      Seq("q1", "q2", "q3").flatMap(name => Seq(name -> true, name -> false))
    val families = recs.groupBy(r => r[String]("MI")).values.map(_.map(r => (r.name, r.firstOfPair)).toSet).toSeq
    families should contain theSameElementsAs Seq(Seq("q1", "q3"), Seq("q2"))
      .map(_.flatMap(name => Seq(name -> true, name -> false)).toSet)
  }

  /** Adds a second read group, with the given ID and library, to the builder's header. */
  private def addReadGroup(builder: SamBuilder, id: String, library: String): Unit = {
    val rg = new SAMReadGroupRecord(id)
    rg.setLibrary(library)
    rg.setSample("Sample")
    builder.header.addReadGroup(rg)
  }

  it should "group reads by library when the strand flags of unmapped mates vary with the Paired strategy" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.header.getReadGroups.get(0).setLibrary("libA")
    addReadGroup(builder, "B", "libB")

    // The unmapped R2's strand flag is set on q1 and q2 and clear on q3 and q4
    Seq("q1" -> Minus, "q3" -> Plus).foreach { case (n, s) =>
      builder.addPair(name=n, start1=100, unmapped2=true, strand2=s, attrs=Map("RX" -> "ACGT-TTTT"))
    }
    Seq("q2" -> Minus, "q4" -> Plus).foreach { case (n, s) =>
      builder.addPair(name=n, start1=100, unmapped2=true, strand2=s, attrs=Map("RX" -> "ACGT-TTTT", "RG" -> "B"))
    }

    val out = Files.createTempFile("umi_grouped.", ".bam")
    new GroupReadsByUmi(input=builder.toTempFile(), output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Paired, edits=1).execute()

    readBamRecs(out).groupBy(r => r[String]("MI")).values.map(_.map(_.name).toSet).toSet shouldBe Set(Set("q1", "q3"), Set("q2", "q4"))
  }

  it should "group reads from several read groups of the same library together" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
    builder.header.getReadGroups.get(0).setLibrary("libA")
    addReadGroup(builder, "A2", "libA")

    Seq("q1", "q3").foreach(n => builder.addPair(name=n, start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT")))
    Seq("q2", "q4").foreach(n => builder.addPair(name=n, start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT", "RG" -> "A2")))

    val out = Files.createTempFile("umi_grouped.", ".bam")
    new GroupReadsByUmi(input=builder.toTempFile(), output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Adjacency, edits=1).execute()

    readBamRecs(out).map(r => r[String]("MI")).distinct should have length 1
  }

  it should "write reverse-strand templates first, then each library in turn, at one position" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
    builder.header.getReadGroups.get(0).setLibrary("libA")
    addReadGroup(builder, "B", "libB")

    // Lower ends at 5' position 100: forward R1s starting at 100, and reverse R1s ending at 100. Names are chosen
    // so that name order disagrees with both strand and library order.
    builder.addPair(name="a", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "AAAA", "RG" -> "B"))
    builder.addPair(name="b", start1=100, start2=300, strand1=Plus,  strand2=Minus, attrs=Map("RX" -> "CCCC"))
    builder.addPair(name="c", start1=1,   start2=300, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "GGGG", "RG" -> "B"))
    builder.addPair(name="d", start1=1,   start2=300, strand1=Minus, strand2=Plus,  attrs=Map("RX" -> "TTTT"))

    val out = Files.createTempFile("umi_grouped.", ".bam")
    new GroupReadsByUmi(input=builder.toTempFile(), output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Identity, edits=0).execute()

    readBamRecs(out).map(_.name).distinct should contain theSameElementsInOrderAs Seq("d", "c", "b", "a")
  }

  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate);
       strategy <- Strategy.values;
       filter <- Seq("PF", "N-in-UMI", "minimum UMI length") if filter != "minimum UMI length" || strategy != Strategy.Paired;
       failingRead <- Seq(1, 2)) {
    it should s"discard whole templates when R$failingRead fails $filter with $strategy and $sortOrder input" in {
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      builder.header.getReadGroups.get(0).setLibrary("libA")
      val umi = if (strategy == Strategy.Paired) "ACGT-TTTT" else "ACGT"
      Seq("q1", "q2", "q3").foreach { name =>
        val pair = builder.addPair(name=name, start1=100, start2=300, attrs=Map("RX" -> umi))
        if (name == "q2") filter match {
          case "PF"                 => pair(failingRead - 1).pf = false
          case "N-in-UMI"           => pair(failingRead - 1)("RX") = umi.replace('G', 'N')
          case "minimum UMI length" => pair(failingRead - 1)("RX") = "ACG"
        }
      }

      val in   = builder.toTempFile()
      val recs = readBamRecs(in)
      if (sortOrder == SamOrder.TemplateCoordinate) {
        val keys = recs.map(SamOrder.TemplateCoordinate.sortkey)
        keys shouldBe keys.sorted
        recs.map(r => (r.name, r.firstOfPair)) shouldBe
          Seq("q1", "q2", "q3").flatMap(name => Seq(name -> true, name -> false))

        if (failingRead == 2) {
          // Removing only R2 would change the template key without changing R1's record sort key.
          val infos = Seq("q1", "q2", "q3").map { name =>
            val retained = recs.iterator.filter(r => r.name == name && (name != "q2" || r.firstOfPair))
            ReadInfo(Template(retained), cellTag="CB")
          }
          infos.head shouldBe infos.last
          infos(1) should not be infos.head
          infos.map(_.start2) shouldBe Seq(399, Int.MaxValue, 399)
        }
      }

      val out     = Files.createTempFile("umi_grouped.", ".bam")
      val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
      val hist    = Files.createTempFile("umi_grouped.", ".histogram.txt")
      val tool = new GroupReadsByUmi(input=in, output=out, strategy=strategy,
        edits=if (strategy == Strategy.Identity) 0 else 1, markDuplicates=true,
        groupingMetrics=Some(metrics), familySizeHistogram=Some(hist),
        minUmiLength=if (filter == "minimum UMI length") Some(4) else None)
      try tool.execute() catch {
        case exception: FailureException => fail(exception.message.value, exception)
      }

      val output = readBamRecs(out)
      output.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
        Seq("q1", "q3").flatMap(name => Seq(name -> true, name -> false))
      output.map(_[String]("MI")).distinct should have size 1
      output.count(_.duplicate) shouldBe 2
      Metric.read[TagFamilySizeMetric](hist).map(m => (m.family_size, m.count)) shouldBe Seq(2 -> 1L)
      Metric.read[UmiGroupingMetric](metrics) shouldBe Seq(UmiGroupingMetric(4,
        discarded_non_pf=if (filter == "PF") 2 else 0,
        discarded_poor_alignment=0,
        discarded_ns_in_umi=if (filter == "N-in-UMI") 2 else 0,
        discarded_umis_to_short=if (filter == "minimum UMI length") 2 else 0))
    }
  }

  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate);
       strategy <- Strategy.values;
       failingRead <- Seq(1, 2)) {
    it should s"retain non-PF R$failingRead and its mate with $strategy and $sortOrder input when included" in {
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      val umi = if (strategy == Strategy.Paired) "ACGT-TTTT" else "ACGT"
      Seq("q1", "q2", "q3").foreach { name =>
        val pair = builder.addPair(name=name, start1=100, start2=300, attrs=Map("RX" -> umi))
        if (name == "q2") pair(failingRead - 1).pf = false
      }
      val out     = Files.createTempFile("umi_grouped.", ".bam")
      val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=strategy,
        edits=if (strategy == Strategy.Identity) 0 else 1,
        includeNonPfReads=true, groupingMetrics=Some(metrics)).execute()

      val recs = readBamRecs(out)
      recs.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
        Seq("q1", "q2", "q3").flatMap(name => Seq(name -> true, name -> false))
      recs.map(_[String]("MI")).distinct should have size 1
      recs.count(!_.pf) shouldBe 1
      Metric.read[UmiGroupingMetric](metrics) shouldBe Seq(UmiGroupingMetric(6, 0, 0, 0, 0))
    }
  }

  for (sortOrder <- Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate);
       filter <- Seq("PF", "N-in-UMI", "minimum UMI length");
       failingRead <- Seq(1, 2)) {
    it should s"filter whole templates using custom tags when R$failingRead fails $filter with $sortOrder input" in {
      val builder = new SamBuilder(readLength=100, sort=Some(sortOrder))
      Seq("q1", "q2", "q3").foreach { name =>
        val pair = builder.addPair(name=name, start1=100, start2=300,
          attrs=Map("UX" -> "ACGT", "XX" -> "AA", "CB" -> "AA"))
        if (name == "q2") filter match {
          case "PF"                 => pair(failingRead - 1).pf = false
          case "N-in-UMI"           => pair(failingRead - 1)("UX") = "ACGN"
          case "minimum UMI length" => pair(failingRead - 1)("UX") = "ACG"
        }
      }
      val out = Files.createTempFile("umi_grouped.", ".bam")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=Strategy.Identity, edits=0,
        rawTag="UX", assignTag="XI", cellTag="XX",
        minUmiLength=if (filter == "minimum UMI length") Some(4) else None).execute()

      val recs = readBamRecs(out)
      recs.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
        Seq("q1", "q3").flatMap(name => Seq(name -> true, name -> false))
      recs.map(_[String]("XI")).distinct should have size 1
      recs.map(_[String]("UX")).toSet shouldBe Set("ACGT")
      recs.map(_[String]("XX")).toSet shouldBe Set("AA")
    }
  }

  Seq(1, 2).foreach { failingRead =>
    it should s"reject non-PF templates before canonicalizing UMIs when R$failingRead fails PF" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))
      builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "ACGT"))
      val pair = builder.addPair(name="q2", start1=100, start2=300)
      pair(failingRead - 1).pf = false
      val out = Files.createTempFile("umi_grouped.", ".bam")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=Strategy.Identity, edits=0).execute()
      readBamRecs(out).map(r => (r.name, r.firstOfPair)) shouldBe Seq("q1" -> true, "q1" -> false)
    }
  }

  Seq(false, true).foreach { includeNonPfReads =>
    it should s"count rejected templates under the first matching filter with includeNonPfReads=$includeNonPfReads" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
      builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "ACGT"))
      val pair = builder.addPair(name="q2", start1=100, start2=300, attrs=Map("RX" -> "ACGT"))
      pair.head.pf = false
      pair.last("RX") = "ACGN"
      val out     = Files.createTempFile("umi_grouped.", ".bam")
      val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=Strategy.Adjacency, edits=1,
        includeNonPfReads=includeNonPfReads, groupingMetrics=Some(metrics)).execute()

      readBamRecs(out).map(r => (r.name, r.firstOfPair)) shouldBe Seq("q1" -> true, "q1" -> false)
      Metric.read[UmiGroupingMetric](metrics) shouldBe Seq(UmiGroupingMetric(2,
        discarded_non_pf=if (includeNonPfReads) 0 else 2,
        discarded_poor_alignment=0,
        discarded_ns_in_umi=if (includeNonPfReads) 2 else 0,
        discarded_umis_to_short=0))
    }
  }

  it should "ignore failing secondary and supplementary records when filtering primary templates" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "ACGT"))
    Seq(false, true).foreach { supplementary =>
      builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "N")).foreach { rec =>
        rec.secondary = !supplementary
        rec.supplementary = supplementary
        rec.pf = false
      }
    }
    val out     = Files.createTempFile("umi_grouped.", ".bam")
    val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
    new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=Strategy.Adjacency, edits=1,
      groupingMetrics=Some(metrics), minUmiLength=Some(4)).execute()

    val recs = readBamRecs(out)
    recs.map(r => (r.name, r.firstOfPair)) shouldBe Seq("q1" -> true, "q1" -> false)
    recs.forall(r => r.pf && !r.secondary && !r.supplementary) shouldBe true
    Metric.read[UmiGroupingMetric](metrics) shouldBe Seq(UmiGroupingMetric(2, 0, 0, 0, 0))
  }

  Seq(true, false).foreach { lowMapqOnR1 =>
    it should s"discard both mates on valid pre-sorted input when one has low MAPQ with lowMapqOnR1=$lowMapqOnR1" in {
      val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
      Seq("q1", "q2", "q3").foreach { name =>
        builder.addPair(name=name, start1=100, start2=300, attrs=Map("RX" -> "ACGT"),
          mapq1=if (name == "q2" && lowMapqOnR1) 0 else 60,
          mapq2=if (name == "q2" && !lowMapqOnR1) 0 else 60)
      }
      val in = builder.toTempFile()
      val pair = readBamRecs(in).filter(_.name == "q2")
      pair.head[Int]("MQ") shouldBe pair.last.mapq
      pair.last[Int]("MQ") shouldBe pair.head.mapq

      val out     = Files.createTempFile("umi_grouped.", ".bam")
      val metrics = Files.createTempFile("umi_grouped.", ".metrics.txt")
      new GroupReadsByUmi(input=in, output=out, strategy=Strategy.Adjacency, edits=1,
        groupingMetrics=Some(metrics), minMapQ=Some(30)).execute()

      val recs = readBamRecs(out)
      recs.map(r => (r.name, r.firstOfPair)) should contain theSameElementsAs
        Seq("q1", "q3").flatMap(name => Seq(name -> true, name -> false))
      recs.map(_[String]("MI")).distinct should have size 1
      Metric.read[UmiGroupingMetric](metrics) shouldBe Seq(UmiGroupingMetric(4, 0, 2, 0, 0))
    }
  }

  it should "fail on pre-sorted input whose libraries are not contiguous at a position" in {
    val builder = new SamBuilder(readLength=100, sort=None)
    SamOrder.TemplateCoordinate.applyTo(builder.header)
    builder.header.getReadGroups.get(0).setLibrary("libA")
    addReadGroup(builder, "B", "libB")

    // The order fgbio 4.1.1 and earlier wrote: by name, with the two libraries interleaved
    builder.addPair(name="q1", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT"))
    builder.addPair(name="q2", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT", "RG" -> "B"))
    builder.addPair(name="q3", start1=100, start2=300, attrs=Map("RX" -> "ACGT-ACGT"))

    val out       = Files.createTempFile("umi_grouped.", ".bam")
    val exception = the[FailureException] thrownBy {
      new GroupReadsByUmi(input=builder.toTempFile(), output=out, rawTag="RX", assignTag="MI", strategy=Strategy.Adjacency, edits=1).execute()
    }
    exception.message.value should include ("library 'libA' are not contiguous")
  }
}
