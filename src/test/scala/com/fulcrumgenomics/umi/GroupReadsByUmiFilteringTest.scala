/**
 * Copyright (c) 2026, Fulcrum Genomics LLC
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
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric

/** UMI filters must reject whole templates rather than leaving orphaned mates for grouping. */
class GroupReadsByUmiFilteringTest extends UnitSpec {
  private val strategies = Seq(Strategy.Identity -> 0, Strategy.Edit -> 0, Strategy.Edit -> 1,
    Strategy.Adjacency -> 0, Strategy.Adjacency -> 1, Strategy.Paired -> 0, Strategy.Paired -> 1)
  private val orders = Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate)

  /** Executes the command and reads its output, counters, and final statistics. */
  private def run(builder: SamBuilder,
                  strategy: Strategy,
                  edits: Int,
                  assignTag: String,
                  minimum: Option[Int] = None,
                  truncate: Boolean = false,
                  rawTag: String = "RX"): (IndexedSeq[SamRecord], UmiGroupingMetric, Seq[String]) = {
    val out = makeTempFile("umi-filtering.", ".bam")
    val metrics = makeTempFile("umi-filtering.", ".metrics.txt")
    val tool = new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=strategy,
      edits=edits, assignTag=assignTag, rawTag=rawTag, minUmiLength=minimum, truncate=truncate,
      groupingMetrics=Some(metrics))
    val logs = executeFgbioTool(tool)
    val counters = Metric.read[UmiGroupingMetric](metrics)
    counters should have size 1
    (readBamRecs(out), counters.head, logs)
  }

  /** Asserts both mates survive exactly once for every expected template. */
  private def assertPairs(records: IndexedSeq[SamRecord], names: Seq[String]): Unit = {
    records should have size (names.size * 2).toLong
    records.forall(_.paired) shouldBe true
    val expected = names.flatMap(name => Seq((name, true, false), (name, false, true)))
    val _ = records.map(r => (r.name, r.firstOfPair, r.secondOfPair)) should contain theSameElementsAs expected
  }

  "GroupReadsByUmi sorting" should "keep retained mates together when their raw UMIs differ" in {
    val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.Coordinate))
    val pair = builder.addPair(name="a", start1=100, start2=300, attrs=Map("RX" -> "ACGTAC"))
    pair.last("RX") = "ACGTACG"
    builder.addPair(name="z", start1=100, start2=300, attrs=Map("RX" -> "ACGTAC"))

    val (recs, counters, _) = run(builder, Strategy.Identity, 0, "MI", minimum=Some(6))
    assertPairs(recs, Seq("a", "z"))
    recs.map(_[String]("MI")).distinct should have size 1
    recs.map(r => (r.name, r.firstOfPair) -> r[String]("RX")).toMap shouldBe Map(
      ("a", true) -> "ACGTAC", ("a", false) -> "ACGTACG",
      ("z", true) -> "ACGTAC", ("z", false) -> "ACGTAC"
    )
    counters shouldBe UmiGroupingMetric(4, 0, 0, 0, 0)
  }

  for ((strategy, edits) <- strategies; order <- orders; assignTag <- Seq("MI", "ZI")) {
    val desc = s"$strategy, edits=$edits, $order and $assignTag"
    val valid = if (strategy == Strategy.Paired) "ACG-TAC" else "ACGTAC"

    // Later-sorting neighbor names expose per-record UMI sorting keys that separate rejected templates' mates.
    for (base <- Seq('n', 'N'); mate <- Seq(0, 1)) {
      "GroupReadsByUmi filtering" should s"discard both mates when mate $mate contains $base with $desc" in {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        val pair = builder.addPair(name="rejected", start1=100, start2=300, attrs=Map("RX" -> valid))
        pair(mate)("RX") = valid.dropRight(1) + base
        builder.addPair(name="z-neighbor", start1=100, start2=300, attrs=Map("RX" -> valid))
        val (recs, counters, logs) = run(builder, strategy, edits, assignTag)
        assertPairs(recs, Seq("z-neighbor"))
        counters shouldBe UmiGroupingMetric(2, 0, 0, 2, 0)
        logs.exists(_.contains("Accepted 2 SAM records for grouping.")) shouldBe true
        logs.exists(_.contains("Filtered out 2 SAM records")) shouldBe true
        logs.filter(_.contains(" | Warning]")) shouldBe empty
      }
    }

    for (mate <- Seq(0, 1); truncate <- (if (strategy == Strategy.Paired) Seq(false) else Seq(false, true))) {
      it should s"discard both mates when mate $mate has too few UMI bases with truncate=$truncate and $desc" in {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        val pair = builder.addPair(name="rejected", start1=100, start2=300, attrs=Map("RX" -> "AC-GTAC"))
        // Six characters but only five bases: the length filter must count bases, not characters.
        pair(mate)("RX") = "AC-GTA"
        builder.addPair(name="z-neighbor", start1=100, start2=300, attrs=Map("RX" -> "AC-GTAC"))
        val (recs, counters, logs) = run(builder, strategy, edits, assignTag, minimum=Some(6), truncate=truncate)
        assertPairs(recs, Seq("z-neighbor"))
        counters shouldBe UmiGroupingMetric(2, 0, 0, 0, 2)
        logs.exists(_.contains("Accepted 2 SAM records for grouping.")) shouldBe true
        logs.exists(_.contains("Filtered out 2 SAM records")) shouldBe true
        logs.filter(_.contains(" | Warning]")) shouldBe empty
      }
    }

    for (mate <- Seq(0, 1)) {
      it should s"count N-containing templates before short templates when mate $mate has n with $desc" in {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        val bothReasons = builder.addPair(name="both-reasons", start1=100, start2=300, attrs=Map("RX" -> "AC-GTA"))
        bothReasons(mate)("RX") = "AC-GTAn"
        builder.addPair(name="short", start1=100, start2=300, attrs=Map("RX" -> "AC-GTA"))
        builder.addPair(name="z-neighbor", start1=100, start2=300, attrs=Map("RX" -> "AC-GTAC"))
        val (recs, counters, _) = run(builder, strategy, edits, assignTag, minimum=Some(6))
        assertPairs(recs, Seq("z-neighbor"))
        counters shouldBe UmiGroupingMetric(2, 0, 0, 2, 2)
      }
    }
  }

  for ((strategy, edits) <- strategies.filterNot(_._1 == Strategy.Paired); order <- orders; assignTag <- Seq("MI", "ZI")) {
    it should s"preserve legitimate fragments while applying template UMI filters with $strategy, edits=$edits, $order and $assignTag" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      Seq("fragment" -> "ACGTAC", "upper-n" -> "ACGTAN", "lower-n" -> "ACGTAn", "short" -> "AC-GTA").foreach {
        case (name, umi) => builder.addFrag(name=name, start=100, attrs=Map("RX" -> umi))
      }
      builder.addPair(name="z-neighbor", start1=100, start2=300, attrs=Map("RX" -> "ACGTAC"))
      val (recs, counters, _) = run(builder, strategy, edits, assignTag, minimum=Some(6))
      recs should have size 3
      recs.filter(_.unpaired).map(_.name) shouldBe Seq("fragment")
      assertPairs(recs.filter(_.paired), Seq("z-neighbor"))
      counters shouldBe UmiGroupingMetric(3, 0, 0, 2, 1)
    }
  }

  for (order <- orders; assignTag <- Seq("MI", "ZI"); reject <- Seq("n", "short")) {
    it should s"apply whole-template $reject filtering to custom raw tags with $order and $assignTag" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      val pair = builder.addPair(name="rejected", start1=100, start2=300, attrs=Map("ZA" -> "ACGTAC"))
      pair.head("ZA") = if (reject == "n") "ACGTAn" else "AC-GTA"
      builder.addPair(name="z-neighbor", start1=100, start2=300, attrs=Map("ZA" -> "ACGTAC"))
      val (recs, counters, _) = run(builder, Strategy.Identity, 0, assignTag, minimum=Some(6), rawTag="ZA")
      assertPairs(recs, Seq("z-neighbor"))
      counters shouldBe (if (reject == "n") UmiGroupingMetric(2, 0, 0, 2, 0) else UmiGroupingMetric(2, 0, 0, 0, 2))
      recs.foreach { r =>
        r[String]("ZA") shouldBe "ACGTAC"
        r.get[String]("RX") shouldBe None
      }
    }
  }
}
