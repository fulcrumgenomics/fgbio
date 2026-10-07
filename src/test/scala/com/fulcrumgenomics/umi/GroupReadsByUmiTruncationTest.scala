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

import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}

/** End-to-end regressions for explicit UMI truncation and its output tags. */
class GroupReadsByUmiTruncationTest extends UnitSpec {
  private val strategies = Strategy.values.filterNot(_ == Strategy.Paired)
  private val orders = Seq(SamOrder.Coordinate, SamOrder.TemplateCoordinate)

  /** Groups reads from a generated input and returns the records and captured tool log. */
  private def run(builder: SamBuilder,
                  strategy: Strategy = Strategy.Identity,
                  length: Option[Int] = Some(6),
                  truncate: Boolean = true,
                  assignTag: String = "MI",
                  rawTag: String = "RX"): (IndexedSeq[SamRecord], Seq[String]) = {
    val out = makeTempFile("umi-truncation.", ".bam")
    val tool = new GroupReadsByUmi(input=builder.toTempFile(), output=out, strategy=strategy,
      edits=if (strategy == Strategy.Identity) 0 else 1, minUmiLength=length, truncate=truncate,
      assignTag=assignTag, rawTag=rawTag)
    val logs = executeFgbioTool(tool)
    (readBamRecs(out), logs)
  }

  /** Requires exactly one R1 and one R2 for every expected template before comparing their tags. */
  private def assertPairIdentities(records: IndexedSeq[SamRecord], names: Seq[String]): Unit = {
    val expected = names.flatMap(name => Seq((name, true, false), (name, false, true)))
    records should have size expected.size.toLong
    val _ = records.map(r => (r.name, r.firstOfPair, r.secondOfPair)) should contain theSameElementsAs expected
  }

  "GroupReadsByUmi assignment" should "accept an empty group with every strategy" in {
    Strategy.values.foreach { strategy =>
      new GroupReadsByUmi(strategy=strategy, edits=0).assignUmiGroups(Seq.empty)
    }
  }

  for (strategy <- strategies; order <- orders; assignTag <- Seq("MI", "ZI")) {
    val desc = s"$strategy, $order and --assign-tag $assignTag"
    for ((length, a, b, merged) <- Seq(
      (8, "-ACGTACGT", "ACGTACGT-", false),
      (10, "ACGT-ACGTAC", "ACGTAC-GTAC", false),
      (9, "ACTACT-GGG", "ACTAC-TGGG", false),
      (4, "ACGT-ACGTAC", "ACGTAC-GTAC", true),
      (8, "ACGTACGT-", "ACGTACGT-AAAA", true)
    )) {
      "GroupReadsByUmi truncation" should s"respect segment boundaries for $a/$b at $length bases with $desc" in {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        Seq(a, b).zipWithIndex.foreach { case (umi, i) =>
          builder.addPair(name=s"r$i", start1=100, start2=100, attrs=Map("RX" -> umi))
        }
        val (recs, _) = run(builder, strategy=strategy, length=Some(length), assignTag=assignTag)
        assertPairIdentities(recs, Seq("r0", "r1"))
        recs.map(_[String](assignTag)).distinct should have size (if (merged) 1 else 2)
        if (merged) recs.map(_[String]("RX")).distinct shouldBe Seq(a.take(length))
      }
    }

    it should s"use a fixed prefix regardless of neighboring UMIs with $desc" in {
      for (neighbor <- Seq(None, Some("ACGT"), Some("CCCC"))) {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        (Seq("ACGTAAAA", "ACGTTTTT") ++ neighbor).zipWithIndex.foreach { case (umi, i) =>
          builder.addPair(name=s"r$i", start1=100, start2=100, attrs=Map("RX" -> umi))
        }
        val (recs, _) = run(builder, strategy=strategy, length=Some(4), assignTag=assignTag)
        assertPairIdentities(recs, if (neighbor.isDefined) Seq("r0", "r1", "r2") else Seq("r0", "r1"))
        recs.filter(r => r.name == "r0" || r.name == "r1").map(_[String](assignTag)).distinct should have size 1
      }
    }

    for (umis <- Seq(Seq("ACGTAC", "ACGTACG"), Seq("AC-GTAC", "ACG-TAC"))) {
      it should s"warn once about filter-only behavior for ${umis.mkString("/")} across full position groups with $desc" in {
        val builder = new SamBuilder(readLength=30, sort=Some(order))
        for (start <- Seq(100, 200); umi <- umis) {
          builder.addPair(start1=start, start2=start, attrs=Map("RX" -> umi))
        }
        val (_, logs) = run(builder, strategy=strategy, truncate=false, assignTag=assignTag)
        val warnings = logs.filter(_.contains(" | Warning]"))
        warnings should have size 1
        warnings.head should include ("--min-umi-length")
        warnings.head should include ("filter")
        warnings.head should include ("--truncate")
      }
    }
  }

  for (strategy <- strategies; order <- orders; codec <- Seq(false, true);
       umis <- Seq(Seq("ACGTAC", "ACGTACG"), Seq("AC-GTAC", "ACG-TAC"))) {
    "GroupReadsByUmi output" should s"support consensus with $strategy, $order, codec=$codec and ${umis.mkString("/")}" in {
      val builder = new SamBuilder(readLength=30, baseQuality=30, sort=Some(order))
      umis.zipWithIndex.foreach { case (umi, i) =>
        builder.addPair(name=s"r$i", start1=100, start2=100, bases1="AC" * 15, bases2="AC" * 15,
          attrs=Map("RX" -> umi, "QX" -> umi.map(c => if (c == '-') ' ' else 'I')))
      }
      val grouped = makeTempFile("umi-grouped.", ".bam")
      new GroupReadsByUmi(input=builder.toTempFile(), output=grouped, strategy=strategy,
        edits=if (strategy == Strategy.Identity) 0 else 1, minUmiLength=Some(6), truncate=true).execute()
      val groups = if (umis.head.contains('-')) 2 else 1
      val sources = readBamRecs(grouped)
      assertPairIdentities(sources, Seq("r0", "r1"))
      sources.map(_[String]("MI")).distinct should have size groups.toLong
      sources.foreach(r => r[String]("QX").length shouldBe r[String]("RX").length)
      val out = makeTempFile("umi-consensus.", ".bam")
      val rejects = makeTempFile("umi-rejects.", ".bam")
      if (codec) new CallCodecConsensusReads(input=grouped, output=out, rejects=Some(rejects)).execute()
      else new CallMolecularConsensusReads(input=grouped, output=out, minReads=1, rejects=Some(rejects)).execute()
      val consensus = readBamRecs(out)
      consensus should have size (groups * (if (codec) 1 else 2)).toLong
      consensus.map(_[String]("RX")).toSet shouldBe (if (groups == 1) Set("ACGTAC") else umis.toSet)
      readBamRecs(rejects) shouldBe empty
    }
  }

  for (order <- orders; assignTag <- Seq("MI", "ZI");
       ox <- Seq(None, Some("earliest")); bz <- Seq(None, Some("stray")); qualities <- Seq(false, true)) {
    "GroupReadsByUmi tags" should s"preserve originals by OX alone with $order, $assignTag, OX=$ox, BZ=$bz, QX=$qualities" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      val pair = builder.addPair(name="originals", start1=100, start2=100, attrs=Map("RX" -> "ac-gtacG"))
      pair.last("RX") = "AC-GTATTT"
      val originals = Seq("ac-gtacG", "AC-GTATTT")
      val originalQuals = Seq("12 34567", "AB CDEFGH")
      pair.zip(originalQuals).foreach { case (r, q) =>
        ox.foreach(r("OX") = _)
        bz.foreach(r("BZ") = _)
        if (qualities) r("QX") = q
      }
      val (recs, _) = run(builder, length=Some(5), assignTag=assignTag)
      assertPairIdentities(recs, Seq("originals"))
      recs.sortBy(_.secondOfPair).zipWithIndex.foreach { case (r, i) =>
        r[String]("RX") shouldBe originals(i).take(6)
        r.get[String]("OX") shouldBe ox.orElse(Some(originals(i)))
        r.get[String]("BZ") shouldBe (if (ox.isDefined || !qualities) bz else Some(originalQuals(i)))
        r.get[String]("QX") shouldBe (if (qualities) Some(originalQuals(i).take(6)) else None)
      }

      // Reprocessing must not replace original values or otherwise change the working tags.
      val again = new SamBuilder(readLength=30, sort=Some(order))
      again ++= recs
      val (reprocessed, _) = run(again, length=Some(5), assignTag=assignTag)
      assertPairIdentities(reprocessed, Seq("originals"))
      reprocessed.sortBy(_.secondOfPair).zip(recs.sortBy(_.secondOfPair)).foreach { case (after, before) =>
        Seq("RX", "QX", "OX", "BZ").foreach(tag => after.get[String](tag) shouldBe before.get[String](tag))
      }
    }
  }

  for (order <- orders; umi <- Seq("ACGTAC", "acgtac"); qualities <- Seq(false, true)) {
    it should s"not invent originals or normalize unchanged RX=$umi with QX=$qualities and $order" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      val pair = builder.addPair(name="unchanged", start1=100, start2=100, attrs=Map("RX" -> umi))
      if (qualities) pair.foreach(_("QX") = "123456")
      val (recs, _) = run(builder)
      assertPairIdentities(recs, Seq("unchanged"))
      recs.foreach { r =>
        r[String]("RX") shouldBe umi
        r.get[String]("OX") shouldBe None
        r.get[String]("BZ") shouldBe None
        r.get[String]("QX") shouldBe (if (qualities) Some("123456") else None)
      }
    }
  }

  for (order <- orders; assignTag <- Seq("MI", "ZI")) {
    it should s"save originals only on the mate whose RX changes with $order and $assignTag" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      val pair = builder.addPair(name="one-changed", start1=100, start2=100, attrs=Map("RX" -> "ACGTAC", "QX" -> "123456"))
      pair.last("RX") = "ACGTACG"
      pair.last("QX") = "1234567"
      val (recs, _) = run(builder, assignTag=assignTag)
      assertPairIdentities(recs, Seq("one-changed"))
      val r1 = recs.find(_.firstOfPair).get
      val r2 = recs.find(_.secondOfPair).get
      r1.get[String]("OX") shouldBe None
      r1.get[String]("BZ") shouldBe None
      r2[String]("OX") shouldBe "ACGTACG"
      r2[String]("BZ") shouldBe "1234567"
      recs.foreach { r =>
        r[String]("RX") shouldBe "ACGTAC"
        r[String]("QX") shouldBe "123456"
      }
    }
  }

  for (order <- orders; umi <- Seq("ACGTAC", "ACGTACG"); mate <- Seq(0, 1)) {
    it should s"reject mismatched RX/QX on mate $mate even if RX=$umi is unchanged with $order" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      val pair = builder.addPair(name="bad-qualities", start1=100, start2=100, attrs=Map("RX" -> umi))
      pair(mate)("QX") = "12345"
      val ex = intercept[FailureException] { run(builder) }
      val message = ex.message.getOrElse("")
      Seq("bad-qualities", "RX", "QX", umi.length.toString, "5").foreach(message should include (_))
    }
  }

  it should "require the standard raw tag for truncation, but allow custom tags for filtering" in {
    val ex = intercept[ValidationException] {
      new GroupReadsByUmi(strategy=Strategy.Identity, edits=0, rawTag="ZA", minUmiLength=Some(6), truncate=true)
    }
    ex.messages.mkString(" ") should include ("--raw-tag RX")
    ex.messages.mkString(" ") should include ("consensus")
    val builder = new SamBuilder(readLength=30, sort=Some(SamOrder.Coordinate))
    Seq("ACGTA", "ACGTAC", "ACGTACG").foreach { umi =>
      builder.addPair(name=umi, start1=100, start2=100, attrs=Map("ZA" -> umi))
    }
    val (recs, _) = run(builder, truncate=false, rawTag="ZA")
    assertPairIdentities(recs, Seq("ACGTAC", "ACGTACG"))
    recs.map(_[String]("MI")).distinct should have size 2
    recs.foreach(_.get[String]("RX") shouldBe None)
  }

  for (order <- orders; scenario <- Seq("no minimum", "truncate", "same lengths", "different positions", "filtered short", "filtered N", "different cells", "different libraries", "different orientations")) {
    "GroupReadsByUmi warnings" should s"not warn for $scenario with $order" in {
      val builder = new SamBuilder(readLength=30, sort=Some(order))
      if (scenario == "different libraries") {
        val rg = builder.header.getReadGroups.get(0)
        rg.setLibrary("one")
        val other = new htsjdk.samtools.SAMReadGroupRecord("other", rg)
        other.setLibrary("two")
        builder.header.addReadGroup(other)
      }
      val secondUmi = scenario match {
        case "same lengths"   => "TCGTAC"
        case "filtered short" => "ACGTA"
        case "filtered N"     => "ACGTACN"
        case _                => "ACGTACG"
      }
      builder.addPair(name="a", start1=100, start2=100, attrs=Map("RX" -> "ACGTAC", "CB" -> "one"))
      val start = if (scenario == "different positions") 200 else 100
      val pair = builder.addPair(name="b", start1=start, start2=start, attrs=Map("RX" -> secondUmi,
        "CB" -> (if (scenario == "different cells") "two" else "one")))
      if (scenario == "different libraries") pair.foreach(_("RG") = "other")
      if (scenario == "different orientations") {
        pair.head.negativeStrand = true
        pair.last.negativeStrand = false
        htsjdk.samtools.SamPairUtil.setMateInfo(pair.head.asSam, pair.last.asSam, true)
      }
      val (_, logs) = run(builder, length=if (scenario == "no minimum") None else Some(6), truncate=scenario == "truncate")
      logs.filter(_.contains(" | Warning]")) shouldBe empty
    }
  }
}
