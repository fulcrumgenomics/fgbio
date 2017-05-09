/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.umi.CorrectUmis.{UmiCorrectionMetrics, UmiMatch}
import com.fulcrumgenomics.util.{Io, Metric}
import com.fulcrumgenomics.sopt.cmdline.ValidationException

class CorrectUmisTest extends UnitSpec {
  private val FixedUmis = Seq("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT")
  private val FixedUmisArray = FixedUmis.toArray
  private val NoBam = makeTempFile("correct_umis.", ".bam")

  "CorrectUmis.findBestMatch" should "find perfect matches" in {
    val corrector = new CorrectUmis(input = NoBam, output = NoBam, maxMismatches = 2, minDistance = 2, umis = FixedUmis)
    val hit1 = corrector.findBestMatch("AAAAAA", FixedUmisArray)
    hit1.matched shouldBe true
    hit1.mismatches shouldBe 0
    hit1.umi shouldBe "AAAAAA"

    val hit2 = corrector.findBestMatch("CCCCCC", FixedUmisArray)
    hit2.matched shouldBe true
    hit2.mismatches shouldBe 0
    hit2.umi shouldBe "CCCCCC"

    val hit3 = corrector.findBestMatch("GGGGGG", FixedUmisArray)
    hit3.matched shouldBe true
    hit3.mismatches shouldBe 0
    hit3.umi shouldBe "GGGGGG"

    val hit4 = corrector.findBestMatch("TTTTTT", FixedUmisArray)
    hit4.matched shouldBe true
    hit4.mismatches shouldBe 0
    hit4.umi shouldBe "TTTTTT"
  }

  it should "match UMIs with up to the maximum allowed mismatches" in {
    val corrector = new CorrectUmis(input = NoBam, output = NoBam, maxMismatches = 2, minDistance = 2, umis = FixedUmis)
    corrector.findBestMatch("AAAAAA", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 0)
    corrector.findBestMatch("AAAAAT", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 1)
    corrector.findBestMatch("AAAACT", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 2)
    corrector.findBestMatch("AAAGCT", FixedUmisArray) shouldBe UmiMatch(matched = false, "AAAAAA", 3)
  }

  it should "not match with mismatches when two UMIs are too similar (e.g. poorly picked UMIs)" in {
    val corrector = new CorrectUmis(input = NoBam, output = NoBam, maxMismatches = 2, minDistance = 2, umis = Seq("AAAG", "AAAT"))
    corrector.findBestMatch("AAAG", Array("AAAG", "AAAT")) shouldBe UmiMatch(matched = false, "AAAG", 0)
  }

  it should "match with lots of mismatches, but only if the next best UMI is far enough away" in {
    val corrector = new CorrectUmis(input = NoBam, output = NoBam, maxMismatches = 3, minDistance = 2, umis = FixedUmis)
    corrector.findBestMatch("AAAAAA", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 0)
    corrector.findBestMatch("AAACGT", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 3)
    corrector.findBestMatch("AAACGT", FixedUmisArray) shouldBe UmiMatch(matched = true, "AAAAAA", 3)
    corrector.findBestMatch("AAACCC", FixedUmisArray) shouldBe UmiMatch(matched = false, "AAAAAA", 3)
    corrector.findBestMatch("AAACCT", FixedUmisArray) shouldBe UmiMatch(matched = false, "AAAAAA", 3) // still too close - three to AAAAAA, and four to CCCCCC
  }

  "CorrectUmis.findUmiPairsWithinDistance" should "find no pairs within the edit distance" in {
    CorrectUmis.findUmiPairsWithinDistance(Seq("AAAA", "TTTT", "CCCC", "GGGG"), 2) shouldBe Seq.empty
  }

  it should "find the right set of pairs" in {
    val umis     = Seq("ACACAC", "CTCTCT", "GAGAGA", "TGTGTG", "ACAGAC", "AGAGAG")
    val expected = Seq(("ACACAC", "ACAGAC", 1),("ACAGAC", "AGAGAG", 2))
    CorrectUmis.findUmiPairsWithinDistance(umis, 2) shouldBe expected
  }

  "CorrectUmis" should "rejects reads that do not have umis, and emit an error" in {
    val builder = new SamRecordSetBuilder(readLength = 10)
    builder.addFrag(name = "q1", start = 1)
    val input = builder.toTempFile()
    val corrected = makeTempFile("corrected.", ".bam")
    val rejects = makeTempFile("rejects.", ".bam")

    val corrector = new CorrectUmis(input = input, output = corrected, rejects = Some(rejects), maxMismatches = 2, minDistance = 2, umis = FixedUmis)
    val logLines = executeFgbioTool(corrector)
    logLines.exists(line => line.contains("Error")) shouldBe true
    readBamRecs(corrected) shouldBe empty
    readBamRecs(rejects) should have size 1
  }

  it should "throw an exception if all the fixed umis are not the same length " in {
    an[ValidationException] shouldBe thrownBy {
      new CorrectUmis(input = NoBam, output = NoBam, maxMismatches = 2, minDistance = 2, umis = Seq("AAAAAA", "CCC")).execute()
    }
  }

  it should "reject reads with incorrect length UMIs and emit an error" in {
    val builder = new SamRecordSetBuilder(readLength = 10)
    builder.addFrag(name = "q1", start = 1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "ACGT"))
    val input = builder.toTempFile()
    val corrected = makeTempFile("corrected.", ".bam")
    val rejects = makeTempFile("rejects.", ".bam")

    val corrector = new CorrectUmis(input = input, output = corrected, rejects = Some(rejects), maxMismatches = 2, minDistance = 2, umis = FixedUmis)
    val logLines = executeFgbioTool(corrector)
    logLines.exists(line => line.contains("Error")) shouldBe true
    readBamRecs(corrected) shouldBe empty
    readBamRecs(rejects) should have size 1
  }

  it should "run end to end and do the right thing for a handful of fragments with single UMIs" in {
    val builder = new SamRecordSetBuilder(readLength = 10)
    builder.addFrag(name="q1",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA"))
    builder.addFrag(name="q2",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA"))
    builder.addFrag(name="q3",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AATAAA"))
    builder.addFrag(name="q4",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AATTAA"))
    builder.addFrag(name="q5",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAATGC"))
    builder.addFrag(name="q6",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "CCCCCC"))
    builder.addFrag(name="q7",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "GGGGGG"))
    builder.addFrag(name="q8",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "TTTTTT"))
    builder.addFrag(name="q9",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "GGGTTT"))
    builder.addFrag(name="q10", start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAACCC"))

    val input     = builder.toTempFile()
    val corrected = makeTempFile("corrected.", ".bam")
    val rejects   = makeTempFile("rejects.", ".bam")
    val metrics   = makeTempFile("metrics.", ".txt")

    val corrector = new CorrectUmis(input=input, output=corrected, rejects=Some(rejects), metrics=Some(metrics), maxMismatches=3, minDistance=2, umis=FixedUmis)
    val logLines = executeFgbioTool(corrector)
    logLines.exists(line => line.contains("Error")) shouldBe false // No errors this time

    readBamRecs(corrected).map(_.getReadName) should contain theSameElementsAs Seq("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8")
    readBamRecs(rejects).map(_.getReadName) should contain theSameElementsAs Seq("q9", "q10")

    val metricsByUmi = Metric.read[UmiCorrectionMetrics](metrics).map(m => m.umi -> m).toMap
    metricsByUmi("AAAAAA") shouldBe UmiCorrectionMetrics(umi="AAAAAA", total_matches=5, perfect_matches=2, one_mismatch_matches=1, two_mismatch_matches=1, other_matches=1, fraction_of_matches = 5/10.0, representation = (5/8.0) / (1/4d))
    metricsByUmi("CCCCCC") shouldBe UmiCorrectionMetrics(umi="CCCCCC", total_matches=1, perfect_matches=1, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 1/10.0, representation = (1/8.0) / (1/4d))
    metricsByUmi("GGGGGG") shouldBe UmiCorrectionMetrics(umi="GGGGGG", total_matches=1, perfect_matches=1, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 1/10.0, representation = (1/8.0) / (1/4d))
    metricsByUmi("TTTTTT") shouldBe UmiCorrectionMetrics(umi="TTTTTT", total_matches=1, perfect_matches=1, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 1/10.0, representation = (1/8.0) / (1/4d))
    metricsByUmi("NNNNNN") shouldBe UmiCorrectionMetrics(umi="NNNNNN", total_matches=2, perfect_matches=0, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 2/10.0, representation = (2/8.0) / (1/4d))
  }

  it should "run end to end and do the right thing for a handful of fragments with duplex UMIs" in {
    val builder = new SamRecordSetBuilder(readLength = 10)
    builder.addFrag(name="q1",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA-CCCCCC")) // ok-ok = ok
    builder.addFrag(name="q2",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAACAA-CCCCAC")) // ok-ok = ok
    builder.addFrag(name="q3",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA-ACTACT")) // ok-ko = reject
    builder.addFrag(name="q4",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "GGGGGG-TTTTTT")) // ok-ok = ok
    builder.addFrag(name="q5",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "GCGCGC-TTTTTT")) // ko-ok = reject

    val input     = builder.toTempFile()
    val corrected = makeTempFile("corrected.", ".bam")
    val rejects   = makeTempFile("rejects.", ".bam")
    val metrics   = makeTempFile("metrics.", ".txt")

    val corrector = new CorrectUmis(input=input, output=corrected, rejects=Some(rejects), metrics=Some(metrics), maxMismatches=3, minDistance=2, umis=FixedUmis)
    val logLines = executeFgbioTool(corrector)
    logLines.exists(line => line.contains("Error")) shouldBe false // No errors this time

    readBamRecs(corrected).map(_.getReadName) should contain theSameElementsAs Seq("q1", "q2", "q4")
    readBamRecs(rejects).map(_.getReadName) should contain theSameElementsAs Seq("q3", "q5")

    val metricsByUmi = Metric.read[UmiCorrectionMetrics](metrics).map(m => m.umi -> m).toMap
    metricsByUmi("AAAAAA") shouldBe UmiCorrectionMetrics(umi="AAAAAA", total_matches=3, perfect_matches=2, one_mismatch_matches=1, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 3/10.0, representation = (3/8.0) / (1/4d))
    metricsByUmi("CCCCCC") shouldBe UmiCorrectionMetrics(umi="CCCCCC", total_matches=2, perfect_matches=1, one_mismatch_matches=1, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 2/10.0, representation = (2/8.0) / (1/4d))
    metricsByUmi("GGGGGG") shouldBe UmiCorrectionMetrics(umi="GGGGGG", total_matches=1, perfect_matches=1, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 1/10.0, representation = (1/8.0) / (1/4d))
    metricsByUmi("TTTTTT") shouldBe UmiCorrectionMetrics(umi="TTTTTT", total_matches=2, perfect_matches=2, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 2/10.0, representation = (2/8.0) / (1/4d))
    metricsByUmi("NNNNNN") shouldBe UmiCorrectionMetrics(umi="NNNNNN", total_matches=2, perfect_matches=0, one_mismatch_matches=0, two_mismatch_matches=0, other_matches=0, fraction_of_matches = 2/10.0, representation = (2/8.0) / (1/4d))
  }

  it should "do the same thing whether the UMIs are on the command line or in a file" in {
    val builder = new SamRecordSetBuilder(readLength = 10)
    builder.addFrag(name="q1",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA"))
    builder.addFrag(name="q2",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAAAAA"))
    builder.addFrag(name="q3",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AATAAA"))
    builder.addFrag(name="q4",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AATTAA"))
    builder.addFrag(name="q5",  start=1).foreach(rec => rec.setAttribute(ConsensusTags.UmiBases, "AAATGC"))

    val input     = builder.toTempFile()
    val corrected = makeTempFile("corrected.", ".bam")
    val rejects   = makeTempFile("rejects.",   ".bam")
    val umiFile   = makeTempFile("umis.",      ".txt")
    val metrics1  = makeTempFile("metrics.",   ".txt")
    val metrics2  = makeTempFile("metrics.",   ".txt")

    Io.writeLines(umiFile, FixedUmis)

    new CorrectUmis(input=input, output=corrected, rejects=Some(rejects), metrics=Some(metrics1), maxMismatches=3, minDistance=2, umis=FixedUmis).execute()
    new CorrectUmis(input=input, output=corrected, rejects=Some(rejects), metrics=Some(metrics2), maxMismatches=3, minDistance=2, umiFiles=Seq(umiFile)).execute()

    val m1 = Metric.read[UmiCorrectionMetrics](metrics1)
    val m2 = Metric.read[UmiCorrectionMetrics](metrics2)
    m1 shouldBe m2
  }
}
