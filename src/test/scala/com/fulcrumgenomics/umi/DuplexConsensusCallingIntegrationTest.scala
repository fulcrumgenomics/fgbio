/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus, Strand}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam}
import com.fulcrumgenomics.commons.util.CaptureSystemStreams

/** Integration tests for calling duplex consensus reads.
  *
  * Runs:
  * - GroupReadsByUmi
  * - CallDuplexConsensusReads
  * - FilterConsensusReads
  * */
class DuplexConsensusCallingIntegrationTest extends UnitSpec with CaptureSystemStreams {

  /** Make a reference file that is 100 lines of 100 As. */
  private val ref = {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("A", 5000) // 5000 bases
    builder.toTempFile()
  }

  /** Creates a read pair where the bases are all As and sets the RX tag to "AAAAAAAA"*/
  private def addPair(builder: SamBuilder, name: String, start1: Int, start2: Int, strand1: Strand=Plus, strand2: Strand=Minus, base: Char = 'A', umiBases: String = ""): Seq[SamRecord] = {
    val rxValue = if (umiBases.nonEmpty) umiBases else {
      if (strand1 == Plus) "ACG-GCT" else "GCT-ACG"
    }
    builder.addPair(name=name, start1=start1, start2=start2, strand1=strand1, strand2=strand2, attrs=Map("RX" -> rxValue)).map { rec =>
      rec.bases = base.toString * builder.readLength
      rec
    }
  }

  /** Create some various test cases. */
  private val (inputBam, inputRecords) = {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.Coordinate))

    // duplex 0: one FR read pair -> duplex from one strand (if pass filters)
    addPair(builder, "duplex0:r1", start1=2, start2=101, strand1=Plus, strand2=Minus)

    // duplex 1: one FR read pair and one RF read pair -> duplex from both strand (if pass filters)
    addPair(builder, "duplex1:r1:FR", start1=3, start2=102, strand1=Plus, strand2=Minus)
    addPair(builder, "duplex1:r1:RF", start1=102, start2=3, strand1=Minus, strand2=Plus)

    // duplex 2: five FR read pairs and one RF read pair -> duplex from both strand (if pass filters)
    Range.inclusive(1, 5).foreach { i => addPair(builder, s"duplex3:r$i:FR", start1=4, start2=103, strand1=Plus, strand2=Minus) }
    addPair(builder, "duplex2:r1:RF", start1=103, start2=4, strand1=Minus, strand2=Plus)

    // duplex 3: five FR read pairs and four RF read pair -> duplex from both strand (if pass filters)
    Range.inclusive(1, 5).foreach { i => addPair(builder, s"duplex3:r$i:FR", start1=5, start2=104, strand1=Plus, strand2=Minus) }
    Range.inclusive(1, 4).foreach { i => addPair(builder, s"duplex3:r$i:RF", start1=104, start2=5, strand1=Minus, strand2=Plus) }


    // duplex 4: same as duplex 3, but different RX (C-C instead of A-C)
    Range.inclusive(1, 5).foreach { i => addPair(builder, s"duplex4:r$i:FR", start1=5, start2=104, strand1=Plus, strand2=Minus, umiBases="CTG-GTA") }
    Range.inclusive(1, 4).foreach { i => addPair(builder, s"duplex4:r$i:RF", start1=104, start2=5, strand1=Minus, strand2=Plus, umiBases="GTA-CTG") }

    // duplex 5: five FR read pairs and four RF read pair but wrong bases -> wont call consensus since bases are different
    Range.inclusive(1, 5).foreach { i =>
      addPair(builder, s"duplex5:r$i:FR", start1=6, start2=105, strand1=Plus, strand2=Minus, base='A')
    }
    Range.inclusive(1, 4).foreach { i =>
      addPair(builder, s"duplex5:r$i:RF", start1=105, start2=6, strand1=Minus, strand2=Plus, base='C')
    }

    (builder.toTempFile(), builder.toIndexedSeq)
  }

  private case class ConsensusPaths
  (
    groupedBam: PathToBam          = makeTempFile("grouped.", ".bam"),
    consensusBam: PathToBam        = makeTempFile("consensus.", ".bam"),
    filteredBam: PathToBam         = makeTempFile("filtered.", ".bam"),
    familySizeHistogram: FilePath  = makeTempFile("grouped.", ".family_size_histogram.txt")
  )

  private def call(minReadsForCalling: Seq[Int], minReadsForFiltering: Seq[Int]): ConsensusPaths = {
    val paths = ConsensusPaths()

    // Execute the tools
    captureItAll { () =>
      new GroupReadsByUmi(input=inputBam, output=paths.groupedBam, familySizeHistogram=Some(paths.familySizeHistogram), strategy=Strategy.Paired).execute()
      new CallDuplexConsensusReads(input=paths.groupedBam, output=paths.consensusBam, minReads=minReadsForCalling).execute()
      new FilterConsensusReads(input=paths.consensusBam, output=paths.filteredBam, minReads=minReadsForFiltering, ref=ref, minBaseQuality=10.toByte).execute()
    }

    paths
  }

  /** Extracts the duplex's duplex name from the record's name.  By convention, this is the first colon-separated segment. */
  private def toDuplexName(rec: SamRecord): String = rec.name.split(':').head

  "Duplex Consensus Calling" should "support producing a final BAM with duplex consensuses from one strand" in {
    // duplex consensuses from groups 0-3, and 6
    val paths = call(minReadsForCalling=Seq(1,1,0), minReadsForFiltering=Seq(1,1,0))

    // grouped BAM should contain all the reads
    {
      val groupedRecs = readBamRecs(paths.groupedBam)
      groupedRecs.length shouldBe this.inputRecords.length
      groupedRecs.map(toDuplexName).toSet should contain theSameElementsAs this.inputRecords.map(toDuplexName).toSet
      val actualMis = groupedRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("0/A", "1/A", "1/B", "2/A", "2/B", "3/A", "3/B", "4/A", "4/B", "5/A", "5/B")
      actualMis should contain theSameElementsInOrderAs expectedMis
    }

    // duplex consensuses from all groups
    {
      val consensusRecs = readBamRecs(paths.consensusBam)
      consensusRecs.length shouldBe 12
      val actualMis = consensusRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("0", "1", "2", "3", "4", "5")
      actualMis should contain theSameElementsInOrderAs expectedMis
      consensusRecs.map(_.name).distinct should contain theSameElementsInOrderAs Seq("A:0", "A:1", "A:2", "A:3", "A:4", "A:5")
      // duplex 5 should have Ns
      consensusRecs.filter(_.apply[String]("MI") == "5").foreach { rec => rec.basesString shouldBe "N" * 100 }
    }

    // duplex consensuses from all groups except duplex 5
    {
      val filteredRecs = readBamRecs(paths.filteredBam)
      filteredRecs.length shouldBe 10
      val actualMis = filteredRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("0", "1", "2", "3", "4")
      actualMis should contain theSameElementsInOrderAs expectedMis
      filteredRecs.map(_.name).distinct should contain theSameElementsInOrderAs Seq("A:0", "A:1", "A:2", "A:3", "A:4")
    }
  }

  it should "support filtering based on min-reads at both the calling and filtering step" in {
    val paths = call(minReadsForCalling=Seq(2,1,1), minReadsForFiltering=Seq(9,5,4))

    // grouped BAM should contain all the reads
    {
      val groupedRecs = readBamRecs(paths.groupedBam)
      groupedRecs.length shouldBe this.inputRecords.length
      groupedRecs.map(toDuplexName).toSet should contain theSameElementsAs this.inputRecords.map(toDuplexName).toSet
      val actualMis = groupedRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("0/A", "1/A", "1/B", "2/A", "2/B", "3/A", "3/B", "4/A", "4/B", "5/A", "5/B")
      actualMis should contain theSameElementsInOrderAs expectedMis
    }

    // duplex consensuses from all groups with both strands (excludes duplex 0), so groups 1-5
    {
      val consensusRecs = readBamRecs(paths.consensusBam)
      consensusRecs.length shouldBe 10
      val actualMis = consensusRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("1", "2", "3", "4", "5")
      actualMis should contain theSameElementsInOrderAs expectedMis
      consensusRecs.map(_.name).distinct should contain theSameElementsInOrderAs Seq("A:1", "A:2", "A:3", "A:4", "A:5")
      // duplex 5 should have Ns
      consensusRecs.filter(_.apply[String]("MI") == "5").foreach { rec => rec.basesString shouldBe "N" * 100 }
    }

    // duplex consensuses from groups 3 and 5
    {
      val filteredRecs = readBamRecs(paths.filteredBam)
      filteredRecs.length shouldBe 4
      val actualMis = filteredRecs.map(_.apply[String]("MI").toString).distinct
      val expectedMis = Seq("3", "4")
      actualMis should contain theSameElementsInOrderAs expectedMis
      filteredRecs.map(_.name).distinct should contain theSameElementsInOrderAs Seq("A:3", "A:4")
    }
  }
}
