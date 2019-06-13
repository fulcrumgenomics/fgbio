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

import java.nio.file.{Path, Paths}
import java.util.Random

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Io, Metric, Rscript}
import htsjdk.samtools.util.{Interval, IntervalList}
import org.apache.commons.math3.distribution.NormalDistribution

import scala.math.{max, min}

class CollectDuplexSeqMetricsTest extends UnitSpec {
  private val MI = "MI"
  private val RX = "RX"

  // Case class to hold pointers to all the outputs of CollectDuplexSeqMetrics
  private case class Outputs(families: Path, duplexFamilies: Path, umis: Path, duplexUmis: Path, yields: Path, plots: Path) {
    import CollectDuplexSeqMetrics._
    lazy val familyMetrics:       Seq[FamilySizeMetric]       = Metric.read[FamilySizeMetric](families)
    lazy val duplexFamilyMetrics: Seq[DuplexFamilySizeMetric] = Metric.read[DuplexFamilySizeMetric](duplexFamilies)
    lazy val umiMetrics:          Seq[UmiMetric]              = Metric.read[UmiMetric](umis)
    lazy val duplexUmiMetrics:    Seq[DuplexUmiMetric]        = if (duplexUmis.toFile.exists()) Metric.read[DuplexUmiMetric](duplexUmis) else Seq.empty
    lazy val yieldMetrics:        Seq[DuplexYieldMetric]      = Metric.read[DuplexYieldMetric](yields)
  }

  // Executes duplex seq metrics and returns paths to the various metrics
  private def exec(builder: SamBuilder, intervals: Option[PathToIntervals] = None, ab: Int = 1, ba: Int = 1, plot: Boolean=false, duplexCounts: Boolean = false): Outputs = {
    val bam = builder.toTempFile()
    val prefix = makeTempFile("duplex.", ".output")
    new CollectDuplexSeqMetrics(input=bam, output=prefix, intervals=intervals, duplexUmiCounts=duplexCounts, minAbReads=ab, minBaReads=ba, generatePlots=plot).execute()
    Outputs(
      families       = Paths.get(prefix.toString + CollectDuplexSeqMetrics.FamilySizeMetricsExt),
      duplexFamilies = Paths.get(prefix.toString + CollectDuplexSeqMetrics.DuplexFamilySizeMetricsExt),
      umis           = Paths.get(prefix.toString + CollectDuplexSeqMetrics.UmiMetricsExt),
      duplexUmis     = Paths.get(prefix.toString + CollectDuplexSeqMetrics.DuplexUmiMetricsExt),
      yields         = Paths.get(prefix.toString + CollectDuplexSeqMetrics.YieldMetricsExt),
      plots          = Paths.get(prefix.toString + CollectDuplexSeqMetrics.PlotsExt)
    )
  }

  // Returns a collector as an option for easy mapping over
  private def collector(duplex: Boolean = false) = Some(new CollectDuplexSeqMetrics(input=Io.DevNull, output=Io.DevNull, duplexUmiCounts=duplex))

  "CollectDuplexSeqMetrics" should "have acceptable CLP annotations" in {
    checkClpAnnotations[CollectDuplexSeqMetrics]
  }

  it should "count UMIs once per read-pair" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "TTT-AAA", MI -> "1/B"))
    val results = exec(builder, duplexCounts=false)

    results.umiMetrics.size shouldBe 2
    results.umiMetrics.foreach { m =>
      m.raw_observations shouldBe 2
      m.unique_observations shouldBe 1
    }

    results.duplexUmiMetrics should have size 0
  }

  it should "error-correct UMIs for counting purposes" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(start1=100, start2=200, strand1=Plus , strand2=Minus, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, strand1=Plus , strand2=Minus, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, strand1=Plus , strand2=Minus, attrs=Map(RX -> "CAA-TTT", MI -> "1/A"))
    builder.addPair(start1=200, start2=100, strand1=Minus, strand2=Plus,  attrs=Map(RX -> "TTT-AAA", MI -> "1/B"))
    builder.addPair(start1=200, start2=100, strand1=Minus, strand2=Plus,  attrs=Map(RX -> "TTT-AAA", MI -> "1/B"))
    builder.addPair(start1=200, start2=100, strand1=Minus, strand2=Plus,  attrs=Map(RX -> "CTT-AAA", MI -> "1/B"))
    val results = exec(builder, duplexCounts=true)

    results.umiMetrics.size shouldBe 2
    results.umiMetrics.foreach { m =>
      m.umi == "AAA" || m.umi == "TTT" shouldBe true
      m.raw_observations shouldBe 6
      m.raw_observations_with_errors shouldBe 1
      m.unique_observations shouldBe 1
    }

    results.duplexUmiMetrics should have size 1
    results.duplexUmiMetrics.head.umi shouldBe "AAA-TTT"
    results.duplexUmiMetrics.head.raw_observations    shouldBe 6
    results.duplexUmiMetrics.head.unique_observations shouldBe 1
    results.duplexUmiMetrics.head.fraction_unique_observations shouldBe 1.0
    results.duplexUmiMetrics.head.fraction_unique_observations_expected shouldBe 0.25
  }

  it should "count unique UMI observations 1-to-1 with tag families" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "TTT-AAA", MI -> "1/B"))

    builder.addPair(start1=150, start2=250, attrs=Map(RX -> "TTT-AAA", MI -> "2/A"))
    builder.addPair(start1=150, start2=250, attrs=Map(RX -> "TTT-AAA", MI -> "2/A"))
    builder.addPair(start1=150, start2=250, attrs=Map(RX -> "NTT-AAA", MI -> "2/A"))

    builder.addPair(start1=250, start2=350, attrs=Map(RX -> "CCC-GGG", MI -> "3/B"))

    val metrics = exec(builder).umiMetrics

    metrics.size shouldBe 4
    val aaa = metrics.find(_.umi == "AAA").get
    aaa.raw_observations shouldBe 5
    aaa.raw_observations_with_errors shouldBe 0
    aaa.unique_observations shouldBe 2

    val ttt = metrics.find(_.umi == "TTT").get
    ttt.raw_observations shouldBe 5
    ttt.raw_observations_with_errors shouldBe 1
    ttt.unique_observations shouldBe 2

    val ccc = metrics.find(_.umi == "CCC").get
    ccc.raw_observations shouldBe 1
    ccc.raw_observations_with_errors shouldBe 0
    ccc.unique_observations shouldBe 1

    val ggg = metrics.find(_.umi == "GGG").get
    ggg.raw_observations shouldBe 1
    ggg.raw_observations_with_errors shouldBe 0
    ggg.unique_observations shouldBe 1
  }

  it should "generate accurate family size counts" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    // Three SS families at the same location
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-TTT", MI -> "1/A"))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "ACG-GGA", MI -> "2/A"))
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "TAT-CGT", MI -> "3/B"))

    // Two duplex tag families at the same location
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "TTT-AAA", MI -> "4/A"))
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "AAA-AAA", MI -> "4/B"))
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "CCC-GGG", MI -> "5/A"))
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "GGG-CCC", MI -> "5/B"))

    // One duplex and one SS family at the same location
    builder.addPair(start1=400, start2=500, attrs=Map(RX -> "GCG-GAA", MI -> "6/A"))
    builder.addPair(start1=400, start2=500, attrs=Map(RX -> "ACG-CCT", MI -> "7/A"))
    builder.addPair(start1=400, start2=500, attrs=Map(RX -> "ACG-CCT", MI -> "7/A"))
    builder.addPair(start1=400, start2=500, attrs=Map(RX -> "CCT-ACG", MI -> "7/B"))

    // Resulting families:
    //   cs: three @ size=4,
    //   ss: eight @ size=1, two @ size=2
    //   ds: three @ size=1, three @ size = 2, one @ size=3

    val metrics = exec(builder).familyMetrics
    metrics should have size 4
    val one = metrics.find(_.family_size == 1).get
    one.cs_count shouldBe 0
    one.ss_count shouldBe 8
    one.ds_count shouldBe 3

    val two = metrics.find(_.family_size == 2).get
    two.cs_count shouldBe 0
    two.ss_count shouldBe 2
    two.ds_count shouldBe 3

    val three = metrics.find(_.family_size == 3).get
    three.cs_count shouldBe 0
    three.ss_count shouldBe 0
    three.ds_count shouldBe 1

    val four = metrics.find(_.family_size == 4).get
    four.cs_count shouldBe 3
    four.ss_count shouldBe 0
    four.ds_count shouldBe 0
  }

  it should "generate accurate duplex family size counts" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))

    // 1/0
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-ACG", MI -> "1/A"))

    // 1/1
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "AAA-ACG", MI -> "2/A"))
    builder.addPair(start1=200, start2=300, attrs=Map(RX -> "ACG-AAA", MI -> "2/B"))

    // 2/1
    (1 to 1).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "GGG-AAC", MI -> "3/A")) }
    (1 to 2).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "AAC-GGG", MI -> "3/B")) }

    // 4/3
    (1 to 4).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "GGG-AAC", MI -> "4/A")) }
    (1 to 3).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "AAC-GGG", MI -> "4/B")) }

    // 7/2
    (1 to 2).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "CAC-TGG", MI -> "5/A")) }
    (1 to 7).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "TGG-CAC", MI -> "5/B")) }

    // 5/5
    (1 to 5).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "AGT-GCT", MI -> "6/A")) }
    (1 to 5).foreach { _ => builder.addPair(start1=400, start2=500, attrs=Map(RX -> "GCT-AGT", MI -> "6/B")) }

    // 2/1
    (1 to 1).foreach { _ => builder.addPair(start1=500, start2=600, attrs=Map(RX -> "GGG-AAC", MI -> "7/A")) }
    (1 to 2).foreach { _ => builder.addPair(start1=500, start2=600, attrs=Map(RX -> "AAC-GGG", MI -> "7/B")) }

    val metrics = exec(builder).duplexFamilyMetrics
    metrics should have size 6
    metrics.forall { m => m.ab_size >= m.ba_size } shouldBe true
    metrics.find(m => m.ab_size == 1 && m.ba_size == 0).get.count shouldBe 1
    metrics.find(m => m.ab_size == 1 && m.ba_size == 1).get.count shouldBe 1
    metrics.find(m => m.ab_size == 2 && m.ba_size == 1).get.count shouldBe 2
    metrics.find(m => m.ab_size == 4 && m.ba_size == 3).get.count shouldBe 1
    metrics.find(m => m.ab_size == 7 && m.ba_size == 2).get.count shouldBe 1
    metrics.find(m => m.ab_size == 5 && m.ba_size == 5).get.count shouldBe 1
  }

  it should "respect the min-ab and min-ba for counting duplexes" in {
    val builder = new SamBuilder(readLength=10, sort=Some(SamOrder.TemplateCoordinate))
    (1 to 1).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "AAA-GGG", MI -> "1/A")) }
    (1 to 1).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "GGG-AAA", MI -> "1/B")) }
    (1 to 1).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "ACT-TTA", MI -> "2/A")) }
    (1 to 2).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "TTA-ACT", MI -> "2/B")) }
    (1 to 2).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "CGA-GGT", MI -> "3/A")) }
    (1 to 2).foreach { _ => builder.addPair(start1=300, start2=400, attrs=Map(RX -> "GGT-CGA", MI -> "3/B")) }

    exec(builder, ab=1, ba=1).yieldMetrics.find(_.fraction == 1).get.ds_duplexes shouldBe 3
    exec(builder, ab=2, ba=1).yieldMetrics.find(_.fraction == 1).get.ds_duplexes shouldBe 2
    exec(builder, ab=2, ba=2).yieldMetrics.find(_.fraction == 1).get.ds_duplexes shouldBe 1
  }

  it should "run end to end and generate plots" in {
    val random = new Random(42)
    val abDist = new NormalDistribution(6, 2)
    val baDist = new NormalDistribution(3, 0.5)
    val isize  = new NormalDistribution(150, 10)
    val builder = new SamBuilder(readLength=75)
    def sampleNatural(dist: NormalDistribution): Int = { val d = dist.sample(); if (d < 0) 0 else d.round.toInt }
    val counter = new SimpleCounter[(Int,Int)]

    (1 to 5000).foreach { i =>
      val start1 = random.nextInt(5000) + 1
      val start2 = start1 + isize.sample().toInt - 75
      val abs = sampleNatural(abDist)
      val bas = sampleNatural(baDist)
      (1 to abs) foreach { _ => builder.addPair(start1=start1, start2=start2, attrs=Map(RX -> "AAA-TTT", MI -> (i + "/A"))) }
      (1 to bas) foreach { _ => builder.addPair(start1=start1, start2=start2, attrs=Map(RX -> "TTT-AAA", MI -> (i + "/B"))) }
      counter.count((max(abs, bas), min(abs, bas)))
    }

    val metrics = exec(builder, plot=Rscript.Available)
    metrics.plots.toFile.exists() shouldBe Rscript.Available

    // Check that the duplex family sizes tie out
    metrics.duplexFamilyMetrics.foreach { m => m.count shouldBe counter.countOf((m.ab_size, m.ba_size)) }

    // Check that the yield metrics are sensible
    metrics.yieldMetrics.grouped(2).foreach { case Seq(m1, m2) =>
        m1.fraction    <  m2.fraction shouldBe true
        m1.cs_families <= m2.cs_families shouldBe true
        m1.ss_families <= m2.ss_families shouldBe true
        m1.ds_families <= m2.ds_families shouldBe true
        m1.ds_duplexes <= m2.ds_duplexes shouldBe true
    }
  }

  it should "only count inserts overlapping one or more intervals when intervals are provided" in {
    val builder = new SamBuilder(readLength=100, sort=Some(SamOrder.TemplateCoordinate))
    (1 to 1).foreach { _ => builder.addPair(contig=0, start1=1000, start2=1100, attrs=Map(RX -> "AAA-GGG", MI -> "1/A")) }
    (1 to 2).foreach { _ => builder.addPair(contig=0, start1=2000, start2=2100, attrs=Map(RX -> "GGG-AAA", MI -> "2/A")) }
    (1 to 3).foreach { _ => builder.addPair(contig=0, start1=3000, start2=3100, attrs=Map(RX -> "ACT-TTA", MI -> "3/A")) }
    (1 to 4).foreach { _ => builder.addPair(contig=1, start1=4000, start2=4100, attrs=Map(RX -> "TTA-ACT", MI -> "4/A")) }
    (1 to 5).foreach { _ => builder.addPair(contig=1, start1=5000, start2=5100, attrs=Map(RX -> "CGA-GGT", MI -> "5/A")) }
    (1 to 6).foreach { _ => builder.addPair(contig=1, start1=6000, start2=6100, attrs=Map(RX -> "GGT-CGA", MI -> "6/A")) }

    val intervals = new IntervalList(builder.header)
    intervals.add(new Interval(builder.dict.getSequence(0).getSequenceName, 900,  1001)) // Captures the first insert by overlapping R1
    intervals.add(new Interval(builder.dict.getSequence(0).getSequenceName, 3150, 3500)) // Captures the third insert by overlapping R1
    intervals.add(new Interval(builder.dict.getSequence(1).getSequenceName, 5050, 6050)) // Captures the fifth and sixth inserts
    val path = makeTempFile("targets.", ".interval_list")
    intervals.write(path.toFile)

    val duplexFamilies = exec(builder, intervals=Some(path)).duplexFamilyMetrics
    duplexFamilies.size shouldBe 4
    duplexFamilies.find(f => f.ab_size == 1 && f.ba_size == 0).get.count shouldBe 1
    duplexFamilies.find(f => f.ab_size == 3 && f.ba_size == 0).get.count shouldBe 1
    duplexFamilies.find(f => f.ab_size == 5 && f.ba_size == 0).get.count shouldBe 1
    duplexFamilies.find(f => f.ab_size == 6 && f.ba_size == 0).get.count shouldBe 1
  }

  "CollectDuplexSeqMetrics.updateUmiMetrics" should "not count duplex umis" in collector(duplex=false).foreach { c =>
    val builder = new SamBuilder(readLength=10)
    builder.addPair(start1=100, start2=200, attrs=Map(RX -> "AAA-CCC", MI -> "1/A"))
    c.updateUmiMetrics(Seq(builder.toSeq))
    val metrics = c.duplexUmiMetrics(c.umiMetrics)
    metrics should have size 0
  }

  it should "count UMIs as if on F1R2 molecules" in collector(duplex=true).foreach { c =>
    val builder = new SamBuilder(readLength=10)
    builder.addPair(start1=100, start2=200, strand1=Plus,  strand2=Minus, attrs=Map(RX -> "AAA-CCC", MI -> "1/A"))
    builder.addPair(start1=200, start2=100, strand1=Minus, strand2=Plus,  attrs=Map(RX -> "CCC-AAA", MI -> "1/B"))
    builder.addPair(start1=300, start2=400, strand1=Plus,  strand2=Minus, attrs=Map(RX -> "CCC-GGG", MI -> "2/A"))
    builder.addPair(start1=900, start2=800, strand1=Minus, strand2=Plus,  attrs=Map(RX -> "TTT-AAA", MI -> "3/A"))

    builder.toSeq.groupBy(r => r[String](MI).takeWhile(_ != '/')).values
      .map(rs => rs.groupBy(r => r[String](MI)).values.toSeq).toSeq
      .foreach(group => c.updateUmiMetrics(group))

    val metrics = c.duplexUmiMetrics(c.umiMetrics)
    metrics should have size 3

    metrics.find(_.umi == "AAA-CCC") shouldBe defined
    metrics.find(_.umi == "AAA-CCC").foreach(m => m.unique_observations shouldBe 1)

    metrics.find(_.umi == "CCC-GGG") shouldBe defined
    metrics.find(_.umi == "CCC-GGG").foreach(m => m.unique_observations shouldBe 1)

    metrics.find(_.umi == "AAA-TTT") shouldBe defined
    metrics.find(_.umi == "AAA-TTT").foreach(m => m.unique_observations shouldBe 1)
  }
}
