package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import htsjdk.samtools.util.{Interval, IntervalList}

import scala.jdk.CollectionConverters.SeqHasAsJava

class CoverageManagerTest extends UnitSpec {

  private val dict = SequenceDictionary(
    SequenceMetadata(name = "chr1", length = 10000),
    SequenceMetadata(name = "chr2", length = 50000)
  ).asSam
  val chr1Index: Int = dict.getSequenceIndex("chr1")
  val chr2Index: Int = dict.getSequenceIndex("chr2")

  val i1: Interval = new Interval("chr1", 100, 200)
  val i2: Interval = new Interval("chr1", 100, 300)
  val i3: Interval = new Interval("chr1", 300, 400)

  val i4: Interval = new Interval("chr2", 300, 400)

  val il: IntervalList = new IntervalList(dict)
  il.addall(List(i1, i2, i3).asJava)
  val covManTwo: CoverageManager = new CoverageManager(il, minMapQ = 30.toByte, 2)
  val covManOne: CoverageManager = new CoverageManager(il, minMapQ = 30.toByte, 1)
  val covManZero: CoverageManager = new CoverageManager(il, minMapQ = 10.toByte, 0)

  "CoverageManager.needsCoverage" should "start being true" in {
    covManOne.needsCoverage(i1) shouldBe true
    covManOne.needsCoverage(i2) shouldBe true
    covManOne.needsCoverage(i3) shouldBe true
  }

  it should "return false when out of range" in {
    covManOne.needsCoverage(i4) shouldBe false
  }

  it should "remain true when partially covered" in {
    covManOne.resetCoverage()
    covManOne.incrementCoverage(i1)
    covManOne.needsCoverage(i2) shouldBe true
    covManOne.needsCoverage(i3) shouldBe true
    covManOne.needsCoverage(i1) shouldBe false
  }

  it should "remain false when out of range" in {
    covManOne.resetCoverage()
    covManOne.incrementCoverage(i1)
    covManOne.needsCoverage(i4) shouldBe false
  }

  it should "get to 2 when adding coverage twice" in {
    covManTwo.resetCoverage()
    covManOne.incrementCoverage(i1)
    covManOne.incrementCoverage(i1)
    covManOne.needsCoverage(i2) shouldBe true
    covManOne.needsCoverage(i3) shouldBe true
    covManOne.needsCoverage(i1) shouldBe false
  }

  it should "increment its coverage by 1 correctly from a template with overlapping reads" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150)
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    covManTwo.resetCoverage()
    covManTwo.processTemplates(templates.iterator).isEmpty shouldBe false

    covManTwo.needsCoverage(new Interval("chr1", 150, 200)) shouldBe true
  }

  it should "Correctly add reads from an easy template" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150)
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    // here we should *not* add anything to the coverage, since the coverage target is zero
    covManZero.resetCoverage()

    covManZero.processTemplates(templates.iterator).isEmpty shouldBe true
    covManZero.needsCoverage(new Interval("chr1", 150, 200)) shouldBe false

    // here we *should* add to the coverage, since the coverage target is 10
    covManOne.resetCoverage()

    covManOne.processTemplates(templates.iterator).isEmpty shouldBe false
    covManOne.needsCoverage(new Interval("chr1", 150, 200)) shouldEqual false
  }

  it should "Correctly only consider coverage from the reads of high mapq (target zero)" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150, mapq1 = PhredScore.MinValue, mapq2 = PhredScore.MaxValue)

    // Here we should *not* add anything to the coverage, since the coverage target is zero.
    covManZero.resetCoverage()
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq
    covManZero.processTemplates(templates.iterator).isEmpty shouldBe true
    covManZero.needsCoverage(new Interval("chr1", 150, 200)) shouldBe false
  }

  it should "Correctly only consider coverage from the reads of high mapq (target 10)" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150, mapq1 = PhredScore.MinValue, mapq2 = PhredScore.MaxValue)

    // Here we *should* add to the coverage, since the coverage target is 10.
    covManOne.resetCoverage()
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq
    covManOne.processTemplates(templates.iterator).isEmpty shouldBe false
    covManOne.needsCoverage(new Interval("chr1", 100, 149)) shouldBe true
    covManOne.needsCoverage(new Interval("chr1", 150, 150))  shouldBe false
    covManOne.needsCoverage(new Interval("chr1", 0, 250)) shouldBe true
  }

  it should "Correctly not add reads that should be filtered " in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)

    builder.addPair(name = "low_mapq", contig = chr1Index, start1 = 100, start2 = 150, mapq1 = 10, mapq2 = 10)

    builder.addPair(name = "secondary", contig = chr1Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.secondary = true)

    builder.addPair(name = "duplicates", contig = chr1Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.duplicate = true)

    builder.addPair(name = "non_pf", contig = chr1Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.pf = false)

    builder.addPair(name = "non-overlapping1", contig = chr1Index, start1 = 400, start2 = 550)
      .iterator.foreach(r => r.pf = false)

    builder.addPair(name = "non-overlapping2", contig = chr2Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.pf = false)

    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    covManOne.resetCoverage()

    // here we *should not* add to the coverage, since non of the reads pass filters
    covManOne.processTemplates(templates.iterator).isEmpty shouldBe true
    covManOne.needsCoverage(new Interval("chr1", 150, 200)) shouldBe true
  }

  it should "Correctly not add reads that should be filtered but include good read from same template" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)

    builder.addPair(name = "secondary", contig = chr1Index, start1 = 300, start2 = 400)
    builder.addPair(name = "secondary", contig = chr2Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.secondary = true)

    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    covManOne.resetCoverage()

    // here we *should not* add to the coverage, since none of the reads pass filters
    val filteredTemplates = covManOne.processTemplates(templates.iterator)

    filteredTemplates.isEmpty shouldBe false
    filteredTemplates.size shouldBe 1
    covManOne.needsCoverage(new Interval("chr1", 100, 200)) shouldBe true
    covManOne.needsCoverage(new Interval("chr1", 300, 400)) shouldBe false
  }
}
