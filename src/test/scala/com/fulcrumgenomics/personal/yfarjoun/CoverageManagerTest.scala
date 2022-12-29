package com.fulcrumgenomics.personal.yfarjoun

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
  val covMan: CoverageManager = new CoverageManager(il, minMapQ = 30.toByte, 10)
  val covManZero: CoverageManager = new CoverageManager(il, minMapQ = 10.toByte, 0)

  "CoverageManager.getMinCoverage" should "start its life with coverage = zero" in {
    covMan.getMinCoverage(i1) shouldEqual 0
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
  }

  it should "return Short.MaxValue when out of range" in {
    covMan.getMinCoverage(i4) shouldEqual Short.MaxValue
  }

  it should "remain coverage zero when only partially covered" in {
    covMan.resetCoverage()
    covMan.incrementCoverage(i1)
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
    covMan.getMinCoverage(i1) shouldEqual 1
  }

  it should "remain Short.MaxValue when out of range" in {
    covMan.getMinCoverage(i4) shouldEqual Short.MaxValue
  }

  it should "get to 2 when adding coverage twice" in {
    covMan.resetCoverage()
    covMan.incrementCoverage(i1)
    covMan.incrementCoverage(i1)
    covMan.getMinCoverage(i2) shouldEqual 0
    covMan.getMinCoverage(i3) shouldEqual 0
    covMan.getMinCoverage(i1) shouldEqual 2
  }

  it should "increment its coverage by 1 correctly from a template with overlapping reads" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150)
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    covMan.resetCoverage()
    templates.foreach(t => covMan.incrementCoverage(t))

    covMan.getMinCoverage(new Interval("chr1", 150, 200)) shouldEqual 1
    covMan.getMinTemplateCoverage(template = templates.head) shouldEqual 1
  }

  it should "Correctly add reads from an easy template" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150)
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    // here we should *not* add anything to the coverage, since the coverage target is zero
    covManZero.resetCoverage()

    covManZero.processTemplates(templates.iterator).isEmpty shouldBe true
    covManZero.getMinTemplateCoverage(template = templates.head) shouldEqual 0
    covManZero.getMinCoverage(new Interval("chr1", 150, 200)) shouldEqual 0

    // here we *should* add to the coverage, since the coverage target is 10
    covMan.resetCoverage()

    covMan.processTemplates(templates.iterator).isEmpty shouldBe false
    covMan.getMinCoverage(new Interval("chr1", 150, 200)) shouldEqual 1
    covMan.getMinTemplateCoverage(template = templates.head) shouldEqual 1
  }

  it should "Correctly only consider coverage from the reads of high mapq (target zero)" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150, mapq1 = PhredScore.MinValue, mapq2 = PhredScore.MaxValue)

    // Here we should *not* add anything to the coverage, since the coverage target is zero.
    covManZero.resetCoverage()
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq
    covManZero.processTemplates(templates.iterator).isEmpty shouldBe true
    covManZero.getMinCoverage(new Interval("chr1", 150, 200)) shouldEqual 0
    covManZero.getMinTemplateCoverage(template = templates.head) shouldEqual 0
  }

  it should "Correctly only consider coverage from the reads of high mapq (target 10)" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)
    builder.addPair(name = "p1", contig = chr1Index, start1 = 100, start2 = 150, mapq1 = PhredScore.MinValue, mapq2 = PhredScore.MaxValue)

    // Here we *should* add to the coverage, since the coverage target is 10.
    covMan.resetCoverage()
    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq
    covMan.processTemplates(templates.iterator).isEmpty shouldBe false
    covMan.getMinCoverage(new Interval("chr1", 100, 149)) shouldEqual 0
    covMan.getMinCoverage(new Interval("chr1", 150, 150)) shouldEqual 1
    covMan.getMinTemplateCoverage(template = templates.head) shouldEqual 1
    covMan.getMinCoverage(new Interval("chr1", 0, 250)) shouldEqual 0

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

    covMan.resetCoverage()

    // here we *should not* add to the coverage, since non of the reads pass filteres
    covMan.processTemplates(templates.iterator).isEmpty shouldBe true
    covMan.getMinCoverage(new Interval("chr1", 150, 200)) shouldEqual 0.toShort
    // since there's no region of interest as all the reads get filtered out...
    covMan.getMinTemplateCoverage(template = templates.head) shouldEqual Short.MaxValue
  }


  it should "Correctly not add reads that should be filtered but include good read from same template" in {
    val builder = new SamBuilder(sort = Some(SamOrder.Queryname))
    builder.header.setSequenceDictionary(dict)

    builder.addPair(name = "secondary", contig = chr1Index, start1 = 300, start2 = 400)

    builder.addPair(name = "secondary", contig = chr2Index, start1 = 100, start2 = 150)
      .iterator.foreach(r => r.secondary = true)

    val templates: Seq[Template] = Bams.templateIterator(builder.toSource).toSeq

    covMan.resetCoverage()

    // here we *should not* add to the coverage, since none of the reads pass filters
    val filteredTemplates = covMan.processTemplates(templates.iterator)

    filteredTemplates.isEmpty shouldBe false
    filteredTemplates.size shouldBe 1
    covMan.getMinCoverage(new Interval("chr1", 100, 200)) shouldEqual 0
    covMan.getMinCoverage(new Interval("chr1", 300, 400)) shouldEqual 1
    covMan.getMinTemplateCoverage(template = templates.head) shouldEqual 1
  }
}
