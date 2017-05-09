/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf

import java.nio.file.{Files, Paths}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.{ErrorLogLevel, UnitSpec, VariantContextSetBuilder}
import com.fulcrumgenomics.util.Metric
import com.fulcrumgenomics.vcf.PhaseCigar.IlluminaSwitchErrors
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.NumericCounter
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.{Interval, IntervalList}
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.variantcontext.{GenotypeBuilder, VariantContext, VariantContextBuilder}
import htsjdk.variant.vcf.{VCFFileReader, VCFHeader}

object AssessPhasingTest {
  def withPhasingSetId(ctx: VariantContext, id: Int): VariantContext = {
    val gBuilder = new GenotypeBuilder(ctx.getGenotype(0))
    gBuilder.attribute("PS", id)
    val ctxBuilder = new VariantContextBuilder(ctx)
    ctxBuilder.genotypes(gBuilder.make())
    ctxBuilder.make()
  }

  private val builderTruth = new VariantContextSetBuilder()
  private val builderCall  = new VariantContextSetBuilder()

  // init builderTruth and builderCall
  {
    // BLOCK #1: positions 1 - 4 (call 1-4, truth 1-4)
    builderTruth.addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Match
    builderCall.addVariant( start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) //   - with previous
    builderTruth.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // TruthOnly
    builderCall.addVariant( start=3, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // CallOnly
    builderTruth.addVariant(start=4, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Match
    builderCall.addVariant( start=4, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) //   - with previous

    // BLOCK #2: positions 11-16 (call 11-15, truth 11-16)
    builderTruth.addVariant(start=11, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Match
    builderCall.addVariant( start=11, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) //   - with previous
    builderTruth.addVariant(start=12, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Mismatch **** POINT ERROR ****
    builderCall.addVariant( start=12, variantAlleles=List("A", "C"), genotypeAlleles=List("C", "A"), phased=true) //   - with previous
    builderTruth.addVariant(start=13, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // NA
    builderCall.addVariant( start=13, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // CallOnly
    builderTruth.addVariant(start=14, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // TruthOnly
    builderCall.addVariant( start=14, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // NA
    builderTruth.addVariant(start=15, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Match
    builderCall.addVariant( start=15, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) //   - with Previous
    builderTruth.addVariant(start=16, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // TruthOnly

    // BLOCK #3: position 21 (call 21-21)
    builderCall.addVariant( start=21, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // CallOnly

    // BLOCK #4: position 30-42
    Range(30, 37).foreach { start =>
      builderTruth.addVariant(start=start, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Match
      builderCall.addVariant( start=start, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) //   - with previous
    }
    Range(37, 43).foreach { start =>
      builderTruth.addVariant(start=start, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true) // Mismatch
      builderCall.addVariant( start=start, variantAlleles=List("A", "C"), genotypeAlleles=List("C", "A"), phased=true) //   - with previous
    }

    // NB: call has blocks lengths 4, 5, 1, and 13; truth has block lengths 4, 6, and 13.
  }

  private def addPhaseSetId(ctx: VariantContext): VariantContext = {
    if (ctx.getStart <= 10) withPhasingSetId(ctx, 1)
    else if (ctx.getStart <= 20) withPhasingSetId(ctx, 11)
    else if (ctx.getStart <= 29) withPhasingSetId(ctx, 21)
    else if (ctx.getStart <= 42) withPhasingSetId(ctx, 30)
    else unreachable("Not defined")
  }

  lazy val TruthVariants: Seq[VariantContext] = {
    // Keep the truth variant position 13 without a phase set
    builderTruth.map { ctx =>
      if (ctx.getStart == 13) ctx
      else addPhaseSetId(ctx)
    }.toSeq
  }

  lazy val CallVariants: Seq[VariantContext] = {
    // Keep the call variant position 14 without a phase set
    builderCall.map { ctx =>
      if (ctx.getStart == 14) ctx
      else addPhaseSetId(ctx)
    }.toSeq
  }

  val Header = builderCall.header
}

/**
  * Tests for AssessPhasing
  */
class AssessPhasingTest extends ErrorLogLevel {
  import AssessPhasingTest._

  private def writeVariants(variants: Seq[VariantContext]): PathToVcf = {
    val path = Files.createTempFile("AssessPhasingTest.", ".vcf.gz")
    path.toFile.deleteOnExit()
    val writer = new VariantContextWriterBuilder()
      .setReferenceDictionary(Header.getSequenceDictionary)
      .setOutputFile(path.toFile)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .build()
    writer.writeHeader(Header)
    variants.foreach(writer.add)
    writer.close()
    path
  }

  private def runEndToEndWithFilter(intervals: Option[PathToIntervals],
                                    expectedPrefix: PathPrefix,
                                    start: Int,
                                    end: Int,
                                    debugVcf: Boolean = false): Unit = {
    runEndToEnd(
      intervals      = intervals,
      expectedPrefix = expectedPrefix,
      truthVariants  = TruthVariants.filter { ctx => start <= ctx.getStart && ctx.getStart <= end },
      callVariants   = CallVariants.filter  { ctx => start <= ctx.getStart && ctx.getStart <= end },
      debugVcf       = debugVcf
    )
  }

  private def runEndToEnd(intervals: Option[PathToIntervals],
                          expectedPrefix: PathPrefix,
                          truthVariants: Seq[VariantContext] = TruthVariants,
                          callVariants: Seq[VariantContext] = CallVariants,
                          debugVcf: Boolean = false
                         ): Unit = {
    // input files
    val truthVcf = writeVariants(truthVariants)
    val callVcf  = writeVariants(callVariants)

    // output files
    val output                 = makeTempFile("AssessPhasingTest.", "prefix")
    val phasingMetricsPath     = PathUtil.pathTo(output + AssessPhasingMetric.MetricExtension)
    val blockLengthMetricsPath = PathUtil.pathTo(output + PhaseBlockLengthMetric.MetricExtension)
    Seq(phasingMetricsPath, blockLengthMetricsPath).foreach(_.toFile.deleteOnExit)

    // run it
    new AssessPhasing(truthVcf=truthVcf, calledVcf=callVcf, knownIntervals=intervals, output=output, debugVcf=debugVcf).execute()

    // read the metrics
    val phasingMetrics     = Metric.read[AssessPhasingMetric](path=phasingMetricsPath)
    val blockLengthMetrics = Metric.read[PhaseBlockLengthMetric](path=blockLengthMetricsPath)

    // get the expected metrics paths
    val expectedPhasingMetricsPath     = PathUtil.pathTo(expectedPrefix + AssessPhasingMetric.MetricExtension)
    val expectedBlockLengthMetricsPath = PathUtil.pathTo(expectedPrefix + PhaseBlockLengthMetric.MetricExtension)

    // read the expected metrics
    val expectedPhasingMetrics     = Metric.read[AssessPhasingMetric](path=expectedPhasingMetricsPath)
    val expectedBlockLengthMetrics = Metric.read[PhaseBlockLengthMetric](path=expectedBlockLengthMetricsPath)

    // compare the metrics
    phasingMetrics should contain theSameElementsInOrderAs expectedPhasingMetrics
    blockLengthMetrics should contain theSameElementsInOrderAs expectedBlockLengthMetrics

    if (debugVcf) {
      val annotatedVcf         = Paths.get(output.toString + AssessPhasing.AnnotatedVcfExtension)
      val expectedAnnotatedVcf = PathUtil.pathTo(expectedPrefix + AssessPhasing.AnnotatedVcfExtension)
      val actualContexts = new VCFFileReader(annotatedVcf.toFile, false).iterator().toIterator.toList
      val expectedContexts = new VCFFileReader(expectedAnnotatedVcf.toFile, false).iterator().toIterator.toList
      actualContexts.length shouldBe expectedContexts.length
      actualContexts.zip(expectedContexts).foreach { case (actualCtx, expectedCtx) =>
          actualCtx.toStringDecodeGenotypes shouldBe expectedCtx.toStringDecodeGenotypes
      }
    }
  }

  private def toIntervalList(start: Int, end: Int): PathToIntervals = {
    val path = makeTempFile("AssessPhasingTest.", ".interval_list")
    val contig = Header.getSequenceDictionary.getSequence(0).getSequenceName
    val header = new SAMFileHeader
    header.setSequenceDictionary(Header.getSequenceDictionary)
    val intervalList = new IntervalList(header)
    intervalList.add(new Interval(contig, 1, 10))
    intervalList.write(path.toFile)
    path
  }

  private val dir = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf")

  // Output prefixes for when we *do not* use an interval list
  private val prefix_none_1_10  = dir.resolve("noIntervalList.1-10")
  private val prefix_none_11_20 = dir.resolve("noIntervalList.11-20")
  private val prefix_none_21_21 = dir.resolve("noIntervalList.21-21")
  private val prefix_none_30_42 = dir.resolve("noIntervalList.30-42")
  private val prefix_none_1_20  = dir.resolve("noIntervalList.1-20")
  private val prefix_none_1_21  = dir.resolve("noIntervalList.1-21")
  private val prefix_none_1_29  = dir.resolve("noIntervalList.1-29")
  private val prefix_none_1_42  = dir.resolve("noIntervalList.1-42")

  // Output prefixes for when we *do* use an interval list

  private val prefix_1_10  = dir.resolve("intervalList.1-10")
  private val prefix_11_20 = dir.resolve("intervalList.11-20")
  private val prefix_21_21 = dir.resolve("intervalList.21-21")
  private val prefix_30_42 = dir.resolve("intervalList.30-42")
  private val prefix_1_20  = dir.resolve("intervalList.1-20")
  private val prefix_1_21  = dir.resolve("intervalList.1-21")

  "AssessPhasing" should "run end-to-end on a single phased block" in {
    runEndToEndWithFilter(intervals = None, prefix_none_1_10, 1, 10) // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_11_20, 11, 20)  // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_21_21, 21, 21)  // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_30_42, 30, 42)  // no intervals
  }

  it should "run end-to-end on multiple phased blocks" in {
    runEndToEndWithFilter(intervals = None, prefix_none_1_20, 1, 20)  // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_1_21, 1, 21)  // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_1_29, 1, 29)  // no intervals
    runEndToEndWithFilter(intervals = None, prefix_none_1_42, 1, 42)  // no intervals
  }

  it should "run end-to-end on a single block with an interval list" in {
    runEndToEnd(intervals=Some(toIntervalList(1,  10)), prefix_1_10)  // block 1
    runEndToEnd(intervals=Some(toIntervalList(11, 20)), prefix_11_20) // block 2
    runEndToEnd(intervals=Some(toIntervalList(21, 21)), prefix_21_21) // block 3
    runEndToEnd(intervals=Some(toIntervalList(30, 42)), prefix_30_42) // block 4
  }

  it should "run end-to-end on multiple phased block with an interval list" in {
    runEndToEnd(intervals=Some(toIntervalList(1,  20)), prefix_1_20) // blocks 1-2
    runEndToEnd(intervals=Some(toIntervalList(1,  21)), prefix_1_21) // blocks 1-3
    runEndToEnd(intervals=Some(toIntervalList(1,  42)), prefix_1_21) // blocks 1-4
  }

  it should "produce an annotated VCF when specified" in {
    runEndToEndWithFilter(intervals = None, prefix_none_1_10, 1, 10, debugVcf=true) // no intervals
  }
}

class AssemblyStatisticsTest extends UnitSpec {
  private def forStats(values: Seq[Long]): AssemblyStatistics = {
    val counter = new NumericCounter[Long]()
    values.foreach(counter.count)
    AssemblyStatistics(counter)
  }

  "AssemblyStatistics" should "compute the N50, N90, and L50" in {
    // simple cases
    forStats(Seq(1)) shouldBe AssemblyStatistics(1, 1, 1)
    forStats(Seq(1, 2, 3, 4, 5)) shouldBe AssemblyStatistics(4, 2, 2)
    forStats(Seq(1, 2, 3, 4, 5, 6)) shouldBe AssemblyStatistics(5, 2, 2)
    forStats(Seq(5, 5, 5, 5, 5)) shouldBe AssemblyStatistics(5, 5, 3)

    // a few large blocks can cover n50/n90/l50
    forStats(Seq(1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5)) shouldBe AssemblyStatistics(5, 5, 5)
    forStats(Seq(1, 1, 1, 1, 1, 5, 5, 5, 5, 5)) shouldBe AssemblyStatistics(5, 1, 3)
    forStats(Seq(1, 1, 1, 1, 1, 5, 5, 5)) shouldBe AssemblyStatistics(5, 1, 2)

    // this tests when we get exactly 50% on n50 or l50 that we compute the smallest block length (n50) or the
    // smallest # of blocks (l50).
    forStats(Seq(1, 4, 4, 5, 5)) shouldBe AssemblyStatistics(5, 4, 2) // below
    forStats(Seq(2, 4, 4, 5, 5)) shouldBe AssemblyStatistics(5, 4, 2) // at margin
    forStats(Seq(3, 4, 4, 5, 5)) shouldBe AssemblyStatistics(4, 3, 3) // above
  }

  it should "return empty if no counts were given" in {
    AssemblyStatistics(new NumericCounter[Long]()) shouldBe AssemblyStatistics(0, 0, 0)
  }
}

class PhaseBlockTest extends ErrorLogLevel {
  import AssessPhasingTest.withPhasingSetId

  "PhaseBlock.toOverlapDetector" should "create an empty detector if no variants are given" in {
    val builder = new VariantContextSetBuilder()
    PhaseBlock.buildOverlapDetector(iterator=builder.iterator, dict=builder.header.getSequenceDictionary).getAll shouldBe 'empty
  }

  it should "create an empty detector when variants do not have the phase set tag" in {
    val builder  = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    PhaseBlock.buildOverlapDetector(iterator=builder.iterator, dict=builder.header.getSequenceDictionary).getAll shouldBe 'empty
  }

  it should "create a detector for a single variant" in {
    val builder  = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val iterator = builder.iterator.map { ctx => withPhasingSetId(ctx, 1) }
    val detector = PhaseBlock.buildOverlapDetector(iterator=iterator, dict=builder.header.getSequenceDictionary)
    detector.getAll should have size 1
    val interval = detector.getAll.toSeq.head
    interval.getStart shouldBe 1
    interval.getEnd shouldBe 1
  }

  it should "create a detector from multiple variants within one block" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builder.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builder.addVariant(start=3, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)

    val iterator = builder.iterator.map { ctx => withPhasingSetId(ctx, 1) }
    val detector = PhaseBlock.buildOverlapDetector(iterator=iterator, dict=builder.header.getSequenceDictionary)
    detector.getAll should have size 1
    val interval = detector.getAll.toSeq.head
    interval.getStart shouldBe 1
    interval.getEnd shouldBe 3
  }

  it should "create a detector from multiple variants across multiple blocks" in {
    val builderBlockOne = new VariantContextSetBuilder()
    builderBlockOne.addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=3, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)

    val builderBlockTwo = new VariantContextSetBuilder()
    builderBlockTwo.addVariant(start=4, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockTwo.addVariant(start=5, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockTwo.addVariant(start=6, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)

    val iterator = builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 1) } ++ builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 4) }

    val detector = PhaseBlock.buildOverlapDetector(iterator=iterator, dict=builderBlockOne.header.getSequenceDictionary)
    detector.getAll should have size 2
    val intervals = detector.getAll.toSeq.sortBy(_.getStart)
    val head = intervals.head
    head.getStart shouldBe 1
    head.getEnd shouldBe 3
    val last = intervals.last
    last.getStart shouldBe 4
    last.getEnd shouldBe 6
  }

  it should "keep the larger block when one block encloses/contains another" in {
    val builderBlockOne = new VariantContextSetBuilder()
    builderBlockOne.addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=3, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)

    // second fully contained in the first
    {
      val builderBlockTwo = new VariantContextSetBuilder()
      builderBlockTwo.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)

      val contexts = (builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 1) } ++ builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 2) }).toSeq

      // Check that if we do not want to modify the blocks we get an exception
      an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

      val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
      detector.getAll should have size 1
      val intervals = detector.getAll.toSeq.sortBy(_.getStart)
      val head = intervals.head
      head.getStart shouldBe 1
      head.getEnd shouldBe 3
    }

    // first fully contained in the second
    {
      val builderBlockTwo = new VariantContextSetBuilder()
      builderBlockTwo.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)

      val contexts = (builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 2) } ++ builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 1) }).toSeq

      // Check that if we do not want to modify the blocks we get an exception
      an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

      val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
      detector.getAll should have size 1
      val intervals = detector.getAll.toSeq.sortBy(_.getStart)
      val head = intervals.head
      head.getStart shouldBe 1
      head.getEnd shouldBe 3
    }

    // the same block!
    {
      val contexts = builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 1) }.toSeq

      val detector = PhaseBlock.buildOverlapDetector(iterator = (contexts ++ contexts).iterator, dict = builderBlockOne.header.getSequenceDictionary)
      detector.getAll should have size 1
      val intervals = detector.getAll.toSeq.sortBy(_.getStart)
      val head = intervals.head
      head.getStart shouldBe 1
      head.getEnd shouldBe 3
    }
  }

  it should "truncate the smaller block when too blocks overlap" in {
    val builderBlockOne = new VariantContextSetBuilder()
    builderBlockOne.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=3, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=4, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)

    // first block is extended
    {
      val builderBlockTwo = new VariantContextSetBuilder()
      builderBlockTwo.addVariant(start = 3, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
      builderBlockTwo.addVariant(start = 4, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
      builderBlockTwo.addVariant(start = 5, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)

      val contexts = (builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 2) } ++ builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 3) }).toSeq

      // Check that if we do not want to modify the blocks we get an exception
      an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

      val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
      detector.getAll should have size 2
      val intervals = detector.getAll.toSeq.sortBy(_.getStart)
      val head = intervals.head
      head.getStart shouldBe 2
      head.getEnd shouldBe 4
      val last = intervals.last
      last.getStart shouldBe 5
      last.getEnd shouldBe 5
    }

    // second block is extended
    {
      val builderBlockTwo = new VariantContextSetBuilder()
      builderBlockTwo.addVariant(start = 3, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
      builderBlockTwo.addVariant(start = 4, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
      builderBlockTwo.addVariant(start = 5, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
      builderBlockTwo.addVariant(start = 6, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)

      val contexts = (builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 2) } ++ builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 3) }).toSeq

      // Check that if we do not want to modify the blocks we get an exception
      an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

      val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
      detector.getAll should have size 2
      val intervals = detector.getAll.toSeq.sortBy(_.getStart)
      val head = intervals.head
      head.getStart shouldBe 2
      head.getEnd shouldBe 2
      val last = intervals.last
      last.getStart shouldBe 3
      last.getEnd shouldBe 6
    }
  }

  it should "resolve three overlapping blocks, such that when the middle one is truncated and now starts after the third, it is resolved" in {
    val builderBlockOne = new VariantContextSetBuilder()
    builderBlockOne.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=10, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val oneIter = builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 2) }

    val builderBlockTwo = new VariantContextSetBuilder()
    builderBlockTwo.addVariant(start = 8, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builderBlockTwo.addVariant(start = 14, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    val twoIter = builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 3) }

    val builderBlockThree = new VariantContextSetBuilder()
    builderBlockThree.addVariant(start = 9, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builderBlockThree.addVariant(start = 13, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    val threeIter = builderBlockThree.iterator.map { ctx => withPhasingSetId(ctx, 4) }

    val contexts = (oneIter ++ twoIter ++ threeIter).toSeq

    // Check that if we do not want to modify the blocks we get an exception
    an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

    val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
    detector.getAll should have size 3
    val intervals = detector.getAll.toSeq.sortBy(_.getStart)

    intervals.length shouldBe 3

    val firstInterval = intervals.head
    firstInterval.getStart shouldBe 2
    firstInterval.getEnd shouldBe 10

    val secondInterval = intervals(1)
    secondInterval.getStart shouldBe 9
    secondInterval.getEnd shouldBe 13

    val thirdInterval = intervals(2)
    thirdInterval.getStart shouldBe 14
    thirdInterval.getEnd shouldBe 14
  }

  it should "resolve three overlapping blocks, such that when the middle one is truncated and now is enclosed in the third, we get two blocks" in {
    val builderBlockOne = new VariantContextSetBuilder()
    builderBlockOne.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    builderBlockOne.addVariant(start=10, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val oneIter = builderBlockOne.iterator.map { ctx => withPhasingSetId(ctx, 2) }

    val builderBlockTwo = new VariantContextSetBuilder()
    builderBlockTwo.addVariant(start = 8, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builderBlockTwo.addVariant(start = 13, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    val twoIter = builderBlockTwo.iterator.map { ctx => withPhasingSetId(ctx, 3) }

    val builderBlockThree = new VariantContextSetBuilder()
    builderBlockThree.addVariant(start = 9, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builderBlockThree.addVariant(start = 13, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    val threeIter = builderBlockThree.iterator.map { ctx => withPhasingSetId(ctx, 4) }

    val contexts = (oneIter ++ twoIter ++ threeIter).toSeq

    // Check that if we do not want to modify the blocks we get an exception
    an[Exception] should be thrownBy PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary, modifyBlocks = false)

    val detector = PhaseBlock.buildOverlapDetector(iterator = contexts.iterator, dict = builderBlockOne.header.getSequenceDictionary)
    detector.getAll should have size 2
    val intervals = detector.getAll.toSeq.sortBy(_.getStart)

    intervals.length shouldBe 2

    val firstInterval = intervals.head
    firstInterval.getStart shouldBe 2
    firstInterval.getEnd shouldBe 10

    val secondInterval = intervals(1)
    secondInterval.getStart shouldBe 9
    secondInterval.getEnd shouldBe 13
  }
}

class PhaseCigarTest extends ErrorLogLevel {
  import AssessPhasingTest._
  import PhaseCigarOp._

  private def toCigar(truth: Seq[VariantContext],
                      call: Seq[VariantContext],
                      header: VCFHeader,
                      skipMismatchingAlleles: Boolean,
                      assumeFixedAlleleOrder: Boolean = true): (Seq[PhaseCigarOp], AssessPhasingMetric) = {
    val dict = header.getSequenceDictionary
    val pairedIterator = JointVariantContextIterator(
      iters       = Seq(truth.iterator, call.iterator),
      dict        = dict
    ).map { case Seq(left, right) => (left, right) }
    val truthPhaseBlockDetector  = PhaseBlock.buildOverlapDetector(truth.iterator, dict)
    val calledPhaseBlockDetector = PhaseBlock.buildOverlapDetector(call.iterator, dict)
    val metric                   = AssessPhasingMetric()

    val cigar = PhaseCigar(
      pairedIterator           = pairedIterator,
      truthPhaseBlockDetector  = truthPhaseBlockDetector,
      calledPhaseBlockDetector = calledPhaseBlockDetector,
      metric                   = metric,
      skipMismatchingAlleles   = skipMismatchingAlleles,
      assumeFixedAlleleOrder   = assumeFixedAlleleOrder
     )

    metric.finalizeMetric()

    (cigar.cigar, metric)
  }

  "Cigar.toCigar" should "create an empty cigar if no variants have a phasing set" in {
    val builder = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val ctx     = builder.head

    // truth variant only
    {
      val (cigar, metric) = toCigar(truth = Seq(ctx), call = Seq.empty, header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, BothEnd)
      metric.num_called shouldBe 0
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 0
    }

    // call variant only
    {
      val (cigar, metric) = toCigar(truth = Seq.empty, call = Seq(ctx), header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 0
      metric.num_truth shouldBe 0
    }
  }

  it should "create a cigar from either a single truth or call variant" in {
    val builder = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val ctx     = withPhasingSetId(builder.head, 1)

    // truth variant only
    {
      val (cigar, metric) = toCigar(truth = Seq(ctx), call = Seq.empty, header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, TruthOnly, BothEnd)
      metric.num_called shouldBe 0
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }

    // call variant only
    {
      val (cigar, metric) = toCigar(truth = Seq.empty, call = Seq(ctx), header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, CallOnly, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 0
    }
  }

  it should "create a cigar when both truth and call variants are present but only of the two are phased" in {
    val builder    = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val ctxPhased  = withPhasingSetId(builder.head, 1)
    val ctxNoPhase = builder.head

    // truth variant phased only
    {
      val (cigar, metric) = toCigar(truth = Seq(ctxPhased), call = Seq(ctxNoPhase), header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, TruthOnly, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 0
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }

    // call variant phased only
    {
      val (cigar, metric) = toCigar(truth = Seq(ctxNoPhase), call = Seq(ctxPhased), header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, CallOnly, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 0
    }
  }

  it should "create a cigar when both truth and call variants are present and both are phased" in {
    val builder  = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val ctx      = withPhasingSetId(builder.head, 1)

    // both variants are phased
    {
      val (cigar, metric) = toCigar(truth = Seq(ctx), call = Seq(ctx), header = builder.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, Match, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }
  }

  it should "create a cigar when both truth and call variants are present and both are phased but mismatch alleles" in {
    val builderTruth = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val builderCall  = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("C", "A"), phased=true)
    val truth = withPhasingSetId(builderTruth.head, 1)
    val call  = withPhasingSetId(builderCall.head, 1)

    // both variants are phased, phase is inverted, and we assume a fixed order, so a mismatch
    {
      val (cigar, metric) = toCigar(truth = Seq(truth), call = Seq(call), header = builderTruth.header, skipMismatchingAlleles = true, assumeFixedAlleleOrder = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, Mismatch, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }
    // both variants are phased, phase is inverted, and we don't assume a fixed order, so a match
    {
      val (cigar, metric) = toCigar(truth = Seq(truth), call = Seq(call), header = builderTruth.header, skipMismatchingAlleles = true, assumeFixedAlleleOrder = false)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, Match, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }
    // first site is a match, second is a mismatch since we inverted after the first
    {
      builderTruth.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
      builderCall.addVariant(start=2, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
      val truthTwo = withPhasingSetId(builderTruth.last, 1)
      val callTwo  = withPhasingSetId(builderCall.last, 1)

      val (cigar, metric) = toCigar(truth = Seq(truth, truthTwo), call = Seq(call, callTwo), header = builderTruth.header, skipMismatchingAlleles = true, assumeFixedAlleleOrder = false)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, Match, Mismatch, BothEnd)
      metric.num_called shouldBe 2
      metric.num_phased shouldBe 2
      metric.num_truth shouldBe 2
      metric.num_truth_phased shouldBe 2
    }
  }

  it should "skip sites where alleles mismatch if specified" in {
    val builderTruth = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "C"), genotypeAlleles=List("A", "C"), phased=true)
    val builderCall  = new VariantContextSetBuilder().addVariant(start=1, variantAlleles=List("A", "G"), genotypeAlleles=List("A", "G"), phased=true)
    val truth = withPhasingSetId(builderTruth.head, 1)
    val call  = withPhasingSetId(builderCall.head, 1)

    // skip the sites that have mismatching alleles
    {
      val (cigar, metric) = toCigar(truth = Seq(truth), call = Seq(call), header = builderTruth.header, skipMismatchingAlleles = true)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, BothEnd)
      metric.num_called shouldBe 0
      metric.num_truth shouldBe 0
    }

    // include the sites that have mismatching alleles
    {
      val (cigar, metric) = toCigar(truth = Seq(truth), call = Seq(call), header = builderTruth.header, skipMismatchingAlleles = false)
      cigar should contain theSameElementsInOrderAs Seq(BothEnd, Mismatch, BothEnd)
      metric.num_called shouldBe 1
      metric.num_phased shouldBe 1
      metric.num_truth shouldBe 1
      metric.num_truth_phased shouldBe 1
    }
  }

  private val threeHaplotypeCigar = Seq(BothEnd, Match, TruthOnly, CallOnly, Match, BothEnd, Match, Mismatch, CallOnly, TruthOnly, Match, TruthOnly, CallEnd, CallOnly, BothEnd)
  it should "create a cigar with three haplotype blocks" in {
    val truth = AssessPhasingTest.TruthVariants.filter(ctx => ctx.getStart < 30)
    val call  = AssessPhasingTest.CallVariants.filter(ctx => ctx.getStart < 30)

    val (cigar, metric) = toCigar(truth = truth, call = call, header = AssessPhasingTest.Header, skipMismatchingAlleles = false)
    cigar should contain theSameElementsInOrderAs threeHaplotypeCigar
    metric.num_called shouldBe 9
    metric.num_phased shouldBe 8
    metric.num_truth shouldBe 9
    metric.num_truth_phased shouldBe 8
    metric.frac_phased shouldBe 8/9d
    metric.num_called_with_truth_phased shouldBe 6
    metric.num_phased_with_truth_phased shouldBe 5
    metric.frac_phased_with_truth_phased shouldBe 5/6d
    metric.num_truth_phased_in_called_block shouldBe 7
    metric.num_both_phased_in_called_block shouldBe 5
    metric.frac_truth_phased_in_called_block shouldBe 7/8d
    metric.frac_phased_with_truth_phased_in_called_block shouldBe 5/7d
  }

  it should "create a cigar with four haplotype blocks" in {
    val fourHaplotypeBlockCigar = threeHaplotypeCigar ++ Seq(Match, Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch, BothEnd)
    val truth = AssessPhasingTest.TruthVariants.filter(ctx => ctx.getStart <= 42)
    val call  = AssessPhasingTest.CallVariants.filter(ctx => ctx.getStart <= 42)

    val (cigar, metric) = toCigar(truth = truth, call = call, header = AssessPhasingTest.Header, skipMismatchingAlleles = false)
    cigar should contain theSameElementsInOrderAs fourHaplotypeBlockCigar
    metric.num_called shouldBe 22
    metric.num_phased shouldBe 21
    metric.num_truth shouldBe 22
    metric.num_truth_phased shouldBe 21
    metric.frac_phased shouldBe 21/22d
    metric.num_called_with_truth_phased shouldBe 19
    metric.num_phased_with_truth_phased shouldBe 18
    metric.frac_phased_with_truth_phased shouldBe 18/19d
    metric.num_truth_phased_in_called_block shouldBe 20
    metric.num_both_phased_in_called_block shouldBe 18
    metric.frac_truth_phased_in_called_block shouldBe 20/21d
    metric.frac_phased_with_truth_phased_in_called_block shouldBe 18/20d
  }

  "Cigar.toShortSwitchErrorIndices" should "find no indices for an empty cigar" in {
    PhaseCigar(cigar=Seq.empty).toShortSwitchErrorIndices()shouldBe 'empty
  }

  it should "find an error at the first or last cigar element if the mismatch is isolated" in {
    PhaseCigar(cigar=Seq(Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(0)
    PhaseCigar(cigar=Seq(Mismatch, Match)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(0)
    PhaseCigar(cigar=Seq(Match, Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1)
    PhaseCigar(cigar=Seq(Match, Match, Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(2)
    PhaseCigar(cigar=Seq(Mismatch, Match, Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(0, 2)
  }

  it should "find an error in the middle of a haplotype block" in {
    PhaseCigar(cigar=Seq(Match, Mismatch, Match)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1)
    Seq(CallEnd, TruthEnd, BothEnd).foreach { endMarker =>
      PhaseCigar(cigar=Seq(Match, Mismatch, Match, endMarker, Match)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1)
      PhaseCigar(cigar=Seq(Match, Mismatch, Match, endMarker, Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1, 4)
      PhaseCigar(cigar=Seq(Match, Mismatch, Match, endMarker, Match, Mismatch)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1, 5)
      PhaseCigar(cigar=Seq(Match, Mismatch, Match, endMarker, Match, Mismatch, Match)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1, 5)
      PhaseCigar(cigar=Seq(Match, Mismatch, Match, endMarker, Match, endMarker, Mismatch, Match)).toShortSwitchErrorIndices()should contain theSameElementsInOrderAs Seq(1, 6)
    }
  }

  it should "ignore variant sites only in the truth or call set" in {
    PhaseCigar(cigar=Seq(CallOnly, Match)).toShortSwitchErrorIndices()shouldBe 'empty
    PhaseCigar(cigar=Seq(TruthOnly, Match)).toShortSwitchErrorIndices()shouldBe 'empty
  }

  it should "not count adjacent mismatches" in {
    PhaseCigar(cigar=Seq(Mismatch, Mismatch)).toShortSwitchErrorIndices()shouldBe 'empty
  }

  "Cigar.toLongSwitchErrorsAndSites" should "find no errors and no sites for an empty cigar" in {
    PhaseCigar(cigar=Seq.empty).toLongSwitchErrorsAndSites() shouldBe (0, 0)
  }

  it should "return no errors if Mismatch is not part of the cigar" in {
    PhaseCigar(cigar=Seq(Match, Match, Match, Match, Match)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Match, CallOnly, Match, TruthOnly, Match, Match, Match)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Match, CallOnly, Match, TruthOnly, CallEnd, Match, Match, Match, CallEnd)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Match, CallOnly, Match, TruthOnly, TruthEnd, Match, Match, Match, TruthEnd)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
  }

  it should "return no errors if there are only isolated mismatches" in {
    PhaseCigar(cigar=Seq(Match,    Match, Mismatch, Match, Match)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Mismatch, Match, Match,    Match, Match)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Match,    Match, Match,    Match, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Mismatch, Match, Match,    Match, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
    PhaseCigar(cigar=Seq(Mismatch, Match, Mismatch, Match, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (0, 5)
  }

  it should "ignore variant sites only in the truth or call set" in {
    PhaseCigar(cigar = Seq(CallOnly, Match, Match, Match, Match)).toLongSwitchErrorsAndSites() shouldBe(0, 4)
    PhaseCigar(cigar = Seq(TruthOnly, Match, Match, Match, Match)).toLongSwitchErrorsAndSites() shouldBe(0, 4)
  }

  it should "return a single error for an adjacent mismatches" in {
    PhaseCigar(cigar=Seq(Mismatch, Mismatch, Match,    Match,    Match)).toLongSwitchErrorsAndSites()    shouldBe (1, 5)
    PhaseCigar(cigar=Seq(Match,    Mismatch, Mismatch, Match,    Match)).toLongSwitchErrorsAndSites()    shouldBe (1, 5)
    PhaseCigar(cigar=Seq(Match,    Match,    Mismatch, Mismatch, Match)).toLongSwitchErrorsAndSites()    shouldBe (1, 5)
    PhaseCigar(cigar=Seq(Match,    Match,    Match,    Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (1, 5)

    // introduce a single mismatch
    PhaseCigar(cigar=Seq(Mismatch, Match, Match, Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (1, 5)

    // all mismatches
    PhaseCigar(cigar=Seq(Mismatch, Mismatch, Mismatch, Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe (1, 5)
  }

  it should "find no errors if there are no adjacent mismatches within a block" in {
    Seq(CallEnd, TruthEnd, BothEnd).foreach { endMarker =>
      PhaseCigar(cigar = Seq(Mismatch, endMarker, Mismatch, Match, Match, Match)).toLongSwitchErrorsAndSites() shouldBe(0, 5)
      PhaseCigar(cigar = Seq(Match, Mismatch, endMarker, Mismatch, Match, Match)).toLongSwitchErrorsAndSites() shouldBe(0, 5)
      PhaseCigar(cigar = Seq(Match, Match, Mismatch, endMarker, Mismatch, Match)).toLongSwitchErrorsAndSites() shouldBe(0, 5)
      PhaseCigar(cigar = Seq(Match, Match, Match, Mismatch, endMarker, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(0, 5)

      // introduce a single mismatch
      PhaseCigar(cigar = Seq(Mismatch, Match, Match, Mismatch, endMarker, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(0, 5)
    }
  }

  it should "find a single error if there is an adjacent mismatch within one block" in {
    Seq(CallEnd, TruthEnd, BothEnd).foreach { endMarker =>
      PhaseCigar(cigar = Seq(Mismatch, Mismatch, Mismatch, Mismatch, endMarker, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(1, 5)
      PhaseCigar(cigar = Seq(Mismatch, Mismatch, Mismatch, Mismatch, endMarker, Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(2, 6)
      PhaseCigar(cigar = Seq(Match, Match, Match, Match, endMarker, Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(1, 6)
      PhaseCigar(cigar = Seq(Match, Match, Match, Mismatch, endMarker, Mismatch, Mismatch)).toLongSwitchErrorsAndSites() shouldBe(1, 6)
    }
  }

  "Cigar.toIlluminaSwitchErrorsHmm" should "return no errors for a single variant site" in {
    PhaseCigar(Seq(Match)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 1)
    PhaseCigar(Seq(Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 1)
  }

  it should "return no errors for a phase block with one phase" in {
    PhaseCigar(Seq(Match, Match, Match, Match, Match)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 5)
    PhaseCigar(Seq(Mismatch, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 5)
  }

  it should "return a point error for a phase block with one phase but a point error" in {
    // point error at the start
    PhaseCigar(Seq(Mismatch, Match, Match, Match, Match)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 5)
    PhaseCigar(Seq(Match, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 5)
    // point error in the middle
    PhaseCigar(Seq(Match, Match, Mismatch, Match, Match)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 5)
    PhaseCigar(Seq(Mismatch, Mismatch, Match, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 5)
    // ponit error at the end
    PhaseCigar(Seq(Match, Match, Match, Match, Match, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 6)
    PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 0, 7)
  }

  it should "return a stretch of point errors up to (transitionPenalty / pointPenalty) point errors" in {
    // transitionPenalty is -5
    {
      // point error, since four point errors have a score of -1, while one long switch error would have a score of -5
      PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(4, 0, 10)
      // long switch error, since five point errors have a score of -5, while one long switch error would have a score of -5
      PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 1, 11)
    }

    // transitionPenalty is -4
    {
      // point error, since four point errors have a score of -1, while one long switch error would have a score of -4
      PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors(transitionPenalty = -4) shouldBe IlluminaSwitchErrors(3, 0, 10)
      // long switch error, since five point errors have a score of -4, while one long switch error would have a score of -4
      PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors(transitionPenalty = -4) shouldBe IlluminaSwitchErrors(0, 1, 11)
    }
  }

  it should "return a long switch error for a phase block that switches haplotypes once" in {
    // long switch error, since six point errors have a score of -6, while one long switch error would have a score of -5
    PhaseCigar(Seq(Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 1, 12)
  }

  it should "return a point and a long switch error for a phase block that switches haplotypes once" in {
    // a point error to start, and a long switch error stretching to the end
    PhaseCigar(Seq(Mismatch, Match, Match, Match, Match, Match, Match, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch, Mismatch)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(1, 1, 13)
  }

  it should "return no errors and no sites with no Match/Mismatch CigarTypes" in {
    PhaseCigar(Seq(CallOnly)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 0)
    PhaseCigar(Seq(CallEnd)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 0)
    PhaseCigar(Seq(TruthOnly)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 0)
    PhaseCigar(Seq(TruthEnd)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 0)
    PhaseCigar(Seq(CallOnly, CallEnd, TruthOnly, TruthEnd)).toIlluminaSwitchErrors() shouldBe IlluminaSwitchErrors(0, 0, 0)
  }

  it should "throw an Exception if the cigar was empty" in {
    an[Exception] should be thrownBy PhaseCigar(Seq.empty).toIlluminaSwitchErrorsHmm()
  }

  "Cigar.toPhasedBlocks" should "partition a cigar into phased blocks" in {

    // single blocks
    PhaseCigar(Seq(CallOnly, CallEnd)).toPhasedBlocks().map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(CallOnly))
    PhaseCigar(Seq(CallOnly, BothEnd)).toPhasedBlocks().map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(CallOnly))
    PhaseCigar(Seq(TruthOnly, TruthEnd)).toPhasedBlocks(isTruth=true).map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(TruthOnly))
    PhaseCigar(Seq(TruthOnly, BothEnd)).toPhasedBlocks(isTruth=true).map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(TruthOnly))

    // two blocks with a single call
    PhaseCigar(Seq(CallOnly, CallEnd, CallOnly, CallEnd)).toPhasedBlocks().map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(CallOnly), Seq(CallOnly))
    PhaseCigar(Seq(CallOnly, CallEnd, TruthOnly, CallEnd)).toPhasedBlocks().map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(CallOnly), Seq(TruthOnly))

    // two blocks but with a truth end for fun
    PhaseCigar(Seq(CallOnly, CallEnd, TruthOnly, TruthEnd, CallOnly, CallEnd)).toPhasedBlocks().map(_.cigar) should contain theSameElementsInOrderAs Seq(Seq(CallOnly), Seq(TruthOnly, TruthEnd, CallOnly))
  }

  "Cigar.contextsToBlockEndOperator/contextsToMatchingOperator" should "should return the cigar operator for two variant contexts" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start = 1, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    val ctxStart   = withPhasingSetId(builder.head, 1)
    val ctxNoStart = withPhasingSetId(builder.last, 1)

    // No truth, call is start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=None, call=Some(ctxStart)) shouldBe Some(CallEnd)
    PhaseCigar.contextsToMatchingOperator(truth=None, call=Some(ctxStart)) shouldBe Some(CallOnly)

    // No truth, call is not the start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=None, call=Some(ctxNoStart)) shouldBe 'empty
    PhaseCigar.contextsToMatchingOperator(truth=None, call=Some(ctxNoStart)) shouldBe Some(CallOnly)

    // Truth is the start of a phase block, no call
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxStart), call=None) shouldBe Some(TruthEnd)
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxStart), call=None) shouldBe Some(TruthOnly)

    // Truth is the start of a phase block, no call
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxNoStart), call=None) shouldBe 'empty
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxNoStart), call=None) shouldBe Some(TruthOnly)

    // Truth and call, both start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxStart), call=Some(ctxStart)) shouldBe Some(BothEnd)
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxStart), call=Some(ctxStart)) shouldBe Some(Match)

    // Truth and call, only call is the start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxNoStart), call=Some(ctxStart)) shouldBe Some(CallEnd)
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxNoStart), call=Some(ctxStart)) shouldBe Some(Match)

    // Truth and call, only truth is the start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxStart), call=Some(ctxNoStart)) shouldBe Some(TruthEnd)
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxStart), call=Some(ctxNoStart)) shouldBe Some(Match)

    // Truth and call, neither are the start of a phase block
    PhaseCigar.contextsToBlockEndOperator(truth=Some(ctxNoStart), call=Some(ctxNoStart)) shouldBe 'empty
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctxNoStart), call=Some(ctxNoStart)) shouldBe Some(Match)
  }

  it should "should return the cigar operator for two variant contexts that disagree on phase" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start = 1, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"), phased = true)
    builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("C", "A"), phased = true)
    val ctx         = withPhasingSetId(builder.head, 1)
    val ctxMismatch = withPhasingSetId(builder.last, 2)

    // Truth and call, both start of a phase block
    PhaseCigar.contextsToMatchingOperator(truth=Some(ctx), call=Some(ctxMismatch)) shouldBe Some(Mismatch)
  }

  it should "be unreachable to provide (None, None)" in {
    an[Exception] should be thrownBy PhaseCigar.contextsToBlockEndOperator(truth=None, call=None)
    an[Exception] should be thrownBy PhaseCigar.contextsToMatchingOperator(truth=None, call=None)
  }

  "Cigar.cigarForVariantContexts" should "return a match if two variant contexts share the same alleles in the same order" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start = 1, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"))
    builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"))
    PhaseCigar.cigarTypeForVariantContexts(builder.head, builder.last) shouldBe Match
    PhaseCigar.cigarTypeForVariantContexts(builder.last, builder.head) shouldBe Match
  }

  it should "return a mismatch if two variant contexts share the same alleles but in the different order" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start = 1, variantAlleles = List("A", "C"), genotypeAlleles = List("C", "A"))
    builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"))
    PhaseCigar.cigarTypeForVariantContexts(builder.head, builder.last) shouldBe Mismatch
    PhaseCigar.cigarTypeForVariantContexts(builder.last, builder.head) shouldBe Mismatch
  }

  it should "return a mismatch if two variant contexts share the different alleles" in {
    val builder = new VariantContextSetBuilder()
    builder.addVariant(start = 1, variantAlleles = List("A", "G"), genotypeAlleles = List("A", "G"))
    builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"))
    PhaseCigar.cigarTypeForVariantContexts(builder.head, builder.last) shouldBe Mismatch
    PhaseCigar.cigarTypeForVariantContexts(builder.last, builder.head) shouldBe Mismatch
  }

  it should "throw an exception if the variant contexts do not have exactly two alleles" in {
    // One variant alleles, one genotype allele
    {
      val builder = new VariantContextSetBuilder()
      builder.addVariant(start = 1, variantAlleles = List("A"), genotypeAlleles = List("A"))
      builder.addVariant(start = 2, variantAlleles = List("A"), genotypeAlleles = List("A"))
      an[Exception] should be thrownBy PhaseCigar.cigarTypeForVariantContexts(builder.head, builder.last)
    }

    // Two variant alleles, one genotype allele
    {
      val builder = new VariantContextSetBuilder()
      builder.addVariant(start = 1, variantAlleles = List("A", "C"), genotypeAlleles = List("C"))
      builder.addVariant(start = 2, variantAlleles = List("A", "C"), genotypeAlleles = List("A", "C"))
      an[Exception] should be thrownBy PhaseCigar.cigarTypeForVariantContexts(builder.head, builder.last)
      an[Exception] should be thrownBy PhaseCigar.cigarTypeForVariantContexts(builder.last, builder.head)
    }
  }
}
