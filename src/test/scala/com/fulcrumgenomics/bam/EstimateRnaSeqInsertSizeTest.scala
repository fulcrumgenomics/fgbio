/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.bam
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.EstimateRnaSeqInsertSize._
import com.fulcrumgenomics.testing.SamRecordSetBuilder._
import com.fulcrumgenomics.testing.{SamRecordSetBuilder, UnitSpec}
import com.fulcrumgenomics.util.GeneAnnotations.{Exon, Gene, Transcript}
import com.fulcrumgenomics.util.{Io, Metric}
import dagr.commons.io.PathUtil
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamPairUtil.PairOrientation
import org.scalatest.OptionValues

class EstimateRnaSeqInsertSizeTest extends UnitSpec with OptionValues {
  /** Calculates the insert size from a gene.  Returns None if the record's span is not enclosed in the gene or if
    * the insert size disagree across transcripts. */
  def testInsertSizeFromGene(rec: SAMRecord,
                             gene: Gene,
                             minimumOverlap: Double): Option[Int] = {
    val mateCigar        = EstimateRnaSeqInsertSize.getAndRequireMateCigar(rec)
    val mateAlignmentEnd = EstimateRnaSeqInsertSize.mateAlignmentEndFrom(mateCigar, rec.getMateAlignmentStart)
    EstimateRnaSeqInsertSize.insertSizeFromGene(
      rec              = rec,
      gene             = gene,
      minimumOverlap   = minimumOverlap,
      recInterval      = EstimateRnaSeqInsertSize.intervalFrom(rec=rec, mateAlignmentEnd=mateAlignmentEnd),
      recBlocks        = rec.getAlignmentBlocks.toList,
      mateBlocks       = EstimateRnaSeqInsertSize.mateAlignmentBlocksFrom(mateCigar, rec.getMateAlignmentStart),
      mateAlignmentEnd = mateAlignmentEnd
    )
  }

  def mateAlignmentEnd(rec: SAMRecord): Int = {
    val mateCigar        = EstimateRnaSeqInsertSize.getAndRequireMateCigar(rec)
    EstimateRnaSeqInsertSize.mateAlignmentEndFrom(mateCigar, rec.getMateAlignmentStart)
  }


  "EstimateRnaSeqInsertSize.numReadBasesOverlappingTranscript" should "return the number of read bases overlapping a transcript" in {
    def estimate(rec: SAMRecord, transcript: Transcript) =  EstimateRnaSeqInsertSize.numReadBasesOverlappingTranscript(rec.getAlignmentBlocks.toList, transcript)

    // many .value calls here, I know
    val transcript = Transcript("", 2, 10, 2, 10, exons=Seq(Exon(2,2), Exon(4,4), Exon(11, 11)))
    val builder    = new SamRecordSetBuilder(readLength=10)

    // simple matches
    builder.addFrag(start=1, strand=Plus)   foreach { rec => estimate(rec, transcript) shouldBe 2 }
    builder.addFrag(start=1, strand=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 2 }
    builder.addFrag(start=2, strand=Plus)   foreach { rec => estimate(rec, transcript) shouldBe 3 }
    builder.addFrag(start=2, strand=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 3 }
    builder.addFrag(start=5, strand=Plus)   foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addFrag(start=5, strand=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addFrag(start=11, strand=Plus)  foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addFrag(start=11, strand=Minus) foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addFrag(start=12, strand=Plus)  foreach { rec => estimate(rec, transcript) shouldBe 0 }
    builder.addFrag(start=12, strand=Minus) foreach { rec => estimate(rec, transcript) shouldBe 0 }

    // some indels
    builder.addFrag(start=2, cigar="1M1D1M6D8M")   foreach { rec => estimate(rec, transcript) shouldBe 3 } // deletions between exons
    builder.addFrag(start=2, cigar="1M9I1M")       foreach { rec => estimate(rec, transcript) shouldBe 1 } // insertions
    builder.addFrag(start=1, cigar="1M20D9M")      foreach { rec => estimate(rec, transcript) shouldBe 0 } // deletion skips gene

    // skips
    builder.addFrag(start=2, cigar="1M1N1M6N8M")   foreach { rec => estimate(rec, transcript) shouldBe 3 } // skips between exons
  }

  "EstimateRnaSeqInsertSize.insertSizeFrom" should "return the number of read bases overlapping a transcript" in {
    def estimate(rec: SAMRecord, transcript: Transcript) = insertSizeFromTranscript(rec, transcript, mateAlignmentEnd(rec))

    // many .value calls here, I know
    val transcript = Transcript("", 2, 10, 2, 10, exons=Seq(Exon(2,2), Exon(4,4), Exon(11, 11)))
    val builder    = new SamRecordSetBuilder(readLength=5)

    // Overlaps all three exons
    builder.addPair(start1=1, start2=7, strand1=Plus, strand2=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 3 }
    builder.addPair(start1=7, start2=1, strand1=Minus, strand2=Plus)  foreach { rec => estimate(rec, transcript) shouldBe 3 }

    // Overlaps the first two exons
    builder.addPair(start1=1, start2=7, strand1=Plus, strand2=Plus)   foreach { rec => estimate(rec, transcript) shouldBe 2 }

    // Overlaps the last exon
    builder.addPair(start1=7, start2=1, strand1=Minus, strand2=Minus) foreach { rec => estimate(rec, transcript) shouldBe 1 }

    // Overlaps all last exon
    builder.addPair(start1=7, start2=12, strand1=Plus, strand2=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addPair(start1=12, start2=7, strand1=Minus, strand2=Plus)  foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addPair(start1=7, start2=12, strand1=Plus, strand2=Plus)   foreach { rec => estimate(rec, transcript) shouldBe 1 }
    builder.addPair(start1=12, start2=7, strand1=Minus, strand2=Minus) foreach { rec => estimate(rec, transcript) shouldBe 1 }

    // No overlap (5' positions are 26 and 12)
    builder.addPair(start1=22, start2=8, strand1=Minus, strand2=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 0 }

    // One base overlap (5' positions are 26 and 11)
    builder.addPair(start1=22, start2=7, strand1=Minus, strand2=Minus)  foreach { rec => estimate(rec, transcript) shouldBe 1 }
  }

  it should "return a value if the record overlaps a gene" in {
    val transcript = Transcript("example_transcript", 10, 20, 10, 20, exons=Seq(Exon(10,10), Exon(14,14), Exon(20,20)))
    val gene       = Gene(contig="chr1", start=10, end=20, negativeStrand=false, name="", transcripts=Seq(transcript))
    val builder    = new SamRecordSetBuilder(readLength=5)

    ///////////////////////////////////////////////////////
    // not enclosed
    ///////////////////////////////////////////////////////

    // too far left
    builder.addPair(start1=9, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=16, start2=9, strand1=Minus, strand2=Plus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=9, start2=16, strand1=Plus, strand2=Plus)   foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=16, start2=9, strand1=Minus, strand2=Minus) foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }

    // too far right
    builder.addPair(start1=10, start2=17, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=17, start2=10, strand1=Minus, strand2=Plus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=10, start2=17, strand1=Plus, strand2=Plus)   foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    builder.addPair(start1=17, start2=10, strand1=Minus, strand2=Minus) foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
    
    ///////////////////////////////////////////////////////
    // enclosed
    ///////////////////////////////////////////////////////

    // enclosed (just barely)
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0).value shouldBe 3 }
    builder.addPair(start1=16, start2=10, strand1=Minus, strand2=Plus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0).value shouldBe 3 }
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Plus)   foreach { rec => testInsertSizeFromGene(rec, gene, 0.0).value shouldBe 2 }
    builder.addPair(start1=16, start2=10, strand1=Minus, strand2=Minus) foreach { rec => testInsertSizeFromGene(rec, gene, 0.0).value shouldBe 2 }
  }

  it should "not return a value if there is too little overlap" in {
    val transcript = Transcript("example_transcript", 10, 20, 10, 20, exons=Seq(Exon(10,10), Exon(20, 20)))
    val gene       = Gene(contig="chr1", start=10, end=20, negativeStrand=false, name="", transcripts=Seq(transcript))
    val builder    = new SamRecordSetBuilder(readLength=5)

    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0).value shouldBe 2 }
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.2).value shouldBe 2 }
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.2001) shouldBe 'empty }
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.3) shouldBe 'empty }
    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 1.0) shouldBe 'empty }
  }

  it should "not return a value if the insert size disagrees across two transcripts" in {
    val transcriptA = Transcript("example_transcript_A", 10, 20, 10, 20, exons=Seq(Exon(10,10), Exon(20, 20)))
    val transcriptB = Transcript("example_transcript_B", 10, 20, 10, 20, exons=Seq(Exon(10,10), Exon(19, 20))) // one longer than A
    val gene        = Gene(contig="chr1", start=10, end=20, negativeStrand=false, name="", transcripts=Seq(transcriptA, transcriptB))
    val builder     = new SamRecordSetBuilder(readLength=5)

    builder.addPair(start1=10, start2=16, strand1=Plus, strand2=Minus)  foreach { rec => testInsertSizeFromGene(rec, gene, 0.0) shouldBe 'empty }
  }

  private val RefFlatFile = {
    val lines = Seq(
      // a run-of-the-mill gene, with one transcript and one exon
      Seq("ACKR4-3", "NM_178445-A", "chr3", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175"),
      // two genes that overlap
      Seq("ACKR4-4-1", "NM_178445-B", "chr4", "+", "133801670", "133804175", "133801931", "133802984", "1", "133801670", "133804175"),
      Seq("ACKR4-4-2", "NM_178445-C", "chr4", "+", "133801671", "133804176", "133801931", "133802984", "1", "133801671", "133804176"),
      // two transcripts that overlap but have different lengths
      Seq("ACKR4-5", "NM_178445-D", "chr5", "+", "133801670", "133804176", "133801931", "133802985", "1", "133801670", "133804175"),
      Seq("ACKR4-5", "NM_178445-E", "chr5", "+", "133801670", "133804176", "133801931", "133802985", "1", "133801670", "133804176"),
      // a transcript with two exons
      Seq("ACKR4-6", "NM_178445-F", "chr6", "+", "133801670", "133804176", "133801931", "133802985", "2", "133801670,133804175", "133801671,133804176")

    ).map(_.mkString("\t"))
    val refFlat = makeTempFile("refFlat.", ".txt")
    Io.writeLines(path=refFlat, lines=lines)
    refFlat
  }

  private val EmptyMetrics = PairOrientation.values().map { po => InsertSizeMetric(po) }

  "EstimateRnaSeqInsertSize" should "run end-to-end" in {
    val builder = new SamRecordSetBuilder()
    // FR = (133804075 + 100) - 133801671 = 2504
    builder.addPair(contig=2, start1=133801671, start2=133804075, strand1=Plus, strand2=Minus) // overlaps ACKR4 by 100%
    builder.addPair(contig=2, start1=133801672, start2=133804074, strand1=Plus, strand2=Minus) //insert is two less

    // RF = (133804075 + 1) - (133801671 + 100 - 1) = 2306
    builder.addPair(contig=2, start1=133801671, start2=133804075, strand1=Minus, strand2=Plus) // overlaps ACKR4 by 100%
    builder.addPair(contig=2, start1=133801672, start2=133804074, strand1=Minus, strand2=Plus) // overlaps ACKR4 by 100%

    // TANDEM = (133804075 + 1) - 133801671 = 2405
    builder.addPair(contig=2, start1=133801671, start2=133804075, strand1=Plus, strand2=Plus) // overlaps ACKR4 by 100%
    builder.addPair(contig=2, start1=133801672, start2=133804074, strand1=Plus, strand2=Plus) // overlaps ACKR4 by 100%

    val bam = builder.toTempFile()
    val out = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricExtension)
    new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile).execute()
    val metrics = Metric.read[InsertSizeMetric](path=out)
    metrics.length shouldBe PairOrientation.values().length

    val expectedMetrics = Seq(
      InsertSizeMetric(
        pair_orientation = PairOrientation.FR,
        read_pairs         = 2,
        standard_deviation = 1.414214,
        mean               = 2503,
        min                = 1,
        max                = 1,
        median             = 2503,
        median_absolute_deviation = 1
      ),
      InsertSizeMetric(
        pair_orientation = PairOrientation.RF,
        read_pairs         = 2,
        standard_deviation = 1.414214,
        mean               = 2305,
        min                = 1,
        max                = 1,
        median             = 2305,
        median_absolute_deviation = 1
      ),
      InsertSizeMetric(
        pair_orientation = PairOrientation.TANDEM,
        read_pairs         = 2,
        standard_deviation = 1.414214,
        mean               = 2404,
        min                = 1,
        max                = 1,
        median             = 2404,
        median_absolute_deviation = 1
      )
    )

    metrics.zip(expectedMetrics).foreach { case (actual, expected) =>
      actual shouldBe expected
    }

    val histogramPath = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricHistogramExtension)
    Io.readLines(path=histogramPath).mkString("\n") shouldBe
      """
        |insert_size	fr	rf	tandem
        |2304	0	1	0
        |2306	0	1	0
        |2403	0	0	1
        |2405	0	0	1
        |2502	1	0	0
        |2504	1	0	0
      """.stripMargin.trim
  }

  /** Developer Note (Nils Homer Jan 19 2017)
    *
    * The tests below are kept here for now for added test coverage, but can be removed at a later date if they become
    * difficult to maintain.
    */

  it should "run end-to-end and ignore reads that overlap multiple genes" in {
    val builder = new SamRecordSetBuilder()
    builder.addPair(contig=3, start1=133801671, start2=133801671) // overlaps ACKR4-4-1 and ACKR4-4-2
    val bam = builder.toTempFile()
    val out = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricExtension)
    new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile, minimumOverlap=0.0).execute()
    val metrics = Metric.read[InsertSizeMetric](path=out)
    metrics.length shouldBe PairOrientation.values().length
    metrics.zip(EmptyMetrics).foreach { case (actual, expected) =>
        actual shouldBe expected
    }
  }

  it should "run end-to-end and ignore a reads that are not fully enclosed in a gene" in {
    val builder = new SamRecordSetBuilder()
    builder.addPair(contig=2, start1=1, start2=133801671) // before ACKR4-3
    builder.addPair(contig=2, start1=133801671, start2=133814175) // after ACKR4-3
    val bam = builder.toTempFile()
    val out = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricExtension)
    new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile, minimumOverlap=0.0).execute()
    val metrics = Metric.read[InsertSizeMetric](path=out)
    metrics.length shouldBe PairOrientation.values().length
    metrics.zip(EmptyMetrics).foreach { case (actual, expected) =>
      actual shouldBe expected
    }
  }

  it should "run end-to-end and ignore reads when the insert size is different across transcripts" in {
    val builder = new SamRecordSetBuilder()
    builder.addPair(contig=4, start1=133801671, start2=133804077) // overlaps ACKR4-5 (multiple transcripts)
    val bam = builder.toTempFile()
    val out = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricExtension)
    new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile, minimumOverlap=0.0).execute()
    val metrics = Metric.read[InsertSizeMetric](path=out)
    metrics.length shouldBe PairOrientation.values().length
    metrics.zip(EmptyMetrics).foreach { case (actual, expected) =>
      actual shouldBe expected
    }
  }

  it should "run end-to-end and ignore reads when there are too few mapped bases overlapping exonic sequence" in {
    val builder = new SamRecordSetBuilder()
    builder.addPair(contig=5, start1=133801671, start2=133804077, strand1=Plus, strand2=Plus) // overlaps ACKR4-6 by 2 bases!
    val bam = builder.toTempFile()
    val out = PathUtil.pathTo(PathUtil.removeExtension(bam) + EstimateRnaSeqInsertSize.RnaSeqInsertSizeMetricExtension)

    // OK
    {
      new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile, minimumOverlap=2/200.0).execute()
      val metrics = Metric.read[InsertSizeMetric](path=out)
      metrics.length shouldBe PairOrientation.values().length
      metrics.zip(EmptyMetrics).foreach { case (actual, expected) =>
        if (actual.pair_orientation == PairOrientation.TANDEM) {
          actual shouldBe InsertSizeMetric(
            pair_orientation = PairOrientation.TANDEM,
            read_pairs         = 1,
            mean               = 1,
            min                = 1,
            max                = 1,
            median             = 1
          )
        }
        else {
          actual shouldBe expected
        }
      }
    }

    // Too few bases overlapping exonic sequence
    {
      new EstimateRnaSeqInsertSize(input=bam, refFlat=RefFlatFile, minimumOverlap=3/200.0).execute()
      val metrics = Metric.read[InsertSizeMetric](path=out)
      metrics.length shouldBe PairOrientation.values().length
      metrics.zip(EmptyMetrics).foreach { case (actual, expected) =>
        actual shouldBe expected
      }
    }
  }
}
