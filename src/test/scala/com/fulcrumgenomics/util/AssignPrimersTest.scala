/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.SamBuilder.Minus
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SamPairUtil
import org.scalatest.OptionValues

class AssignPrimersTest extends UnitSpec with OptionValues {
  private def checkMatch(rec: SamRecord, amplicon: Option[Amplicon] = None, mateAmplicon: Option[Amplicon] = None): Unit = {
    amplicon match {
      case None =>
        rec.contains(AssignPrimers.PrimerCoordinateTag) shouldBe false
        rec.contains(AssignPrimers.AmpliconIdentifierTag) shouldBe false
      case Some(amp) =>
        if (rec.positiveStrand) {
          rec.get[String](AssignPrimers.PrimerCoordinateTag).value shouldBe amp.leftPrimerLocation
        }
        else {
          rec.get[String](AssignPrimers.PrimerCoordinateTag).value shouldBe amp.rightPrimerLocation
        }
        rec.get[String](AssignPrimers.AmpliconIdentifierTag).value shouldBe amp.identifier
    }

    mateAmplicon match {
      case None =>
        rec.contains(AssignPrimers.MatePrimerCoordinateTag) shouldBe false
        rec.contains(AssignPrimers.MateAmpliconIdentifierTag) shouldBe false
      case Some(amp) =>
        rec.paired shouldBe true
        if (rec.matePositiveStrand) {
          rec.get[String](AssignPrimers.MatePrimerCoordinateTag).value shouldBe amp.leftPrimerLocation
        }
        else {
          rec.get[String](AssignPrimers.MatePrimerCoordinateTag).value shouldBe amp.rightPrimerLocation
        }
        if (amplicon.contains(amp)) { // same amplicon
          rec.get[String](AssignPrimers.MateAmpliconIdentifierTag).value shouldBe "="
        }
        else {
          rec.get[String](AssignPrimers.MateAmpliconIdentifierTag).value shouldBe amp.identifier
        }
    }
  }
  private def checkMatch(rec: SamRecord, amplicon: Amplicon): Unit = checkMatch(rec=rec, amplicon=Some(amplicon))
  private def checkMatch(rec: SamRecord, amplicon: Amplicon, mateAmplicon: Amplicon): Unit = {
    checkMatch(rec=rec, amplicon=Some(amplicon), mateAmplicon=Some(mateAmplicon))
  }

  {
    val ampliconsPath = makeTempFile("amplicons.", ".txt")
    val outputBam     = makeTempFile("output.", ".bam")
    val metricsPath   = makeTempFile("metrics.", ".txt")

    // The amplicons to match against
    val amplicons = Seq(
      Amplicon("chr1", 100, 120, 180, 200),
      Amplicon("chr1", 500, 520, 580, 600),
      Amplicon("chr5", 100, 120, 180, 200)
    )

    // The reads to assign
    val builder = new SamBuilder()
    builder.addFrag(start=100, cigar="100M") // matches amplicon #1
    builder.addFrag(start=100, cigar="50S50M", strand=SamBuilder.Minus) // no match, wrong strand
    builder.addPair(start1=500, cigar1="100M", start2=501, cigar2="100M") // match amplicon #2 R1->F, R2->R
    builder.addPair(start1=501, cigar1="100M", strand1=SamBuilder.Minus, start2=500, cigar2="100M", strand2=SamBuilder.Plus) // match amplicon #2 R1->R, R2->F
    builder.addPair(start1=100, start2=501) // match amplicon #1 and #2, R1->1F, R2->2R
    builder.addPair(start1=100, start2=10000) // match amplicon #1 for R1->1F, no match for R2
    builder.addPair(contig=4, start1=100, start2=101) // match amplicon #5 , R1->F, R2->R
    // Reads to test --annotate-all, where we have a primary and supplementary read
    // R1->1F, R2_primary->2R, R2_supplementary->NA.  If --annotate-all is used, then R2_supplementary should have 2R
    val Seq(matePrimary, _) = builder.addPair(start1=100, start2=551, cigar1="100M", cigar2="50S50M")
    builder.addFrag(name=matePrimary.name, contig=2, start=1, cigar="50M50S", strand=Minus).foreach { supplemental =>
      supplemental.paired        = true
      supplemental.firstOfPair   = false
      supplemental.secondOfPair  = true
      supplemental.supplementary = true
      SamPairUtil.setMateInformationOnSupplementalAlignment(supplemental.asSam, matePrimary.asSam, true)
    }

    Metric.write(ampliconsPath, amplicons)

    // Test with/without --annotate-all set
    Seq(false, true).foreach { annotateAll =>
      val label = if (annotateAll) " with --annotate-all" else ""
      "AssignPrimers" should s"assign primers to reads$label" in {
        val tool = new AssignPrimers(
          input       = builder.toTempFile(),
          primers     = ampliconsPath,
          output      = outputBam,
          metrics     = metricsPath,
          annotateAll = annotateAll
        )

        tool.execute()

        // check the output BAM
        val outputRecs = readBamRecs(bam=outputBam)
        // The output BAM may be resorted, so the following methods retrieve records by read name (since they're numerically indexed)
        def recs(nameIndex: Int): Seq[SamRecord] = outputRecs.filter(_.name.toInt == nameIndex)
        def rec(nameIndex: Int): SamRecord = recs(nameIndex) match {
          case Seq(rec) => rec
          case records  => throw new IllegalStateException(s"Found ${records.length}")
        }

        // frag
        checkMatch(rec=rec(0), amplicons(0))
        // frag
        checkMatch(rec=rec(1))
        // pair #1
        checkMatch(rec=recs(2).head, amplicons(1), amplicons(1))
        checkMatch(rec=recs(2).last, amplicons(1), amplicons(1))
        // pair #2
        checkMatch(rec=recs(3).head, amplicons(1), amplicons(1))
        checkMatch(rec=recs(3).last, amplicons(1), amplicons(1))
        // pair #3
        checkMatch(rec=recs(4).head, amplicons(0), amplicons(1))
        checkMatch(rec=recs(4).last, amplicons(1), amplicons(0))
        // pair #4
        checkMatch(rec=recs(5).head, amplicons(0))
        checkMatch(rec=recs(5).last, None, Some(amplicons(0)))
        // pair #5
        checkMatch(rec=recs(6).head, amplicons(2), amplicons(2))
        checkMatch(rec=recs(6).last, amplicons(2), amplicons(2))
        // pair #6
        checkMatch(rec=recs(7).head, amplicons(0), amplicons(1))
        checkMatch(rec=recs(7)(1), amplicons(1), amplicons(0))
        if (annotateAll) { // the amplicon tag (ma) **is** set in this case
          checkMatch(rec=recs(7)(2), amplicons(1), amplicons(0))
        } else { // the amplicon tag (ma) **is not** set in this case
          checkMatch(rec=recs(7)(2), None, Some(amplicons(0)))
        }

        // check the metrics
        val metrics       = Metric.read[AssignPrimersMetric](metricsPath)
        val allIdentifier = AssignPrimersMetric.AllAmpliconsIdentifier
        val expected      = IndexedSeq(
          AssignPrimersMetric(identifier=amplicons(0).identifier, left=4, right=0, r1s=4, r2s=0, pairs=0),
          AssignPrimersMetric(identifier=amplicons(1).identifier, left=2, right=4, r1s=2, r2s=4, pairs=2),
          AssignPrimersMetric(identifier=amplicons(2).identifier, left=1, right=1, r1s=1, r2s=1, pairs=1),
          AssignPrimersMetric(identifier=allIdentifier,           left=7, right=5, r1s=7, r2s=5, pairs=3)
        ).map(_.finalize(total=outputRecs.length))

        metrics.length shouldBe amplicons.length + 1
        metrics.zip(expected).foreach { case (act, exp) =>
          // compare one-by-one to make it easier to debug
          act.finalize(0) shouldBe exp.finalize(0)  // just the values
          act shouldBe exp // now the fractions
        }
        metrics.length shouldBe expected.length

      }
    }
  }
}
