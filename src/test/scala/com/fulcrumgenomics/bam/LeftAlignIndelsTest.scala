/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

import com.fulcrumgenomics.alignment.{Cigar, CigarElem}
import com.fulcrumgenomics.bam.api.SamOrder.Coordinate
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.util.SequenceUtil
import org.scalatest.OptionValues

class LeftAlignIndelsTest extends UnitSpec with OptionValues {

  private def testIndelPaddedAlignment(query: String, target: String, cigar: String, cigarX: String): Unit = {
    val alignment = IndelPaddedAlignment(query, target)
    withClue (f"Hint: query: $query target: $target ::") {
      alignment.cigar(useEqualsAndX = false).toString shouldBe cigar
      alignment.cigar(useEqualsAndX = true).toString shouldBe cigarX
    }
  }

  "IndelPaddedAlignment.cigar" should "return the cigar for a simple alignment" in {
    testIndelPaddedAlignment(
      "AAAAAAAAAA",
      "AAAAATAAAA",
      "10M", "5=1X4="
    )
  }

  it should "return the cigar for internal deletion" in {
    testIndelPaddedAlignment(
      "AA--AAAAAA",
      "AATTTAAAAA",
      "2M2D6M", "2=2D1X5=")
  }

  it should "return the cigar for an internal insertion" in {
    testIndelPaddedAlignment(
      "AATTAAAAAA",
      "AA--AAAAAG",
      "2M2I6M", "2=2I5=1X")
  }

  it should "return the cigar for leading and trailing soft-clips" in {
    testIndelPaddedAlignment(
      "AAAAAAAAAA",
      "XXXAAAAAAX",
      "3S6M1S", "3S6=1S")
  }

  it should "return the cigar for a combination of insertions, deletions, soft-clips, and matches" in {
    testIndelPaddedAlignment(
      "AAAAAAAAA-AA",
      "XXX-AAAAAAXX",
      "3S1I5M1D2S", "3S1I5=1D2S")
  }

  it should "return the cigar for an alignment without matches" in {
    testIndelPaddedAlignment(
      "AAAA-AA",
      "XXX-AXX",
      "3S1I1D2S", "3S1I1D2S")
  }

  private def testLeftAlign(bases1: String, bases2: String, cigar: String): Unit = {
    withClue (f"Hint: left: $bases1 right: $bases2 ::") {
      IndelPaddedAlignment(bases1, bases2).leftAligned.cigar().toString shouldBe cigar

      if (!cigar.contains(CigarOperator.S.toString)) {
        val rightIsQueryCigar = {
          val elems = Cigar(cigar).map {
            case CigarElem(CigarOperator.D, length) => CigarElem(CigarOperator.I, length)
            case CigarElem(CigarOperator.I, length) => CigarElem(CigarOperator.D, length)
            case elem => elem
          }
          Cigar(elems.toIndexedSeq)
        }
        IndelPaddedAlignment(bases2, bases1).leftAligned.cigar().toString shouldBe rightIsQueryCigar.toString
      }
    }
  }

  "IndelPaddedAlignment.leftAligned" should "not left align an alignment that has no indels" in {
    testLeftAlign(
      "AAAAAAAAAA",
      "AAAAATAAAA",
      "10M",
    )
  }

  it should "not modify an alignment that has a leading left-aligned indel" in {
    testLeftAlign(
      "AAAAAAAAAA",
      "-AAAAAAAAA",
      "1I9M",
    )
  }

  it should "left align a 1bp indel" in {
    testLeftAlign(
      "AAAAAAAAAA",
      "AAAAA-AAAA",
      "1I9M",
    )
  }

  it should "left align a trailing 1bp indel" in {
    testLeftAlign(
      "AAAAAAAAAA",
      "AAAAAAAAA-",
      "1I9M"
    )
  }

  it should "left align and join two 1bp indel" in {
    testLeftAlign(
      "AAAAAAAAAA",
      "AA-AAAAAA-",
      "2I8M"
    )
  }

  it should "not modify an already left aligned 4bp motif indel" in {
    testLeftAlign(
      "TTTTCGCGCGCGTTTT",
      "TTTT----CGCGTTTT",
      "4M4I8M"
    )
  }

  it should "left align a 4bp motif indel" in {
    testLeftAlign(
      "TTTTTACGTACGTTTT",
      "TTTTTACG----TTTT",
      "4M4I8M"
    )
  }

  it should "left align and join two seperated 4bp motif indels" in {
    testLeftAlign(
      "TTTTTACGTACGTACGTACGTTTT",
      "TTTTTACGT----ACGT----TTT",
      "4M8I12M"
    )
  }

  it should "left align indels with leading and trailing soft-clips" in {
    testLeftAlign(
      "AAAAAAAAAAAA",
      "XA-AAAAAA-XX",
      "1S2I7M2S"
    )
  }

  private def testApply(cigar: String, query: String, target: String, cigarM: String, cigarX: String): Unit = {
    val alignment = IndelPaddedAlignment(Cigar(cigar), query, target)
    withClue (f"Hint: query: $query target: $target ::") {
      alignment.cigar(useEqualsAndX = false).toString shouldBe cigarM
      alignment.cigar(useEqualsAndX = true).toString shouldBe cigarX
    }
  }

  "IndelPaddedAlignment.apply(Cigar, String, String)" should "build a padded alignment" in {
    testApply("5M", "AAAAA", "AATAA", "5M", "2=1X2=")
    testApply("8H2S4M5I4M5D10M3S4H", "A"*28, "A"*23, "2S4M5I4M5D10M3S", "2S4=5I4=5D10=3S")
  }

  it should "build a padded alignment for a complicated cigar" in {
    def s(str: String): String = str.filter(b => SequenceUtil.isValidBase(b.toByte))
    val query   = s("HHHHHAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT----------GGGGGGGGGG----------CCCCCCCCCCTTTHH")
    val target  = s("HHHHHXXXXTTTTTTTTTT----------TTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCXXXHH")
    val cigar   = "5H4S10M10I10M10D10M10D10M3S2H"
    testApply(cigar, query, target, "4S10M10I10M10D10M10D10M3S", "4S10=10I10=10D10=10D10=3S")
  }

  "LeftIndelAligner.align" should "do nothing for a 10M alignment" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="10M", start=100).value
    aligner.align(rec, "A"*10)
    rec.cigar.toString shouldBe "10="
    rec.start shouldBe 100
  }

  it should "do nothing with a leading insertion" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="4I6M", start=100).value
    aligner.align(rec, "A"*6)
    rec.cigar.toString shouldBe "4I6="
    rec.start shouldBe 100
  }

  it should "do nothing with a leading deletion" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="4D10M", start=100).value
    aligner.align(rec, "A"*14)
    rec.cigar.toString shouldBe "4D10="
    rec.start shouldBe 100
  }

  it should "left-align an insertion to the start" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="6M4I", start=100).value
    aligner.align(rec, "A"*6)
    rec.cigar.toString shouldBe "4I6="
    rec.start shouldBe 100
  }

  it should "left-align an deletion to the start" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="10M4D", start=100).value
    aligner.align(rec, "A"*14)
    rec.cigar.toString shouldBe "4D10="
    rec.start shouldBe 100
  }

  it should "preserve hard-clips" in {
    val builder = new SamBuilder(readLength=10)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)
    val rec     = builder.addFrag(bases="A"*10, cigar="5H10M4H", start=100).value
    aligner.align(rec, "A"*10)
    rec.cigar.toString shouldBe "5H10=4H"
    rec.start shouldBe 100
  }

  it should "left-align a complicated case" in {
    def s(str: String): String = str.filter(b => SequenceUtil.isValidBase(b.toByte))
    val query   = s("HHHHHAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT----------GGGGGGGGGG----------CCCCCCCCCCTTTHH")
    val target  = s("HHHHHXXXXTTTTTTTTTT----------TTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCXXXHH")
    val cigar   = "5H4S10M10I10M10D10M10D10M3S2H"
    val builder = new SamBuilder(readLength=query.length)
    val aligner = new LeftIndelAligner(useEqualsAndX=true)

    val rec     = builder.addFrag(bases=query, cigar=cigar, start=100).value
    aligner.align(rec, target)
    rec.cigar.toString shouldBe "5H4S10I20=20D20=3S2H"
    rec.start shouldBe 100
  }

  "LeftAlignIndels" should "run end to end" in { }
  Seq(SamOrder.Coordinate, SamOrder.Queryname, SamOrder.RandomQuery, SamOrder.Random).foreach { samOrder =>
    Seq(true, false).foreach { useRef =>
      val testName = f"on a $samOrder sorted BAM " + (if (useRef) "using the ref" else "using the MD tag")
      it should f"run end to end $testName" in {
        val samBuilder       = new SamBuilder(readLength=100, sort=Some(samOrder))
        val output           = makeTempFile("test.", ".bam")

        // Build the reference
        val referenceBuilder = new ReferenceSetBuilder()
        val refSeq = referenceBuilder.add("chr1").add("A", 1000)
        val ref = if (!useRef) None else Some(referenceBuilder.toTempFile())

        // Build the reads
        {
          // frags
          samBuilder.addFrag(bases="A"*100, contig=0, start=10, cigar="100M")
          samBuilder.addFrag(bases="A"*100, contig=0, start=11, cigar="10M10D10M10I70M")
          // pairs
          samBuilder.addPair(contig=0, bases1="A"*100, start1=12, cigar1="10M10D10M10I70M", bases2="A"*100, start2=13, cigar2="100M")
          // add the MD tag
          val refBases = refSeq.bases.toArray.map(_.toByte)
          samBuilder.foreach { rec =>
            SequenceUtil.calculateMdAndNmTags(rec.asSam, refBases, true, false)
          }
        }

        val tool = new LeftAlignIndels(input=samBuilder.toTempFile(), output=output, ref=ref, useEqualsAndX=true)
        executeFgbioTool(tool)

        val records = readBamRecs(output).sortBy(Coordinate.sortkey)

        def checkFrag(rec: SamRecord, start: Int, cigar: String): Unit = {
          withClue(f"${rec.name} start(${rec.start}, $start) cigar(${rec.cigar}, $cigar)"){
            rec.start shouldBe start
            rec.cigar.toString shouldBe cigar
          }
        }

        def checkPair(rec: SamRecord, start: Int, cigar: String, mateStart: Int, mateCigar: String): Unit = {
          checkFrag(rec, start, cigar)
          rec.mateStart shouldBe mateStart
          rec.mateCigar.value.toString shouldBe mateCigar
        }

        checkFrag(records(0), 10, "100=")
        checkFrag(records(1), 11, "10D10I90=")
        checkPair(records(2), 12, "10D10I90=", 13, "100=")
        checkPair(records(3), 13, "100=", 12, "10D10I90=")
      }
    }
  }
}
