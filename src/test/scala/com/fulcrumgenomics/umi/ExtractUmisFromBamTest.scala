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
 *
 */

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.bam.api.{SamRecord, SamWriter}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.ReadStructure
import org.scalatest.OptionValues

/**
  * Tests for ExtractUmisFromBam.
  */
class ExtractUmisFromBamTest extends UnitSpec with OptionValues {

  def newBam = makeTempFile("extract_umis_from_bam_test.", ".bam")

  /** An implicit conversion to read stucture to make test writing easier below by allowing the use of Strings
    * directly instead of having to write ReadStructure("75T") everywhere.
    */
  private implicit def stringToReadStructure(rs: String): ReadStructure = ReadStructure(rs)

  def annotateRecordFragment: SamRecord = {
    new SamBuilder(readLength=100).addFrag(name="Frag", start=1).map {rec => rec.bases = "A" * 100; rec }.get
  }

  "ExtractUmisFromBam.annotateRecord" should "should not annotate with no molecular indices" in {
    val frag = annotateRecordFragment
    ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("100T"), Seq.empty) shouldBe 'empty
    frag.length shouldBe 100
  }

  it should "throw an exception when fewer molecular indices in the read structure than molecular SAM tags are given then " in {
    val tags = Seq("AB", "CD")
    val frag = annotateRecordFragment
    an[Exception] should be thrownBy ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("10M80T"), tags)
  }

  it should "throw an exception when fewer molecular SAM tags are given than molecular indices in the read structure" in {
    val tags = Seq("AB", "CD")
    val frag = annotateRecordFragment
    an[Exception] should be thrownBy ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("10M10M10M70T"), tags)
  }

  "ExtractUmisFromBam" should "should annotate a single molecular index for a read pair" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec.bases = base * 100
    }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("A1")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec[String]("A1") shouldBe "AAAAAAAAAA-CCCCCCCCCC"
      rec.length shouldBe 90
      rec.basesString shouldBe base * 90
    }
  }

  it should "should annotate a single molecular index for a fragment read" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1, unmapped=true).map { rec => rec.bases = "A" * 100; rec }
    builder.addFrag(name="Frag2", start=1, unmapped=true).map { rec => rec.bases = "C" * 100; rec }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T"), molecularIndexTags=Seq("A1")).execute()
    val recs = readBamRecs(output)

    recs.head[String]("A1") shouldBe "AAAAAAAAAA"
    recs.head.length shouldBe 90
    recs.head.basesString shouldBe "A" * 90

    recs.last[String]("A1") shouldBe "CCCCCCCCCC"
    recs.last.length shouldBe 90
    recs.last.basesString shouldBe "C" * 90
  }

  it should "accept a single molecular index SAM tag" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("R1")).execute()
  }

  it should "should annotate a two molecular index" in {
    val builder = new SamBuilder(readLength=100)
    val frag = annotateRecordFragment
    ExtractUmisFromBam.annotateRecord(record=frag, ReadStructure("90T5M5M"), Seq("R1", "R2")) shouldBe "A" * 10
    frag.length shouldBe 90
  }

  it should "not accept a single read" in {
    val output = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept mapped fragments" in {
    val output = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1)
    builder.addFrag(name="Frag2", start=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T")).execute()
  }

  it should "not accept mapped reads" in {
    val output = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addPair(name="Pair", start1=1, start2=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }


  it should "not accept read pairs when only one read structure is given" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1)
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T")).execute()
  }

  it should "only not accept paired records with the first end before the second end" in {
    val input   = newBam
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addPair(name="Pair", start1=1, start2=1)

    val writer = SamWriter(input, builder.header)
    writer ++= builder.iterator.toSeq.reverse
    writer.close()

    an[Exception] should be thrownBy new ExtractUmisFromBam(input=input, output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept read pairs with the different names" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1)
    records.head.name = "PairNot"
    an[Exception] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T")).execute()
  }

  it should "not accept multiple molecular index SAM tags and have a different number of molecular indices in the read structure" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    an[ValidationException] should be thrownBy new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("A1", "B1", "C1")).execute()
  }

  it should "accept duplicate molecular index SAM tags and concatenate the values" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec.bases = base * 100
    }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("A1", "A1")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec[String]("A1") shouldBe "AAAAAAAAAA-CCCCCCCCCC"
      rec.length shouldBe 90
      rec.basesString shouldBe base * 90
    }
  }

  it should "accept duplicate molecular index SAM tags and concatenate the values for fragment reads" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1, unmapped=true).map { rec => rec.bases = "A" * 100; rec }
    builder.addFrag(name="Frag2", start=1, unmapped=true).map { rec => rec.bases = "C" * 100; rec }

    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M80T10M"), molecularIndexTags=Seq("A1", "A1")).execute()
    val recs = readBamRecs(output)

    recs.head[String]("A1") shouldBe "AAAAAAAAAA-AAAAAAAAAA"
    recs.head.length shouldBe 80
    recs.head.basesString shouldBe "A" * 80

    recs.last[String]("A1") shouldBe "CCCCCCCCCC-CCCCCCCCCC"
    recs.last.length shouldBe 80
    recs.last.basesString shouldBe "C" * 80
  }

  it should "extract the molecular indices and annotate them for a read pair in tags" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec => rec.bases = "A" * 100 }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "1M99T"), molecularIndexTags=Seq("P1", "P2")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      if (rec.firstOfPair) {
        rec[String]("P1") shouldBe "A" * 10
        rec.length shouldBe 90
        rec.basesString shouldBe "A" * 90
      }
      else {
        rec[String]("P2") shouldBe "A"
        rec.length shouldBe 99
        rec.basesString shouldBe "A" * 99
      }
    }
  }

  it should "extract the molecular indices and annotate them for a read pair in tags and in the names" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec => rec.bases = "A" * 100 }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "1M99T"), molecularIndexTags=Seq("P1", "P2"), annotateReadNames=true).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      if (rec.firstOfPair) {
        rec[String]("P1") shouldBe "A" * 10
        rec.length shouldBe 90
        rec.basesString shouldBe "A" * 90
        rec.name shouldBe "Pair+" + "A" * 11
      }
      else {
        rec[String]("P2") shouldBe "A"
        rec.length shouldBe 99
        rec.basesString shouldBe "A" * 99
        rec.name shouldBe "Pair+" + "A" * 11
      }
    }
  }

  it should "extract the molecular indices and annotate them for a fragment reads in tags and in the names" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addFrag(name="Frag1", start=1, unmapped=true).map { rec => rec.bases = "A" * 100; rec }
    builder.addFrag(name="Frag2", start=1, unmapped=true).map { rec => rec.bases = "C" * 100; rec }

    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("1M99T"), molecularIndexTags=Seq("P1"), annotateReadNames=true).execute()
    val recs = readBamRecs(output)

    recs.head[String]("P1") shouldBe "A"
    recs.head.length shouldBe 99
    recs.head.basesString shouldBe "A" * 99
    recs.head.name shouldBe "Frag1+A"

    recs.last[String]("P1") shouldBe "C"
    recs.last.length shouldBe 99
    recs.last.basesString shouldBe "C" * 99
    recs.last.name shouldBe "Frag2+C"
  }

  it should "extract the molecular indices and annotate them for a read pair in tags and in a single tag" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    builder.addPair(name="q", start1=1, start2=1, unmapped1=true, unmapped2=true).foreach(r => r.bases = "A" * 100)
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("4M96T", "6M94T"), molecularIndexTags=Seq("ZA", "ZB"), singleTag=Some("RX")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      rec[String]("ZA") shouldBe "AAAA"
      rec[String]("ZB") shouldBe "AAAAAA"
      rec[String]("RX") shouldBe "AAAA-AAAAAA"
    }
  }

  it should "not allow invalid usages of the single-tag option" in {
    val (in, out) = (newBam, newBam)

    an[ValidationException] shouldBe thrownBy {
      new ExtractUmisFromBam(input=in, output=out, molecularIndexTags=Seq("RX"), singleTag=Some("RX"), readStructure=Seq("5M5T", "10T"))
    }

    an[ValidationException] shouldBe thrownBy {
      new ExtractUmisFromBam(input=in, output=out, molecularIndexTags=Seq("RX"), singleTag=Some("RX"), readStructure=Seq("5M5T", "5M5T"))
    }

    an[ValidationException] shouldBe thrownBy {
      new ExtractUmisFromBam(input=in, output=out, molecularIndexTags=Seq("RX", "RZ"), singleTag=Some("RX"), readStructure=Seq("5M5T", "5M5T"))
    }

    an[ValidationException] shouldBe thrownBy {
      new ExtractUmisFromBam(input=in, output=out, molecularIndexTags=Seq("RX", "RZ"), singleTag=Some("RXX"), readStructure=Seq("5M5T", "5M5T"))
    }
  }

  it should "should annotate a single molecular index and update clipping information" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec.bases = base * 100
      rec("XT") =  20
    }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("A1"), clippingAttribute=Some("XT")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec[String]("A1") shouldBe "AAAAAAAAAA-CCCCCCCCCC"
      rec.length shouldBe 90
      rec.basesString shouldBe base * 90
      rec[Int]("XT") shouldBe 10
    }
  }

  it should "should annotate a single molecular index and update clipping information, even with skips" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec.bases = base * 100
      rec("XT") = 30
    }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M10S80T", "10M10S80T"), molecularIndexTags=Seq("A1"), clippingAttribute=Some("XT")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec[String]("A1") shouldBe "AAAAAAAAAA-CCCCCCCCCC"
      rec.length shouldBe 80
      rec.basesString shouldBe base * 80
      rec[Int]("XT") shouldBe 10
    }
  }

  it should "should annotate a single molecular index and update clipping information on only the first of pair" in {
    val output  = newBam
    val builder = new SamBuilder(readLength=100)
    val records = builder.addPair(name="Pair", start1=1, start2=1, unmapped1=true, unmapped2=true)
    records.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec.bases = base * 100
      if (rec.firstOfPair) rec("XT") = 20
    }
    new ExtractUmisFromBam(input=builder.toTempFile(), output=output, readStructure=Seq("10M90T", "10M90T"), molecularIndexTags=Seq("A1"), clippingAttribute=Some("XT")).execute()
    val recs = readBamRecs(output)
    recs.foreach { rec =>
      val base = if (rec.firstOfPair) "A" else "C"
      rec[String]("A1") shouldBe "AAAAAAAAAA-CCCCCCCCCC"
      rec.length shouldBe 90
      rec.basesString shouldBe base * 90
      if (rec.firstOfPair) rec[Int]("XT") shouldBe 10
      else rec.get("XT") shouldBe None
    }
  }

  "ExtractUmisFromBam.updateClippingInformation" should "update the clipping information for non-template bases" in {
    val builder = new SamBuilder(readLength=100)
    val record = builder.addFrag(contig=0, start=1).value

    // no clipping information on the record
    ExtractUmisFromBam.updateClippingInformation(record, Some("XT"), ReadStructure("100T"))
    record.get("XT") shouldBe None

    // no clipping tag given
    record("XT") = 100
    ExtractUmisFromBam.updateClippingInformation(record, None, ReadStructure("100T"))
    record[Int]("XT") shouldBe 100

    // no non-template bases
    record("XT") = 80
    ExtractUmisFromBam.updateClippingInformation(record, Some("XT"), ReadStructure("100T"))
    record[Int]("XT") shouldBe 80

    // 10-bases are skipped, so the read will be 90 bases long, and the clipping position should be shifted accordingly
    record("XT") = 80
    ExtractUmisFromBam.updateClippingInformation(record, Some("XT"), ReadStructure("20T10S70T"))
    record[Int]("XT") shouldBe 70

    // same as previous, but the clipping position is in the skip
    record("XT") = 80
    ExtractUmisFromBam.updateClippingInformation(record, Some("XT"), ReadStructure("75T10S15T"))
    record[Int]("XT") shouldBe 75
  }
}
