/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.testing.UnitSpec

class VcfHeaderTest extends UnitSpec {
  private val header = VcfHeader(
    contigs = IndexedSeq(VcfContigHeader(index=0, name="chr1", length=Some(100))),
    infos   = Seq(VcfInfoHeader("I1", VcfCount.Fixed(1), VcfFieldType.Integer, "info 1"), VcfInfoHeader("I2", VcfCount.Fixed(1), VcfFieldType.Integer, "info 2")),
    formats = Seq(VcfFormatHeader("F1", VcfCount.Fixed(1), VcfFieldType.Integer, "format 1"), VcfFormatHeader("F2", VcfCount.Fixed(1), VcfFieldType.Integer, "format 2")),
    filters = Seq(VcfFilterHeader("X1", "filter 1")),
    others  = Seq(VcfGeneralHeader("ALT", "DEL")),
    samples = IndexedSeq("s1")
  )

  private val newInfo   = VcfInfoHeader("I1", VcfCount.Unknown, VcfFieldType.String, "new info 1")
  private val newFormat = VcfFormatHeader("F3", VcfCount.Fixed(1), VcfFieldType.Float, "format 3")
  private val newFilter = VcfFilterHeader("X1", "new filter 1")

  "VcfHeader.withEntries" should "add INFO, FORMAT and FILTER entries, replacing any of the same kind and ID" in {
    val updated = header.withEntries(Seq(newInfo, newFormat, newFilter))
    updated.infos.map(_.description) should contain theSameElementsInOrderAs Seq("info 2", "new info 1")
    updated.formats.map(_.id) should contain theSameElementsInOrderAs Seq("F1", "F2", "F3")
    updated.filters should contain theSameElementsInOrderAs Seq(newFilter)
    updated.contigs shouldBe header.contigs
    updated.others shouldBe header.others
    updated.samples shouldBe header.samples
    updated.info("I1") shouldBe newInfo
  }

  it should "not replace an entry of a different kind with the same ID" in {
    val format = VcfFormatHeader("I1", VcfCount.Fixed(1), VcfFieldType.Integer, "format I1")
    val updated = header.withEntries(Seq(format))
    updated.info("I1").description shouldBe "info 1"
    updated.format("I1") shouldBe format
  }

  it should "fail if given a contig or other entry" in {
    Seq(VcfContigHeader(index=1, name="chr2"), VcfGeneralHeader("ALT", "INS")).foreach { entry =>
      an[IllegalArgumentException] should be thrownBy header.withEntries(Seq(entry))
      an[IllegalArgumentException] should be thrownBy header.withEntriesIfAbsent(Seq(entry))
    }
  }

  it should "fail if given two entries of the same kind with the same ID" in {
    val other = VcfInfoHeader("I1", VcfCount.Fixed(1), VcfFieldType.Float, "another info 1")
    Seq[VcfHeader => VcfHeader](_.withEntries(Seq(newInfo, other)), _.withEntriesIfAbsent(Seq(newInfo, other))).foreach { f =>
      intercept[IllegalArgumentException] { f(header) }.getMessage should include ("ID I1")
    }
    // the same ID on entries of different kinds is fine
    header.withEntries(Seq(newInfo, VcfFormatHeader("I1", VcfCount.Fixed(1), VcfFieldType.Integer, "format I1"))).format.contains("I1") shouldBe true
  }

  "VcfHeader.withEntriesIfAbsent" should "add only the entries whose ID is not already used by an entry of the same kind" in {
    val updated = header.withEntriesIfAbsent(Seq(newInfo, newFormat, newFilter))
    updated.infos shouldBe header.infos
    updated.filters shouldBe header.filters
    updated.formats.map(_.id) should contain theSameElementsInOrderAs Seq("F1", "F2", "F3")
  }

  "VcfContigHeader" should "convert to and from a sequence in a sequence dictionary" in {
    val contig = VcfContigHeader(index=3, name="chr4", length=Some(100), assembly=Some("hg38"), md5=Some("abc"), species=Some("human"), url=Some("http://x"))
    val metadata = contig.toSequenceMetadata
    metadata.index shouldBe 3
    metadata.length shouldBe 100
    metadata.uri shouldBe Some("http://x")
    VcfContigHeader.fromSequenceMetadata(metadata) shouldBe contig

    val noLength = VcfContigHeader(index=0, name="chr1")
    noLength.toSequenceMetadata.length shouldBe SequenceMetadata.UnknownSequenceLength
    VcfContigHeader.fromSequenceMetadata(noLength.toSequenceMetadata) shouldBe noLength
  }
}
