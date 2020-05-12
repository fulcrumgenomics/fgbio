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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.fasta.Converters.{FromSAMSequenceDictionary, FromSAMSequenceRecord, ToSAMSequenceDictionary, ToSAMSequenceRecord}
import com.fulcrumgenomics.testing.UnitSpec
import org.scalatest.OptionValues
import com.fulcrumgenomics.fasta.SequenceMetadata.{AlternateLocus, Keys}
import com.fulcrumgenomics.fasta.Topology.{Circular, Linear}
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}
import com.fulcrumgenomics.FgBioDef._

class SequenceDictionaryTest extends UnitSpec with OptionValues {

  private val validSequence = SequenceMetadata(
    name             = "chr1",
    length           = 100,
    aliases          = Seq( "chr1_1", "chr1_2", "chr1_3"),
    assembly         = Some("assembly"),
    description      = Some("description"),
    md5              = Some("123"),
    species          = Some("species"),
    customAttributes = Map(
      "key1" -> "value1",
      "key2" -> "value2"
    )
  )

  private val emptySequence  = SequenceMetadata(name="chr1", length=123)

  "SequenceMetadata" should "fail if not a validated sequence" in {
    an[Exception] should be thrownBy SequenceMetadata(name="(((242*(@&$)))", length=12)
  }

  it should "fail if the length is less than zero" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=(-1))
  }

  it should "fail if the sequence name tag is given in the attributes" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=1, customAttributes=Map("SN" -> "chr1"))
  }

  it should "fail if the sequence length tag is given in the attributes" in {
    an[Exception] should be thrownBy SequenceMetadata(name="chr1", length=1, customAttributes=Map("LN" -> "2"))
  }

  it should "implement apply, get, and contains" in {
    // apply
    validSequence("key1") shouldBe "value1"
    an[Exception] should be thrownBy validSequence("key3")

    // get
    validSequence.get("key1").value shouldBe "value1"
    validSequence.get("key3") shouldBe 'empty

    // contains
    validSequence.contains("key1") shouldBe true
    validSequence.contains("key3") shouldBe false
  }

  "SequenceMetadata.allNames" should "return all names including aliases" in {
    validSequence.allNames should contain theSameElementsInOrderAs Seq("chr1", "chr1_1", "chr1_2", "chr1_3")
  }

  "SequenceMetadata.aliases" should "return all names including aliases" in {
    validSequence.aliases should contain theSameElementsInOrderAs Seq("chr1_1", "chr1_2", "chr1_3")
  }

  "SequenceMetadata.isAlternate" should "return if the sequence is an alternate locus" in {
    val noAltLocus  = SequenceMetadata(name="chr1", length=123)
    val altStar     = new SequenceMetadata(name="chr1", length=123, index=0, attributes=Map("AH" -> "*"))
    val locusEquals = new SequenceMetadata(name="chr1", length=123, index=0, attributes=Map("AH" -> "=:1-2"))
    val locusFull   = new SequenceMetadata(name="chr1", length=123, index=0, attributes=Map("AH" -> "chr4:2-3"))

    noAltLocus.isAlternate shouldBe false
    noAltLocus.alternate shouldBe 'empty

    altStar.isAlternate shouldBe false
    altStar.alternate shouldBe 'empty

    locusEquals.isAlternate shouldBe true
    locusEquals.alternate.value shouldBe AlternateLocus(refName="chr1", start=1, end=2)

    locusFull.isAlternate shouldBe true
    locusFull.alternate.value shouldBe AlternateLocus(refName="chr4", start=2, end=3)
  }

  "SequenceMetadata.md5" should "return the md5 checksum if present" in {
    validSequence.md5.value shouldBe "123"
    validSequence.md5Int.value shouldBe BigInt("123", 16)
    emptySequence.md5 shouldBe 'empty
    emptySequence.md5Int shouldBe 'empty
  }

  "SequenceMetadata.assembly" should "return the assembly if present" in {
    validSequence.assembly.value shouldBe "assembly"
    emptySequence.assembly shouldBe 'empty
  }

  "SequenceMetadata.species" should "return the species if present" in {
    validSequence.species.value shouldBe "species"
    emptySequence.species shouldBe 'empty
  }

  "SequenceMetadata.description" should "return the description if present" in {
    validSequence.description.value shouldBe "description"
    emptySequence.description shouldBe 'empty
  }

  "SequenceMetadata.topology" should "return the topology if present" in {
    val linear   = SequenceMetadata(name="chr1", length=123, topology=Some(Topology.Linear))
    val circular = SequenceMetadata(name="chr1", length=123,  topology=Some(Topology.Circular))

    linear.topology.value shouldBe Linear
    circular.topology.value shouldBe Circular
    emptySequence.topology shouldBe 'empty
  }

  "SequenceMetadata.sameAs" should "return false if not reference names are shared" in {
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr2", length=0) shouldBe false
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr2", length=0, aliases=Seq("chr3")) shouldBe false
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr2", length=0, aliases=Seq("chr3", "chr4")) shouldBe false
  }

  it should "return false if the lengths are different" in {
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr1", length=1) shouldBe false
  }


  it should "return false if the md5s are different" in {
    SequenceMetadata(name="chr1", length=0, md5=Some("1")) sameAs SequenceMetadata(name="chr1", length=0, md5=Some("2")) shouldBe false
  }

  it should "return true if the primary name is the same" in {
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr1", length=0) shouldBe true
  }

  it should "return true if the primary name in one matches an alias in the other" in {
    SequenceMetadata(name="chr1", length=0) sameAs SequenceMetadata(name="chr2", length=0, aliases=Seq("chr1")) shouldBe true
    SequenceMetadata(name="chr2", length=0, aliases=Seq("chr1")) sameAs SequenceMetadata(name="chr1", length=0) shouldBe true
  }

  it should "return true if they share an alias" in {
    SequenceMetadata(name="chrX", length=0, aliases=Seq("chr1")) sameAs SequenceMetadata(name="chrY", length=0, aliases=Seq("chr1")) shouldBe true
    SequenceMetadata(name="chrX", length=0, aliases=Seq("chr1", "chr2")) sameAs SequenceMetadata(name="chrY", length=0, aliases=Seq("chr2", "chr3")) shouldBe true
  }

  it should "return true if the name, length, and md5 match" in {
    SequenceMetadata(name="chr1", length=0, md5=Some("1")) sameAs SequenceMetadata(name="chr1", length=0, md5=Some("1")) shouldBe true
  }

  "SequenceDictionary" should "fail to build two sequence metadatas share a name (including aliases)" in {
    val one   = SequenceMetadata(name="chr1", length=0)
    val two   = SequenceMetadata(name="chr2", length=0, aliases=Seq("chr1"))
    val three = SequenceMetadata(name="chr3", length=0, aliases=Seq("chr0", "chr1", "chrX"))

    an[Exception] should be thrownBy SequenceDictionary(one, one)
    an[Exception] should be thrownBy SequenceDictionary(one, two)
    an[Exception] should be thrownBy SequenceDictionary(two, three)
  }

  it should "implement apply, get, contains, and iterator" in {
    val one   = SequenceMetadata(name="1", length=0)
    val two   = SequenceMetadata(name="2", length=0, aliases=Seq("1_alt"))
    val three = SequenceMetadata(name="3", length=0, aliases=Seq("3_alt_1", "3_alt_2","3_alt_3"))
    val dict  = SequenceDictionary(one, two, three)

    dict.infos.length shouldBe 3
    dict.infos.zip(Seq(one, two, three)).foreach { case (left, right) => left sameAs right shouldBe true }
    dict.iterator.zip(Iterator(one, two, three)).foreach { case (left, right) => left sameAs right shouldBe true }
    dict.zip(Seq(one, two, three)).foreach { case (left, right) => left sameAs right shouldBe true }


    dict.infos.zipWithIndex.foreach { case (seq, index) =>
      seq.index shouldBe index
      dict(index) shouldBe seq
      seq.allNames.foreach { name =>
        dict(name) shouldBe seq
        dict.get(name).value shouldBe seq
        dict.contains(name) shouldBe true
      }
    }

    an[IndexOutOfBoundsException] should be thrownBy dict.apply(3)
    an[Exception] should be thrownBy dict("4")
    dict.get("4") shouldBe 'empty
    dict.contains("4") shouldBe false
  }

  "Converters.ToSequenceRecord" should "convert from a [[SequenceMetadata]] to a [[SAMSequenceRecord]]" in {
    val record: SAMSequenceRecord = new ToSAMSequenceRecord(validSequence).asSam
    record.getSequenceName shouldBe validSequence.name
    record.getSequenceLength shouldBe validSequence.length
    record.getSequenceIndex shouldBe validSequence.index
    record.getAssembly shouldBe validSequence.assembly.value
    record.getDescription shouldBe validSequence.description.value
    record.getMd5 shouldBe validSequence.md5.value
    record.getSpecies shouldBe validSequence.species.value

    validSequence.attributes.foreach { case (key, value) => record.getAttribute(key) shouldBe value }
  }

  "Converters.FromSequenceRecord" should "convert from a [[SAMSequenceRecord]] to a [[SequenceMetadata]]" in {
    // Note: relies on the above test that converts from a [[SequenceMetadata]] to a [[SAMSequenceRecord]]
    val record: SequenceMetadata = new FromSAMSequenceRecord(new ToSAMSequenceRecord(validSequence).asSam).fromSam
    record shouldBe validSequence
    validSequence.attributes.foreach { case (key, value) => record(key) shouldBe value }
  }

  "Converters.ToSequenceDictionary" should "convert from a [[SequenceDictionary]] to a [[SAMSequenceDictionary]]" in {
    val scalaDict: SequenceDictionary = SequenceDictionary(
      validSequence,
      SequenceMetadata(name="chrX", length=234),
      SequenceMetadata(name="chrY", length=345)
    )
    val javaDict: SAMSequenceDictionary = new ToSAMSequenceDictionary(scalaDict).asSam

    javaDict.getReferenceLength shouldBe scalaDict.infos.map(_.length).sum
    javaDict.getSequences.length shouldBe scalaDict.length
    javaDict.getSequences.iterator().zip(scalaDict.iterator).foreach {
      case (javaRecord: SAMSequenceRecord, scalaRecord: SequenceMetadata) =>
        // Note: relies on the above test that converts from a [[SequenceMetadata]] to a [[SAMSequenceRecord]]
        new FromSAMSequenceRecord(javaRecord).fromSam shouldBe scalaRecord
    }
  }

  "Converters.FromSAMSequenceDictionary" should "convert from a [[SAMSequenceDictionary]] to a [[SequenceDictionary]]" in {
    // Note: relies on the above tests that converts from a [[SequenceDictionary]] to a [[SAMSequenceDictionary]]
    val sourceDict: SequenceDictionary = SequenceDictionary(
      validSequence,
      SequenceMetadata(name="chrX", length=234),
      SequenceMetadata(name="chrY", length=345)
    )
    val targetDict: SequenceDictionary = new FromSAMSequenceDictionary(new ToSAMSequenceDictionary(sourceDict).asSam).fromSam

    sourceDict shouldBe targetDict
  }
}
