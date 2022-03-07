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
 *
 */

package com.fulcrumgenomics.testing

import java.nio.file.Files

import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import org.scalatest.OptionValues

class ReferenceSetBuilderTest extends UnitSpec with OptionValues {
  "ReferenceSetBuilder" should "write a simple FASTA file, with accompanying .dict and .fai files" in {
    val builder = new ReferenceSetBuilder
    builder.add("chr1").add("AAAAAAAAAA", 500) // 5000 bases
    builder.add("chr2").add("CCCCCCCCCC", 500) // 5000 bases
    val fasta = builder.toTempFile()

    // Check the .dict file
    val dictIn = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta)
    Files.exists(dictIn) shouldBe true
    val dict = SequenceDictionary.extract(dictIn)

    // Check the .fai file
    val fai =  ReferenceSequenceFileFactory.getFastaIndexFileName(fasta)
    Files.exists(fai) shouldBe true

    // Read it back in
    val ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta)
    ref.isIndexed shouldBe true // since the .fai file exists
    ref.getSequenceDictionary.getSequences.size shouldBe 2
    ref.getSequence("chr1").getBases.count(_ == 'A') shouldBe 5000
    ref.getSequence("chr2").getBases.count(_ == 'C') shouldBe 5000
    ref.getSequenceDictionary.fromSam.sameAs(dict)
    ref.close()
  }

  it should "write an assembly name into the dictionary if given one" in {
    val Seq(plus, minus) = Seq(Some("hg19"), None).map { assembly =>
      val builder = new ReferenceSetBuilder(assembly=assembly)
      builder.add("chr1").add("A", 100)
      val path = builder.toTempFile()
      val dict = SequenceDictionary.extract(path)
      dict.length shouldBe 1
      dict
    }

    plus("chr1").assembly.value shouldBe "hg19"
    minus("chr1").assembly shouldBe empty
  }

  it should "write species name into the dictionary if given one" in {
    val Seq(plus, minus) = Seq(Some("human"), None).map { species =>
      val builder = new ReferenceSetBuilder(species=species)
      builder.add("chr1").add("A", 100)
      val path = builder.toTempFile()
      val dict = SequenceDictionary.extract(path)
      dict.length shouldBe 1
      dict
    }

    plus("chr1").species.value shouldBe "human"
    minus("chr1").species shouldBe empty
  }

  it should "calculate MD5s for sequences only if asked" in {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1").add("A", 100)
    val Seq(plus, minus) = Seq(true, false).map(md5 => builder.toTempFile(calculateMds5=md5)).map(SequenceDictionary.extract)
    plus("chr1").md5 shouldBe defined
    minus("chr1").md5 shouldBe empty
  }
}
