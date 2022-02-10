/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec

class SortSequenceDictionaryTest extends UnitSpec {
  private val dir        = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/fasta")
  private val chr1 = SequenceMetadata(
    name             = "chr1",
    length           = 100,
    aliases          = Seq( "chr1_1", "chr1_2", "chr1_3"),
  )

  private val chr2 = SequenceMetadata(
    name             = "chr2",
    length           = 100,
    aliases          = Seq( "chr2_1", "chr2_2", "chr2_3"),
  )
  private val chrM = SequenceMetadata(
    name             = "chrM",
    length           = 100,
    aliases          = Seq( "chrM_1", "chrM_2", "chrM_3"),
  )

  private val chr1Alias = SequenceMetadata(
    name   = "chr1_1",
    length = 100,
  )
  private val chr2Alias = SequenceMetadata(
    name   = "chr2_1",
    length = 100,
  )
  private val chrMAlias = SequenceMetadata(
    name   = "chrM_1",
    length = 100,
  )

  private val dictChr1Chr2ChrM = SequenceDictionary(chr1, chr2, chrM)

  "SortSequenceDictionary" should "keep output identical to input if the sort order is the same" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    dictChr1Chr2ChrM.write(inputDictionaryPath)
    dictChr1Chr2ChrM.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input         = inputDictionaryPath,
      sortDictionary = sortDictionaryPath,
      output        = outputDictionaryPath,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(dictChr1Chr2ChrM) shouldBe true
  }
  
  it should "append metadata row if contig is present in input and not sortDictionary, and `skip-missing-contigs` is false" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    val subsetDictionary   = SequenceDictionary(chr2, chrM)
    val expectedDictionary = SequenceDictionary(chr2, chrM, chr1)

    dictChr1Chr2ChrM.write(inputDictionaryPath)
    subsetDictionary.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input              = inputDictionaryPath,
      sortDictionary     = sortDictionaryPath,
      output             = outputDictionaryPath,
      skipMissingContigs = false,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(expectedDictionary) shouldBe true
  }

  it should "skip metadata row if contig is present in input and not sortDictionary, and `--skip-missing-contigs` is true" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    val subsetDictionary = SequenceDictionary(chr2, chrM)

    dictChr1Chr2ChrM.write(inputDictionaryPath)
    subsetDictionary.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input              = inputDictionaryPath,
      sortDictionary     = sortDictionaryPath,
      output             = outputDictionaryPath,
      skipMissingContigs = true,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(subsetDictionary) shouldBe true
  }

  it should "skip metadata row if contig is present in sortDictionary and not input" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    val subsetDictionary = SequenceDictionary(chr2, chrM)

    subsetDictionary.write(inputDictionaryPath)
    dictChr1Chr2ChrM.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input          = inputDictionaryPath,
      sortDictionary = sortDictionaryPath,
      output         = outputDictionaryPath,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(subsetDictionary) shouldBe true
  }

  it should "reorder contigs if sortDictionary is provided in a different order" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    val reorderedDictionary = SequenceDictionary(chrM, chr1, chr2)

    dictChr1Chr2ChrM.write(inputDictionaryPath)
    reorderedDictionary.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input          = inputDictionaryPath,
      sortDictionary = sortDictionaryPath,
      output         = outputDictionaryPath,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(reorderedDictionary) shouldBe true
  }

  it should "reorder contigs even if contigs in input are an alias of contigs in sort dictionary" in {
    val aliasDictionary      = SequenceDictionary(chr1Alias, chr2Alias, chrMAlias)
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

    val reorderedDictionary = SequenceDictionary(chrM, chr1, chr2)
    val expectedDictionary  = SequenceDictionary(chrMAlias, chr1Alias, chr2Alias)


    aliasDictionary.write(inputDictionaryPath)
    reorderedDictionary.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input          = inputDictionaryPath,
      sortDictionary = sortDictionaryPath,
      output         = outputDictionaryPath,
    )
    executeFgbioTool(tool)

    SequenceDictionary(outputDictionaryPath).sameAs(expectedDictionary) shouldBe true
  }

  it should "raise Error over duplicate entries that point to the same contig in the input" in {
    val inputDictionaryPath  = makeTempFile("test.", ".1.dict")
    val sortDictionaryPath   = makeTempFile("test.", ".2.dict")
    val outputDictionaryPath = makeTempFile("test.", ".3.dict")

  val chr1IncorrectAlias = SequenceMetadata(
    name    = "chr1_2",
    length  = 100,
    aliases = Seq("chr1_3"),
  )

    val chrMIncorrectAlias = SequenceMetadata(
    name    = "MT",
    length  = 100,
    aliases = Seq( "chr1_1"),
  )

    val duplicateAliasDictionary = SequenceDictionary(chr1IncorrectAlias, chr2Alias, chrMIncorrectAlias)

    dictChr1Chr2ChrM.write(inputDictionaryPath)
    duplicateAliasDictionary.write(sortDictionaryPath)

    val tool = new SortSequenceDictionary(
      input          = inputDictionaryPath,
      sortDictionary = sortDictionaryPath,
      output         = outputDictionaryPath,
    )
    val ex = intercept[Exception] {executeFgbioTool(tool) }
    ex.getMessage should include ("Found duplicate sequence name:")
  }
}
