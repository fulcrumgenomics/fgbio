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

import com.fulcrumgenomics.FgBioDef.PathToSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec

class UpdateDelimitedFileContigNamesTest extends UnitSpec {

  private def toSequenceMetadata(name: String, alias: String*): SequenceMetadata = {
    SequenceMetadata(name=name, length=0, aliases=alias)
  }

  private val dict: SequenceDictionary = {
    SequenceDictionary(
      toSequenceMetadata(name="1-new", "1-old"),
      toSequenceMetadata(name="2-new", "2-old"),
      toSequenceMetadata(name="3-new", "3-old")
    )
  }

  private val pathToSequenceDictionary: PathToSequenceDictionary = {
    val path = makeTempFile("test.", "in.dict")
    dict.write(path)
    path
  }

  private def runTest(delimiter: Char,
                      columns: Seq[Int],
                      actual: Seq[String],
                      expected: Seq[String],
                      outputFirstNumLines: Int = 0,
                      comment: String = "#",
                      skipMissing: Boolean = false,
                      sortOrder: SortOrder = SortOrder.Unsorted,
                      contig: Option[Int] = None,
                      position: Option[Int] = None,
                      maxObjectsInRam: Int = 1e6.toInt): Unit = {
    val input    = makeTempFile("test.", "in.txt")
    val output   = makeTempFile("test.", "out.txt")

    Io.writeLines(input, actual)

    val tool = new UpdateDelimitedFileContigNames(
      input               = input,
      dict                = pathToSequenceDictionary,
      columns             = columns,
      delimiter           = delimiter,
      comment             = comment,
      output              = output,
      outputFirstNumLines = outputFirstNumLines,
      skipMissing         = skipMissing,
      sortOrder           = sortOrder,
      contig              = contig,
      position            = position,
      maxObjectsInRam     = maxObjectsInRam
    )

    executeFgbioTool(tool)

    Io.readLines(output).toSeq should contain theSameElementsInOrderAs expected
  }

  "UpdateDelimitedFileContigNames" should "update a delimited file" in {
    Seq(',', '\t', ':').foreach { delimiter =>
      runTest(
        delimiter = delimiter,
        columns   = Seq(2),
        actual    = Seq(Seq("na", "na", "1-old", "na", "na").mkString(delimiter.toString)),
        expected  = Seq(Seq("na", "na", "1-new", "na", "na").mkString(delimiter.toString)),
      )
    }
  }

  it should "update multiple columns in a delimited file" in {
    Seq(',', '\t', ':').foreach { delimiter =>
      runTest(
        delimiter = delimiter,
        columns   = Seq(2, 4),
        actual    = Seq(Seq("3-old", "na", "1-old", "na", "2-old").mkString(delimiter.toString)),
        expected  = Seq(Seq("3-old", "na", "1-new", "na", "2-new").mkString(delimiter.toString))
      )
    }
  }

  it should "skip the first n-lines" in {
    Seq(0, 1, 2).foreach { outputFirstNumLines =>
      val actual = Seq(
        "1-old,na",
        "2-old,na",
        "3-old,na"
       )
      val allModified = Seq(
        "1-new,na",
        "2-new,na",
        "3-new,na"
      )
      val expected = actual.take(outputFirstNumLines) ++ allModified.drop(outputFirstNumLines)
      runTest(
        delimiter           = ',',
        columns             = Seq(0),
        actual              = actual,
        expected            = expected,
        outputFirstNumLines = outputFirstNumLines
      )
    }
  }

  it should "skip the comment lines" in {
    Seq(0, 1, 2).foreach { outputFirstNumLines =>
      val actual = Seq(
        "$1-old,na",
        "2-old,na",
        "$2-old,na",
        "3-old,na"
      )
      val expected = Seq(
        "$1-old,na",
        "2-new,na",
        "$2-old,na",
        "3-new,na"
      )
      runTest(
        delimiter = ',',
        columns   = Seq(0),
        actual    = actual,
        expected  = expected,
        comment   = "$"
      )
    }
  }

  it should "fail if contig names are not found in the dict" in {
    val exception = intercept[Exception] {
      runTest(
        delimiter = '\t',
        columns   = Seq(2),
        actual    = Seq(Seq("na", "na", "4-old", "na", "na").mkString("\t")),
        expected  = Seq.empty,
      )
    }
    exception.getMessage.contains("Did not find contig") shouldBe true
  }

  it should "skip lines if not all contig names can be updated with --skip-missing" in {
    val actual = Seq(
      "1-old,4-old", // skipped, second can't be updated
      "2-old,3-old", // updated
      "4-old,5-old", // skipped, both can't be updated
      "5-old,1-old", // skipped, first can't be updated
    )
    val expected = Seq(
      "2-new,3-new"
    )
    runTest(
      delimiter   = ',',
      columns     = Seq(0, 1),
      actual      = actual,
      expected    = expected,
      skipMissing = true
    )
  }

  it should "sort by contig name only" in {
    val actual = Seq(
      "1-old,0",
      "3-old,4",
      "3-old,3",
      "2-old,2"
    )
    val expected = Seq(
      "1-new,0",
      "2-new,2",
      "3-new,4",
      "3-new,3"
    )
    // --contig not specified
    runTest(
      delimiter   = ',',
      columns     = Seq(0),
      actual      = actual,
      expected    = expected,
      sortOrder   = SortOrder.ByContigOnly
    )
    // --contig specified
    runTest(
      delimiter   = ',',
      columns     = Seq(0),
      actual      = actual,
      expected    = expected,
      sortOrder   = SortOrder.ByContigOnly,
      contig      = Some(0)
    )
  }

  it should "sort by coordinate" in {
    val actual = Seq(
      "1-old,0,2-old,1",
      "3-old,4,2-old,2",
      "3-old,1,2-old,3",
      "2-old,2,2-old,4",
      "2-old,2,2-old,5"
    )
    // sorts by col0, col1, then line number (col3)s
    val expectedCol0 = Seq(
      "1-new,0,2-new,1",
      "2-new,2,2-new,4",
      "2-new,2,2-new,5",
      "3-new,1,2-new,3",
      "3-new,4,2-new,2"
    )
    // col2 has the same value for all, so sorting by position in col1
    val expectedCol2 = Seq(
      "1-new,0,2-new,1",
      "3-new,1,2-new,3",
      "2-new,2,2-new,4",
      "2-new,2,2-new,5",
      "3-new,4,2-new,2"
    )
    // --contig is not specified, so defaults to column 0
    runTest(
      delimiter   = ',',
      columns     = Seq(0, 2),
      actual      = actual,
      expected    = expectedCol0,
      sortOrder   = SortOrder.ByCoordinate,
      position    = Some(1)
    )
    // --contig is specified as column 0
    runTest(
      delimiter   = ',',
      columns     = Seq(0, 2),
      actual      = actual,
      expected    = expectedCol0,
      sortOrder   = SortOrder.ByCoordinate,
      contig      = Some(0),
      position    = Some(1)
    )
    // --contig is not specified, so defaults to column 0
    runTest(
      delimiter   = ',',
      columns     = Seq(0, 2),
      actual      = actual,
      expected    = expectedCol2,
      sortOrder   = SortOrder.ByCoordinate,
      contig      = Some(2),
      position    = Some(1)
    )
    // tests that when all values are the same, we sort by line number
    val clones = Seq.range(0, 100).map(_ => "1-new,0,2-new,1")
    runTest(
      delimiter   = ',',
      columns     = Seq(0, 2),
      actual      = clones,
      expected    = clones,
      sortOrder   = SortOrder.ByCoordinate,
      contig      = Some(2),
      position    = Some(1)
    )
  }

  it should "spill to disk when sorting" in {
    val actual   = Seq.range(0, 1000).map { i => f"1-old,0,2-old,$i" }
    val expected = Seq.range(0, 1000).map { i => f"1-new,0,2-new,$i" }
    runTest(
      delimiter       = ',',
      columns         = Seq(0, 2),
      actual          = actual,
      expected        = expected,
      sortOrder       = SortOrder.ByCoordinate,
      position        = Some(1),
      maxObjectsInRam = 10
    )
  }
}
