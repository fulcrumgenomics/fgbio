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
                      comment: String = "#"): Unit = {
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
}
