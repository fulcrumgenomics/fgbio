/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics
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

import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.Interval
import htsjdk.tribble.annotation.Strand
import org.scalactic.{Equality, Explicitly}
import org.scalatest.OptionValues._

import java.io.ByteArrayInputStream
import scala.io.Source

/** Unit tests for [[IntervalSource]]. */
class IntervalSourceTest extends UnitSpec with Explicitly {

  /** Equality helper for intervals that fully compares them by contig, start, end, name, and strand. */
  private val equalityByLocatableNameAndStrand: Equality[Interval] = (a: Interval, b: Any) => b match {
    case expected: Interval => a.contigsMatch(expected) &&
      a.getStart        == expected.getStart &&
      a.getEnd          == expected.getEnd &&
      a.getStrand       == expected.getStrand &&
      Option(a.getName) == Option(expected.getName)
    case _ => false
  }

  /** A sequence dictionary for unit testing. */
  private val Dict: SequenceDictionary = SequenceDictionary(
    SequenceMetadata(name = "chr1", length = 50000000),
    SequenceMetadata(name = "chr2", length = 50000)
  )

  /** Contents of a BED file without a header for testing. */
  private val BedWithoutHeader: String =
    """chr19 49302000 49302300 -1.0
      |chr19 49302300 49302600 -0.75
      |chr19 49302600 49302900 -0.50
      |chr19 49302900 49303200 -0.25
      |chr19 49303200 49303500 0.0
      |chr19 49303500 49303800 0.25
      |chr19 49303800 49304100 0.50
      |chr19 49304100 49304400 0.75
      |chr20 10       100      1.00
    """.stripMargin.trim

  /** Deserialized version of [[BedWithoutHeader]] for testing. */
  private val BedWithoutHeaderExpected: Seq[Interval] = {
    Seq(
      new Interval("chr19", 49302001, 49302300, false, "-1.0"),
      new Interval("chr19", 49302301, 49302600, false, "-0.75"),
      new Interval("chr19", 49302601, 49302900, false, "-0.50"),
      new Interval("chr19", 49302901, 49303200, false, "-0.25"),
      new Interval("chr19", 49303201, 49303500, false, "0.0"),
      new Interval("chr19", 49303501, 49303800, false, "0.25"),
      new Interval("chr19", 49303801, 49304100, false, "0.50"),
      new Interval("chr19", 49304101, 49304400, false, "0.75"),
      new Interval("chr20", 11,       100,      false, "1.00"),
    )
  }

  /** Contents of a BED file (actually a bedGraph) with a complex header for testing. */
  private val BedWithHeader: String =
    """browser position chr19:49302001-49304701
      |browser hide all
      |browser pack refGene encodeRegions
      |browser full altGraph
      |#	300 base wide bar graph, autoScale is on by default == graphing
      |#	limits will dynamically change to always show full range of data
      |#	in viewing window, priority = 20 positions this as the second graph
      |#	Note, zero-relative, half-open coordinate system in use for bedGraph format
      |track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20
      |chr19 49302000 49302300 -1.0
      |chr19 49302300 49302600 -0.75
      |chr19 49302600 49302900 -0.50
      |chr19 49302900 49303200 -0.25
      |chr19 49303200 49303500 0.0
      |chr19 49303500 49303800 0.25
      |chr19 49303800 49304100 0.50
      |chr19 49304100 49304400 0.75
      |chr19 49304400 49304700 1.00
    """.stripMargin.trim

  /** Deserialized version of [[BedWithHeader]] for testing. */
  private val BedWithHeaderExpected: Seq[Interval] = {
    Seq(
      new Interval("chr19", 49302001, 49302300, false, "-1.0"),
      new Interval("chr19", 49302301, 49302600, false, "-0.75"),
      new Interval("chr19", 49302601, 49302900, false, "-0.50"),
      new Interval("chr19", 49302901, 49303200, false, "-0.25"),
      new Interval("chr19", 49303201, 49303500, false, "0.0"),
      new Interval("chr19", 49303501, 49303800, false, "0.25"),
      new Interval("chr19", 49303801, 49304100, false, "0.50"),
      new Interval("chr19", 49304101, 49304400, false, "0.75"),
      new Interval("chr19", 49304401, 49304700, false, "1.00"),
    )
  }

  /** Contents of an Interval List header for testing. */
  private val IntervalListData: String =
    "@HD\tVN:1.0\n" + "@SQ\tSN:chr1\tLN:50000000\n" + "@SQ\tSN:chr2\tLN:50000\n" +
      "chr1\t49302000\t49302300\t+\tname1\n" +
      "chr1\t49302300\t49302600\t+\tname2\n" +
      "chr1\t49302600\t49302900\t+\tname3\n" +
      "chr1\t49302900\t49303200\t+\tname4\n" +
      "chr1\t49303200\t49303500\t+\tname5\n" +
      "chr1\t49303500\t49303800\t+\tname6\n" +
      "chr1\t49303800\t49304100\t+\tname7\n" +
      "chr1\t49304100\t49304400\t+\tname8\n" +
      "chr2\t10\t100\t-\tname9\n"

  /** Deserialized version of [[IntervalListData]] for testing. */
  private val IntervalListExpected: Seq[Interval] = Seq(
    new Interval("chr1", 49302000, 49302300, false, "name1"),
    new Interval("chr1", 49302300, 49302600, false, "name2"),
    new Interval("chr1", 49302600, 49302900, false, "name3"),
    new Interval("chr1", 49302900, 49303200, false, "name4"),
    new Interval("chr1", 49303200, 49303500, false, "name5"),
    new Interval("chr1", 49303500, 49303800, false, "name6"),
    new Interval("chr1", 49303800, 49304100, false, "name7"),
    new Interval("chr1", 49304100, 49304400, false, "name8"),
    new Interval("chr2", 10,       100,      true,  "name9"),
  )

  "IntervalSource" should "read a single simple BED data record into an interval" in {
    val actual   = IntervalSource("chr1 100\n".linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read no intervals at all from an empty iterator" in {
    val actual = IntervalSource(Iterator.empty, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe None
    actual.toList shouldBe empty
  }

  it should "return BED intervals from a list of line records" in {
    val actual = IntervalSource(BedWithoutHeader.linesIterator, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs BedWithoutHeaderExpected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "return BED intervals, skipping a header, from a list of line records" in {
    val actual = IntervalSource(BedWithHeader.linesIterator, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs BedWithHeaderExpected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read a single simple Interval List data record into an interval" in {
    val data   = "@HD\tVN:1.0\n" + "@SQ\tSN:chr1\tLN:50000000\n" + "chr1\t100\t100\t+\t.\n"
    val dict   = SequenceDictionary(SequenceMetadata("chr1", length = 50000000))
    val header = new SAMFileHeader()
    header.setAttribute("VN", "1.0")
    header.setSequenceDictionary(dict.asSam)

    val source   = IntervalSource(data.linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 100, 100, false, "."))

    source.dict.value shouldBe dict
    source.header.value shouldBe header
    (source.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read multiple intervals from an interval list" in {
    val source = IntervalSource(IntervalListData.linesIterator, dict = None)
    val header = new SAMFileHeader()
    header.setAttribute("VN", "1.0")
    header.setSequenceDictionary(Dict.asSam)

    source.dict.value shouldBe Dict
    source.header.value shouldBe header
    (source.toList should contain theSameElementsInOrderAs IntervalListExpected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read an interval that has a name intentionally set from a BED source" in {
    val actual   = IntervalSource("chr1 100 101 interval-name\n".linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 101, 101, false, "interval-name"))
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read an interval that has the strand set to positive from a BED source" in {
    val actual   = IntervalSource("chr1 100 101 interval-name 500 +\n".linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 101, 101, false, "interval-name"))
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read an interval that has the strand set to the unknown value from a BED source" in {
    val actual   = IntervalSource("chr1 100 101 interval-name 500 .\n".linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 101, 101, false, "interval-name"))
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "read an interval that has the strand set to negative from a BED source" in {
    val actual   = IntervalSource("chr1 100 101 interval-name 500 -\n".linesIterator, dict = None)
    val expected = Seq(new Interval("chr1", 101, 101, true, "interval-name"))
    actual.dict shouldBe None
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "assert the provided sequence dictionary matches the found sequence dictionary for interval list input" in {
    val invalid = SequenceDictionary(SequenceMetadata("chrX", length =2))
    val caught  = intercept[IllegalArgumentException] { IntervalSource(IntervalListData.linesIterator, dict = Some(invalid)) }
    caught.getMessage should include ("Provided sequence dictionary does not match the input's dict header!")
  }

  it should "pass through the sequence dictionary when supplied for BED input" in {
    val dict   = SequenceDictionary(
      SequenceMetadata(name = "chr19", length = 50000000),
      SequenceMetadata(name = "chr20", length = 50000)
    )
    val actual = IntervalSource(BedWithHeader.linesIterator, dict = Some(dict))
    actual.dict.value shouldBe dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs BedWithHeaderExpected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "assert the records match the provided sequence dictionary for BED input" in {
    val dict   = SequenceDictionary(
      SequenceMetadata(name = "chr19", length = 50000000),
      SequenceMetadata(name = "chr20", length = 50000)
    )
    val caught = intercept[NoSuchElementException] { IntervalSource("chr1 100".linesIterator, dict = Some(dict)).toList }
    caught.getMessage should include ("Contig does not exist within dictionary for locatable")
    caught.getMessage should include ("Failed on line number: 1")
  }

  "IntervalSource.apply" should "allow sourcing intervals from an iterable of string data" in {
    val actual   = IntervalSource(Seq("chr1 100"), dict = Some(Dict))
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict.value shouldBe Dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "allow sourcing interval from an input stream of string data" in {
    val stream   = new ByteArrayInputStream("chr1 100\n".getBytes(java.nio.charset.StandardCharsets.UTF_8.name))
    val actual   = IntervalSource(stream, dict = Some(Dict))
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict.value shouldBe Dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "allow sourcing intervals from a source of string data" in {
    val source   = Source.fromString("chr1 100\n")
    val actual   = IntervalSource(source, dict = Some(Dict))
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict.value shouldBe Dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
    actual.close()
  }

  it should "allow sourcing intervals from a file of string data" in {
    val path = Io.makeTempFile(getClass.getSimpleName, ".bed")
    Io.writeLines(path, Seq("chr1 100"))
    val actual   = IntervalSource(path.toFile, dict = Some(Dict))
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict.value shouldBe Dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }

  it should "allow sourcing intervals from a path of string data" in {
    val path = Io.makeTempFile(getClass.getSimpleName, ".bed")
    Io.writeLines(path, Seq("chr1 100"))
    val actual   = IntervalSource(path, dict = Some(Dict))
    val expected = Seq(new Interval("chr1", 101, 101))
    actual.dict.value shouldBe Dict
    actual.header shouldBe None
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableNameAndStrand)
  }
}
