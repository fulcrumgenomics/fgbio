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

import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec
import htsjdk.tribble.bed.{BEDFeature, FullBEDFeature, SimpleBEDFeature}
import org.scalactic.{Equality, Explicitly}
import org.scalatest.OptionValues._

import java.io.ByteArrayInputStream
import scala.io.Source

/** Unit tests for [[BedSource]]. */
class BedSourceTest extends UnitSpec with Explicitly {

  /** A sequence dictionary for unit testing. */
  private val Dict: SequenceDictionary = SequenceDictionary(SequenceMetadata("chr1", length = 10000))

  /** Convenience method for building a Full BED Feature from a contig, start, end, and name. */
  private def bed(contig: String, start: Int, end: Int, name: String): FullBEDFeature = {
    val feature = new FullBEDFeature(contig, start, end)
    feature.setName(name)
    feature
  }

  /** Equality helper for BED features that only compares contig, start, end, and name. */
  private val equalityByLocatableAndName: Equality[BEDFeature] = (a: BEDFeature, b: Any) => b match {
    case expected: BEDFeature => a.contigsMatch(expected) &&
      a.getStart == expected.getStart &&
      a.getEnd   == expected.getEnd &&
      a.getName  == expected.getName
    case _ => false
  }

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
  private val BedWithoutHeaderExpected: Seq[FullBEDFeature] = {
    Seq(
      bed("chr19", 49302001, 49302300, "-1.0"),
      bed("chr19", 49302301, 49302600, "-0.75"),
      bed("chr19", 49302601, 49302900, "-0.50"),
      bed("chr19", 49302901, 49303200, "-0.25"),
      bed("chr19", 49303201, 49303500, "0.0"),
      bed("chr19", 49303501, 49303800, "0.25"),
      bed("chr19", 49303801, 49304100, "0.50"),
      bed("chr19", 49304101, 49304400, "0.75"),
      bed("chr20", 11,       100,      "1.00"),
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
  private val BedWithHeaderExpected: Seq[FullBEDFeature] = {
    Seq(
      bed("chr19", 49302001, 49302300, "-1.0"),
      bed("chr19", 49302301, 49302600, "-0.75"),
      bed("chr19", 49302601, 49302900, "-0.50"),
      bed("chr19", 49302901, 49303200, "-0.25"),
      bed("chr19", 49303201, 49303500, "0.0"),
      bed("chr19", 49303501, 49303800, "0.25"),
      bed("chr19", 49303801, 49304100, "0.50"),
      bed("chr19", 49304101, 49304400, "0.75"),
      bed("chr19", 49304401, 49304700, "1.00"),
    )
  }

  "BedSource" should "read a single simple BED data record" in {
    val actual   = BedSource("chr1 100\n".linesIterator, dict = None)
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict shouldBe None
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
  }

  it should "read no BED features from an empty iterator" in {
    val actual   = BedSource(Iterator.empty, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe empty
    actual.toList shouldBe empty
  }

  it should "return BED features from a list of line records" in {
    val actual = BedSource(BedWithoutHeader.linesIterator, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs BedWithoutHeaderExpected) (decided by equalityByLocatableAndName)
  }

  it should "return BED features from a list of line records and keep a copy of the header" in {
    val actual = BedSource(BedWithHeader.linesIterator, dict = None)
    actual.dict shouldBe None
    actual.header shouldBe Seq(
      "browser position chr19:49302001-49304701",
      "browser hide all",
      "browser pack refGene encodeRegions",
      "browser full altGraph",
      "#\t300 base wide bar graph, autoScale is on by default == graphing",
      "#\tlimits will dynamically change to always show full range of data",
      "#\tin viewing window, priority = 20 positions this as the second graph",
      "#\tNote, zero-relative, half-open coordinate system in use for bedGraph format",
      "track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20",
    )
    (actual.toList should contain theSameElementsInOrderAs BedWithHeaderExpected) (decided by equalityByLocatableAndName)
  }

  it should "raise an exception when the underlying BED codec cannot read the data line" in {
    an[IllegalStateException] shouldBe thrownBy { BedSource("\n".linesIterator, dict = None).toList }
    an[IllegalStateException] shouldBe thrownBy { BedSource("    \n".linesIterator, dict = None).toList }
    an[IllegalStateException] shouldBe thrownBy { BedSource("chr1\n".linesIterator, dict = None).toList }
  }

  it should "have a sequence dictionary if one is passed to the constructor" in {
    val source = BedSource(Iterator.empty, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    source.toList shouldBe empty
  }

  it should "raise a NoSuchElementException if a BED record has a contig not in the optional sequence dictionary" in {
    val source = BedSource("chr2 100\n".linesIterator, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    val caught = intercept[NoSuchElementException] { source.toList }
    caught.getMessage should include ("Contig does not exist within dictionary for locatable")
    caught.getMessage should include ("Failed on line number: 1")
  }

  it should "raise a IllegalArgumentException if a BED record has a start value less than 1 (1-based)" in {
    val source = BedSource("chr1 -1 50\n".linesIterator, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    val caught = intercept[IllegalArgumentException] { source.toList }
    caught.getMessage should include ("Start is less than 1 for locatable")
    caught.getMessage should include ("Failed on line number: 1")
  }

  it should "raise a IllegalArgumentException if a BED record has an end value greater than the contig length" in {
    val source = BedSource("chr1 0 10001\n".linesIterator, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    val caught = intercept[IllegalArgumentException] { source.toList }
    caught.getMessage should include ("End is beyond the reference contig length for locatable")
    caught.getMessage should include ("Failed on line number: 1")
  }

  it should "raise a IllegalArgumentException if a BED record has a start value greater than an end value" in {
    val source = BedSource("chr1 100 50\n".linesIterator, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    val caught = intercept[IllegalArgumentException] { source.toList }
    caught.getMessage should include ("Start is greater than end for locatable")
    caught.getMessage should include ("Failed on line number: 1")
  }

  it should "know which line of input triggered a validation exception" in {
    val source = BedSource("chr1 100\nchr2 100".linesIterator, dict = Some(Dict))
    source.dict.value shouldBe Dict
    source.header shouldBe empty
    val caught = intercept[NoSuchElementException] { source.toList }
    caught.getMessage should include ("Contig does not exist within dictionary for locatable")
    caught.getMessage should include ("Failed on line number: 2")
  }

  "BedSource.apply" should "allow sourcing BED features from an iterable of string data" in {
    val actual   = BedSource(Seq("chr1 100"), dict = Some(Dict))
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict.value shouldBe Dict
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
  }

  it should "allow sourcing BED features from an input stream of string data" in {
    val stream   = new ByteArrayInputStream("chr1 100\n".getBytes(java.nio.charset.StandardCharsets.UTF_8.name))
    val actual   = BedSource(stream, dict = Some(Dict))
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict.value shouldBe Dict
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
  }

  it should "allow sourcing BED features from a source of string data" in {
    val source   = Source.fromString("chr1 100\n")
    val actual   = BedSource(source, dict = Some(Dict))
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict.value shouldBe Dict
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
    actual.close()
  }

  it should "allow sourcing BED features from a file of string data" in {
    val path = Io.makeTempFile(getClass.getSimpleName, ".bed")
    Io.writeLines(path, Seq("chr1 100"))
    val actual   = BedSource(path.toFile, dict = Some(Dict))
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict.value shouldBe Dict
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
  }

  it should "allow sourcing BED features from a path of string data" in {
    val path = Io.makeTempFile(getClass.getSimpleName, ".bed")
    Io.writeLines(path, Seq("chr1 100"))
    val actual   = BedSource(path, dict = Some(Dict))
    val expected = Seq(new SimpleBEDFeature(101, 101, "chr1"))
    actual.dict.value shouldBe Dict
    actual.header shouldBe empty
    (actual.toList should contain theSameElementsInOrderAs expected) (decided by equalityByLocatableAndName)
  }
}
