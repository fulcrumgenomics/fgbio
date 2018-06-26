/**
  * Copyright (c) 2016, Fulcrum Genomics LLC
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *
  * 1. Redistributions of source code must retain the above copyright notice,
  * this list of conditions and the following disclaimer.
  *
  * 2. Redistributions in binary form must reproduce the above copyright notice,
  * this list of conditions and the following disclaimer in the documentation
  * and/or other materials provided with the distribution.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */

package com.fulcrumgenomics.illumina

import java.nio.file.Path

import scala.collection.Set
import scala.io.Source

object SampleSheet {

  /** Required header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
  val SampleId: String       = "Sample_ID"
  val SampleName: String     = "Sample_Name"
  val LibraryId: String      = "Library_ID"
  /** Optional header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
  val SampleProject: String  = "Sample_Project"
  val Description: String    = "Description"
  val Lane: String           = "Lane"
  val I7Bases: String        = "Index"
  val I5Bases: String        = "Index2"

  /** The header names that are standard fields in upper case. */
  val HeaderNames = Set(SampleId, SampleName, LibraryId, SampleProject, Description, Lane, I7Bases, I5Bases).map(_.toUpperCase)

  /** Attempts to clean special characters from a sample id or name as Illumina's software does. */
  def cleanMiseqSampleId(id: String): String = {
    id.map {
      case '#' => ' '
      case '_' | '+' | ' ' => '-'
      case c => c
    }
  }

  /**
    * Parses an Illumina Experiment Manager Sample Sheet (typically a MiSeq).
    *
    * If library id is not found, it will be set to the sample identifier.
    *
    * If the lane number is given, returns only the samples that have that lane specified.
    *
    * NB: the lookup in the columns in the Sample Sheet is case insensitive.
    */
  def apply(sampleSheet: Path, lane: Option[Int] = None): SampleSheet = {
    apply(lines=Source.fromFile(sampleSheet.toFile).getLines(), lane=lane)
  }

  /**
    * Parses an Illumina Experiment Manager Sample Sheet (typically a MiSeq).
    *
    * If library id is not found, it will be set to the sample identifier.
    *
    * If the lane number is given, returns only the samples that have that lane specified.
    *
    * NB: the lookup in the columns in the Sample Sheet is case insensitive.
    */
  def apply(lines: Iterator[String], lane: Option[Int]): SampleSheet = {
    val samples = parseSampleData(lines=lines, lane=lane).zipWithIndex.map {
      case (sampleDatum, sampleOrdinal) => makeSample(sampleOrdinal+1, sampleDatum)
    }
    new SampleSheet(samples=samples)
  }

  /** Reads in the sample data from an Illumina Experiment Manager sample sheet.  For each sample, a map of column name
    * to value is returned. If the lane number is given, returns only the samples that have that lane specified. */
  protected def parseSampleData(lines: Iterator[String], lane: Option[Int] = None): List[Map[String, String]] = {
    // ignore data pre-"[Data]"
    val (preData, postData) = lines.toList.span(!_.startsWith("[Data]"))

    // get the header (keys) (skip "[Data]")
    val header = (postData.drop(1).headOption match {
      case Some(line) => line.split(SplitRegex, -1) // NB: include trailing empty strings
      case None => throw new IllegalArgumentException("Could not find the header for sample data.")
    }).map(_.toUpperCase.trim)

    // get the rows (values), so skip "[Data]" and the header
    val lineNumber = preData.size + 2 // 0-based
    val sampleDatums = postData.drop(2).zipWithIndex.map { case (line, rowNumber) =>
        val values = line.split(SplitRegex, -1)
        // check we have the correct # of columns
        if (values.size != header.length) {
          throw new IllegalArgumentException(s"Found a line with a mismatching number of columns on line #$lineNumber: " + line)
        }
        header.zip(values.map(_.trim)).toMap
      }
    lane match {
      case Some(l: Int) => sampleDatums.filter { SampleSheet.getIntField(_, Lane).exists(_ == l) }
      case None         => sampleDatums
    }
  }

  /** Creates a sample from the given row data. */
  protected def makeSample(sampleOrdinal: Int, sampleDatum: Map[String, String]): Sample = {
    val sampleName: String          = SampleSheet.getStringField(sampleDatum, SampleName)    getOrElse (throw new IllegalArgumentException(s"Missing: $SampleName"))
    val sampleId: String            = SampleSheet.getStringField(sampleDatum, SampleId)      getOrElse (throw new IllegalArgumentException(s"Missing: $SampleId"))
    val libraryId: String           = SampleSheet.getStringField(sampleDatum, LibraryId)     getOrElse sampleId
    val project: Option[String]     = SampleSheet.getStringField(sampleDatum, SampleProject)
    val description: Option[String] = SampleSheet.getStringField(sampleDatum, Description)
    val extendedAttributes          = sampleDatum.filterNot { case (columnName, value) =>  HeaderNames.contains(columnName) }
    new Sample(
      sampleOrdinal        = sampleOrdinal,
      sampleId             = sampleId,
      sampleName           = sampleName,
      libraryId            = libraryId,
      project              = project,
      description          = description,
      lane                 = SampleSheet.getIntField(sampleDatum, Lane),
      i7IndexBases         = SampleSheet.getStringField(sampleDatum, I7Bases),
      i5IndexBases         = SampleSheet.getStringField(sampleDatum, I5Bases),
      extendedAttributes   = extendedAttributes
    )
  }

  /** Gets the trimmed value of the given key (`name`) from the map.  If empty or not found, returns None */
  protected def getStringField(sampleDatum: Map[String, String], name: String): Option[String] = {
    sampleDatum.get(name.toUpperCase) match {
      case None => None
      case Some(str) if str.isEmpty => None
      case Some(str) => Some(str.trim)
    }
  }

  protected def getIntField(sampleDatum: Map[String, String], name: String): Option[Int] = {
    getStringField(sampleDatum, name).map(_.toInt)
  }

  protected[illumina] val SplitRegex = ",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)"

  /** Validate that the sample ids are unique, and that the combination of sample name and library id are unique. */
  private def validateSamplesFromTheSameLane(samples: Seq[Sample]): Unit = {
    // All samples should either be from the same lane or have no lane
    require(samples.flatMap(_.lane).distinct.lengthCompare(1) <= 0)
    /** We include the lane if present in the sample key for debugging purposes*/
    def sampleKey(s: Sample): String = s"${s.sampleName} ${s.libraryId}" + s.lane.map(" " + _).getOrElse("")

    // Validate some properties for a sample sheet
    require(samples.map(_.sampleId).toSet.size == samples.length,
      "Sample identifiers were not unique: " + samples.groupBy(_.sampleId).filter(_._2.lengthCompare(1) > 0).keys.mkString(", "))
    require(samples.map(sampleKey).toSet.size == samples.length,
      "Sample name and library identifier combinations were not unique: " +
        samples.groupBy(sampleKey).filter(_._2.lengthCompare(1) > 0).keys.mkString(", "))
  }
}

/**
  * Stores information about samples from an Illumina Experiment Manager Sample Sheet (typically a MiSeq).  The
  * samples may also include derived information.
  *
  * Optional fields include support for specifying the expected sample barcode for each sample.  The sample barcode
  * can be present as a sub-sequence (or sub-sequences) in the i7 or i5 read.  If additional bases are found in the i7
  * or i5 read, such as molecular barcodes, the should be included as Ns.  It is up to the developer to obtain the
  * correct read structure elsewhere to infer which bases are sample barcode and which bases are not (ex. molecular
  * identifiers).
  *
  * The sample identifier should be unique across samples within the sample sheet, and the combination of sample name
  * and library identifier should also be unique across samples within the sample sheet.
  */
class SampleSheet(samples: Seq[Sample]) extends Iterable[Sample] {

  // Validate samples
  samples.flatMap(_.lane).toSet[Int].toList.foreach { lane =>
    SampleSheet.validateSamplesFromTheSameLane(samples.filter(_.lane.contains(lane)))
  }
  SampleSheet.validateSamplesFromTheSameLane(samples.filter(_.lane.isEmpty))

  def iterator: Iterator[Sample] = this.samples.iterator

  override def size: Int = this.samples.size

  def get(index: Int): Sample = this.samples(index)
}
