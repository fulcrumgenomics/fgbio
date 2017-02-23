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

package com.fulcrumgenomics.basecalling

import com.fulcrumgenomics.FgBioDef.{DirPath, FilePath, PathToBam, unreachable}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.illumina.{Sample, SampleSheet}
import com.fulcrumgenomics.util.Io
import dagr.commons.io.PathUtil
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}

import scala.collection.mutable.ListBuffer


@clp(group=ClpGroups.Basecalling, description=
  """Extracts sample and library information from an sample sheet for a given lane.
    |
    |The sample sheet should be an Illumina Experiment Manager sample sheet. The tool writes two files to the output
    |directory: a barcode parameter file and a library parameter file.
    |
    |The barcode parameter file is used by Picard's ExtractIlluminaBarcodes and CollectIlluminaBasecallingMetrics to
    |determine how to match sample barcodes to each read.  The parameter file will be written to the output directory
    |with name "barcode_params.<lane>.txt".
    |
    |The library parameter file is used by Picard's IlluminaBasecallsToSam to demultiplex samples and name the output
    |BAM file path for each sample output BAM file.  The parameter file will be written to the output directory with name
    |"library_params.<lane>.txt".  The path to each sample's BAM file will be specified in the library parameter
    |file.  Each BAM file will have path "<output>/<sample-name>.<barcode-sequence>.<lane>.bam".
  """)
class ExtractBasecallingParamsForPicard
(
  @arg(flag="i", doc="The input sample sheet.") val input: FilePath,
  @arg(flag="o", doc="The output folder to where per-lane parameter files should be written.") val output: DirPath,
  @arg(flag="b", doc="Optional output folder to where per-lane BAM files should be written, otherwise the output directory will be used.") val bamOutput: Option[DirPath] = None,
  @arg(flag="l", doc="The lane(s) (1-based) for which to write per-lane parameter files.") val lanes: Seq[Int]
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertWritableDirectory(output)
  lanes.foreach { lane => validate(0 < lane, "--lane must be greater than zero") }

  override def execute(): Unit = {
    val params = BasecallingParams.from(sampleSheet=SampleSheet(sampleSheet=input), lanes=lanes, output=output, bamOutput=bamOutput)
    params.foreach { param =>
      logger.info(s"Barcode params for lane ${param.lane} written to: ${param.barcodeFile}")
      logger.info(s"Library params for lane ${param.lane} written to: ${param.libraryParamsFile}")
    }
  }
}


/** A small class to store where the basecalling parameter and BAM files are located. */
case class BasecallingParams(lane: Int, barcodeFile: FilePath, libraryParamsFile: FilePath, bams: Seq[PathToBam])

/** Contains methods useful for extractin basecalling paramter information. */
object BasecallingParams {

  /** Extracts sample and library information from an sample sheet for a given lane.  The resulting files are used
    * to run Picard's ExtractIlluminaBarcodes, CollectIlluminaBasecallingMetrics, and IlluminaBasecallsToSam. See
    * [[ExtractBasecallingParamsForPicard]] for more information.
    *
    * @param sampleSheet the Illumina Experiment Manager sample sheet.
    * @param lanes the lane number(s) (1-based).
    * @param output the output directory where the barcode, library param, and BAM files should be written.
    * @param bamOutput the optional output directory to where per-lane BAM files will be written, otherwise use output.
    */
  def from(sampleSheet: SampleSheet, lanes: Seq[Int], output: DirPath, bamOutput: Option[DirPath] = None): Seq[BasecallingParams] = lanes.map { lane =>
    require(0 < lane, s"lane '$lane' was not greater than zero (1-based).")

    val barcodeLines      = new ListBuffer[String]()
    val libraryParamLines = new ListBuffer[String]()
    val samples           = sampleSheet.filter(_.lane.forall(_ == lane)).toSeq
    val barcodeFile       = barcodeFileFrom(output=output, lane=lane)
    val libraryParamsFile = libraryParamsFileFrom(output=output, lane=lane)

    // Dual vs single indexed
    val dualIndexed = samples.exists(_.i5IndexBases.nonEmpty)
    if (dualIndexed) require(samples.forall(_.i5IndexBases.nonEmpty), "Not all samples were dual-indexed")
    else require(samples.forall(_.i5IndexBases.isEmpty), "Not all samples were single-indexed")
    require(samples.forall(_.i7IndexBases.nonEmpty), "Missing i7 bases for some samples")

    // Check for description on the sample
    val includeDescription    = samples.exists(_.description.nonEmpty)
    val standardLibraryParams = Seq("OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME") ++ (if (includeDescription) Seq("DS") else Seq.empty)

    // Headers
    if (dualIndexed) {
      barcodeLines += Seq("barcode_sequence_1", "barcode_sequence_2", "barcode_name", "library_name").mkString("\t")
      libraryParamLines += (Seq("BARCODE_1", "BARCODE_2") ++ standardLibraryParams).mkString("\t")
    }
    else {
      barcodeLines += Seq("barcode_sequence_1", "barcode_name", "library_name").mkString("\t")
      libraryParamLines += (Seq("BARCODE_1") ++ standardLibraryParams).mkString("\t")
    }

    // Samples
    val bams = ListBuffer[PathToBam]()
    samples.foreach { sample =>
      // Get the library id
      val libraryId = sample.libraryId

      // Get the barcodes
      val barcodes = Seq(sample.i7IndexBases, sample.i5IndexBases).flatten
      require(barcodes.nonEmpty, s"No sample barcodes found for sample: ${sample.sampleName}")

      // Get the description
      val description = if (includeDescription) Seq(sample.description).flatten else Seq.empty

      // Add the barcode line
      barcodeLines += (barcodes ++ Seq(barcodes.mkString, libraryId)).mkString("\t")

      // Add the library param line
      val bam            = bamFileFrom(output=bamOutput.getOrElse(output), sample=sample, lane=lane)
      val sampleAlias    = sample.sampleName
      libraryParamLines += (barcodes ++ Seq(bam, sampleAlias, libraryId) ++ description).mkString("\t")

      bams.append(bam)
    }

    // Ensure that the BAM file names are unique
    require(bams.toSet.size == bams.size, "BAM file names collide: ." + bams.groupBy(_.toString).filter(_._2.length > 1).keys.mkString(", "))

    // Add the unmatched (only for library params!)
    val description = if (includeDescription) Seq("unmatched") else Seq.empty
    val barcodes    = if (dualIndexed) Seq("N", "N") else Seq("N")
    libraryParamLines += (barcodes ++ Seq(bamOutput.getOrElse(output).resolve(s"unmatched.$lane.bam"), "unmatched", "unmatched") ++ description).mkString("\t")

    // Output
    Io.writeLines(path=barcodeFile, lines=barcodeLines)
    Io.writeLines(path=libraryParamsFile, lines=libraryParamLines)

    BasecallingParams(lane=lane, barcodeFile=barcodeFile, libraryParamsFile=libraryParamsFile, bams=bams)
  }

  /** The lane-specific file containing information necessary to run Picard's ExtractIlluminaBarcodes and
    * CollectIlluminaBasecallingMetrics.
    * @param output the base output directory.
    * @param lane the lane number.
    * @return the path to the barcode file, named "barcode_data.<lane>.txt".
    */
  def barcodeFileFrom(output: DirPath, lane: Int): FilePath = {
    require(0 < lane)
    output.resolve(s"barcode_params.$lane.txt")
  }

  /** The lane-specific file containing information necessary to run Picard's IlluminaBasecallsToSam.
    * @param output the base output directory.
    * @param lane the lane number.
    * @return the path to the library params file, named "library_params.<lane>.txt".
    */
  def libraryParamsFileFrom(output: DirPath, lane: Int): FilePath = {
    require(0 < lane)
    output.resolve(s"library_params.$lane.txt")
  }

  /** Produces the path to the BAM file for the given sample and lane.
    *
    * The BAM file name will <sample-name>.<barcode-seq>.<lane>.bam.  The resulting file name will be sanitized to
    * remove problematic characters.
    *
    * @param output the base output directory.
    * @param sample the sample.
    * @param lane the lane number.
    * @return the path to the BAM file for the given sample and lane.
    */
  def bamFileFrom(output: DirPath, sample: Sample, lane: Int): PathToBam = {
    require(0 < lane)
    val barcodeSeq = Seq(sample.i7IndexBases, sample.i5IndexBases).flatten.mkString
    output.resolve(PathUtil.sanitizeFileName(s"${sample.sampleName}.$barcodeSeq.$lane.bam"))
  }
}
