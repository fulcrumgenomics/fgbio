/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam

import java.util.concurrent.atomic.AtomicInteger

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.util.Iso8601Date

import scala.collection.mutable

@clp(group = ClpGroups.SamOrBam, description =
  """
    |Adds read groups to a BAM file for a single sample by parsing the read name.
    |
    |Will add one or read groups by parsing the read names.  Assumes the read names are of the form:
    |  <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos> <read>:<is filtered>:<control number>:<barcode sequence>
    |
    |Each unique combination of <instrument>:<run number>:<flowcell ID>:<lane> will be its own read group. The ID of the
    |read group will be an integer, the platform unit will be <flowcell-id>.<lane>, and the platform model if not given
    |will be <instrument.run-number>.
    |
    |The input is assumed to contain reads for one sample and library.  Therefore, the sample and library must be given
    |and will be applied to all read groups.
    |
    |Two passes will be performed on the input: first to gather all the read groups, and second to write the output BAM
    |file.
  """
)
class AddReadGroupsByName
(@arg(doc="Input SAM or BAM file") val input: PathToBam,
 @arg(doc="Output SAM or BAM file") val output: PathToBam,
 @arg(doc="The sample to insert into the read group header") val sample: String,
 @arg(doc="The library to insert into the read group header") val library: String,
 @arg(doc="The platform type (e.g. illumina, solid) to insert into the read group header") val platform: Option[String] = Some("ILLUMINA"),
 @arg(doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
 @arg(doc="Predicted median insert size, to insert into the read group header") val predictedInsertSize: Option[Integer] = None,
 @arg(doc="Program group to insert into the read group header.") val programGroup: Option[String] = None,
 @arg(doc="Platform model to insert into the group header (free-form text providing further details of the platform/technology used)") val platformModel: Option[String] = None,
 @arg(doc="Description inserted into the read group header") val description: Option[String] = None,
 @arg(doc="Date the run was produced, to insert into the read group header") val runDate: Option[Iso8601Date] = None,
 @arg(doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
 @arg(doc="The sort order for the output sam/bam file.") val sortOrder: SortOrder = SortOrder.queryname
) extends FgBioTool with LazyLogging {

  import scala.collection.JavaConversions.{asScalaIterator, seqAsJavaList}

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  private val nextId = new AtomicInteger(0)
  private val readGroups = new mutable.HashMap[RunInfo, SAMReadGroupRecord]()

  override def execute(): Unit = {

    // Gather all the read groups
    {
      val progress = new ProgressLogger(logger, verb = "read", unit = 5e6.toInt)
      val in = SamReaderFactory.make().open(input.toFile)

      in.iterator().foreach { record =>
        val runInfo       = RunInfo(name=record.getReadName)
        val platformModel = this.platformModel.getOrElse(runInfo.platformModel)

        if (!readGroups.contains(runInfo)) {
          val readGroup = new SAMReadGroupRecord(nextId.incrementAndGet().toString)

          readGroup.setSample(sample)
          readGroup.setLibrary(library)
          readGroup.setPlatformUnit(runInfo.platformUnit)
          readGroup.setPlatformModel(platformModel)
          platform.foreach(readGroup.setPlatform)
          sequencingCenter.foreach(readGroup.setSequencingCenter)
          predictedInsertSize.foreach(readGroup.setPredictedMedianInsertSize)
          description.foreach(readGroup.setDescription)
          runDate.foreach(readGroup.setRunDate)

          readGroups(runInfo) = readGroup
        }
        else {
          // Check that the read information matches the read group
          val readGroup = readGroups(runInfo)

          // Used in the lookup of the read group, should never happen!
          if (readGroup.getPlatformUnit != runInfo.platformUnit) unreachable(s"PU mismatch: s'${readGroup.getPlatformModel} != ${runInfo.platformUnit}")
          if (readGroup.getPlatformModel != platformModel) unreachable(s"PM mismatch: s'${readGroup.getPlatformModel} != $platformModel")

          // These should be set across all reads, so should never happen!
          if (platform.exists(_ != readGroup.getPlatform)) unreachable(s"PL mismatch: s'${readGroup.getPlatform} != ${platform.get}")
          if (sequencingCenter.exists(_ != readGroup.getSequencingCenter)) unreachable(s"CN mismatch: s'${readGroup.getSequencingCenter} != sequencingCenter.get}")
          if (predictedInsertSize.exists(_ != readGroup.getPredictedMedianInsertSize)) unreachable(s"PI mismatch: s'${readGroup.getPredictedMedianInsertSize} != ${predictedInsertSize.get}")
          if (description.exists(_ != readGroup.getDescription)) unreachable(s"DS mismatch: s'${readGroup.getDescription} != ${description.get}")
          if (runDate.exists(_ != readGroup.getRunDate)) unreachable(s"DT mismatch: s'${readGroup.getRunDate} != ${runDate.get}")
        }
        progress.record(record)
      }
      in.safelyClose()
    }

    // Write them all out
    {
      val progress = new ProgressLogger(logger, verb = "written", unit = 5e6.toInt)
      val in = SamReaderFactory.make().open(input.toFile)
      val header = in.getFileHeader
      val outHeader = header.clone()
      outHeader.setSortOrder(sortOrder)
      outHeader.setReadGroups(readGroups.values.toSeq)
      comments.foreach(outHeader.addComment)
      val out = new SAMFileWriterFactory().makeWriter(outHeader, sortOrder == header.getSortOrder, output.toFile, null)

      in.iterator().foreach { record =>
        val runInfo = RunInfo(name = record.getReadName)
        val readGroupId = readGroups(runInfo).getReadGroupId
        record.setAttribute(SAMTag.RG.name(), readGroupId)
        out.addAlignment(record)
        progress.record(record)
      }

      in.safelyClose()
      out.close()
    }
  }

  /** Generates a read group identifier for the read: <sample>.<platform-unit>. */
  private def toReadGroupId(runInfo: RunInfo): String = s"$sample.${runInfo.platformUnit}"
}

/** Stores parsed run-level information about a read. */
private[bam] case class RunInfo(instrument: String,
                                runNumber: String,
                                flowcellId: String,
                                lane: Int)
{
  /** The platform unit: <flowcell.lane> */
  val platformUnit: String = s"$flowcellId.$lane"

  /** The platform mode: <instrument.run-number> */
  def platformModel: String = s"$instrument.$runNumber"
}

private[bam] object RunInfo {
  /** Parse a read name to get run-level information. */
  def apply(name: String): RunInfo = {
    val tokens = name.split(":").toSeq
    if (tokens.length != 7) throw new IllegalArgumentException(s"Expected seven colon-delimited fields for read with name: $name")
    new RunInfo(
      instrument = tokens(0),
      runNumber  = tokens(1),
      flowcellId = tokens(2),
      lane       = tokens(3).toInt
    )
  }
}
