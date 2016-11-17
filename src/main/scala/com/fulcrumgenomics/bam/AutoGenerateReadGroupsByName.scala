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
    |Adds read groups to a BAM file for a single sample by parsing the read names.
    |
    |Will add one or more read groups by parsing the read names.  The read names should be of the form:
    |  <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos>
    |
    |Each unique combination of <instrument>:<run number>:<flowcell ID>:<lane> will be its own read group. The ID of the
    |read group will be an integer and the platform unit will be <flowcell-id>.<lane>.
    |
    |The input is assumed to contain reads for one sample and library.  Therefore, the sample and library must be given
    |and will be applied to all read groups.  Read groups will be replaced if present.
    |
    |Two passes will be performed on the input: first to gather all the read groups, and second to write the output BAM
    |file.
  """
)
class AutoGenerateReadGroupsByName
(@arg(flag="i", doc="Input SAM or BAM file") val input: PathToBam,
 @arg(flag="o", doc="Output SAM or BAM file") val output: PathToBam,
 @arg(flag="s", doc="The sample to insert into the read groups") val sample: String,
 @arg(flag="l", doc="The library to insert into the read groups") val library: String,
 @arg(doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
 @arg(doc="Predicted median insert size, to insert into the read groups") val predictedInsertSize: Option[Integer] = None,
 @arg(doc="Program group to insert into the read groups") val programGroup: Option[String] = None,
 @arg(doc="Platform model to insert into the groups (free-form text providing further details of the platform/technology used)") val platformModel: Option[String] = None,
 @arg(doc="Description inserted into the read groups") val description: Option[String] = None,
 @arg(doc="Date the run was produced (ISO 8601: YYYY-MM-DD), to insert into the read groups") val runDate: Option[Iso8601Date] = None,
 @arg(doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
 @arg(doc="The sort order for the output sam/bam file.") val sortOrder: Option[SortOrder] = None
) extends FgBioTool with LazyLogging {

  import scala.collection.JavaConversions.asScalaIterator

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

        if (!readGroups.contains(runInfo)) {
          val readGroup = new SAMReadGroupRecord(nextId.incrementAndGet().toString)

          readGroup.setSample(sample)
          readGroup.setLibrary(library)
          readGroup.setPlatformUnit(runInfo.platformUnit)
          platformModel.foreach(readGroup.setPlatformModel)
          readGroup.setPlatform("ILLUMINA")
          sequencingCenter.foreach(readGroup.setSequencingCenter)
          predictedInsertSize.foreach(readGroup.setPredictedMedianInsertSize)
          programGroup.foreach(readGroup.setProgramGroup)
          description.foreach(readGroup.setDescription)
          runDate.foreach(readGroup.setRunDate)

          readGroups(runInfo) = readGroup
        }
        progress.record(record)
      }
      in.safelyClose()
    }

    // Write them all out
    {
      val progress  = new ProgressLogger(logger, verb = "written", unit = 1e6.toInt)
      val in        = SamReaderFactory.make().open(input.toFile)
      val inHeader  = in.getFileHeader
      val outHeader = inHeader.clone()
      sortOrder.foreach(outHeader.setSortOrder)
      outHeader.setReadGroups(java.util.Arrays.asList(readGroups.values.toSeq:_*))
      comments.foreach(outHeader.addComment)
      val out = new SAMFileWriterFactory().makeWriter(outHeader, sortOrder.forall(_ == inHeader.getSortOrder), output.toFile, null)

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
}

/** Stores parsed run-level information about a read. */
private[bam] case class RunInfo(instrument: String, runNumber: String, flowcellId: String, lane: Int)
{
  /** The platform unit: <flowcell.lane> */
  val platformUnit: String = s"$flowcellId.$lane"
}

private[bam] object RunInfo {
  /** Parse a read name to get run-level information. */
  def apply(name: String): RunInfo = {
    val tokens = name.split(':').toSeq
    if (tokens.length != 7) throw new IllegalArgumentException(s"Expected seven colon-delimited fields for read with name: $name")
    new RunInfo(
      instrument = tokens(0),
      runNumber  = tokens(1),
      flowcellId = tokens(2),
      lane       = tokens(3).toInt
    )
  }
}
