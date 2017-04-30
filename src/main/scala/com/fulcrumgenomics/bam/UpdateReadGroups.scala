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
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import htsjdk.samtools._

import scala.collection.JavaConverters._

/**
  * Updates one or more read groups.
  *
  * @author Nils Homer
  */
@clp(description =
  """
      |Updates one or more read groups and their identifiers.
      |
      |This tool will replace each read group with a new read group, including a new read group identifier.  If the read
      |group identifier is not to be changed, it is recommended to use `samtools reheader` or Picard's
      |`ReplaceSamHeader` instead as in this case only the header needs modification. If all read groups are to be
      |assigned to one read group, it is recommended to use Picard's `AddOrReplaceReadGroups`.  Nonetheless, if the read
      |group identifier also needs to be changed, use this tool.
      |
      |Each read group in the input file will be mapped to one and only one new read group identifier, unless
      |`--ignore-missing-read-groups` is set.  A SAM header file should be given with the new read groups and the ID field
      |foreach read group containing the new read group identifier.  An additional attribute (`FR`) should be provided
      |that gives the original read group identifier (`ID`) to which this new read group corresponds.
      |
      |If `--keep-read-group-attributes` is true, then any read group attribute not replaced will be kept in the new read
      |group.  Otherwise, only the attributes in the provided SAM header file will be used.
    """,
  group = ClpGroups.SamOrBam)
class UpdateReadGroups
(@arg(flag='i', doc = "Input BAM file.") val input: PathToBam,
 @arg(flag='o', doc = "Output BAM file.") val output: PathToBam,
 @arg(flag='r', doc = "A SAM header file with the replacement read groups (see detailed usage).") val readGroupsFile: PathToBam,
 @arg(flag='k', doc = "Keep all read group attributes that are not replaced.") val keepReadGroupAttributes: Boolean = false,
 @arg(flag='g', doc = "Keep all read groups not found in the replacement header, otherwise throw an error.") val ignoreMissingReadGroups: Boolean = false
) extends FgBioTool with LazyLogging {

  private val fromReadGroupKey = "FR"

  Io.assertReadable(input)
  Io.assertReadable(readGroupsFile)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in              = SamSource(input)
    val toReadGroupsMap = getReadGroupMap(readGroupsFile, in.header)
    val outHeader = in.header.clone()
    outHeader.setReadGroups(toReadGroupsMap.values.toList.sortBy(_.getId).asJava)
    val out = SamWriter(output, outHeader)

    // main loop
    in.foreach { record =>
      updateRecord(record, outHeader, toReadGroupsMap)
      out += record
    }

    in.safelyClose()
    out.close()
  }

  private def updateRecord(record: SamRecord,
                           header: SAMFileHeader,
                           readGroupMap: Map[String, SAMReadGroupRecord]): SamRecord = {
    val fromReadGroup = record.readGroup match {
      case null => fail(s"Record '${record.name}' has no read group")
      case rg => rg
    }
    readGroupMap.get(fromReadGroup.getId) match {
      case None => unreachable(s"Read group id '${record.readGroup.getId}' from record '${record.name}' not found in the SAM header.")
      case Some(readGroup: SAMReadGroupRecord) =>
        record(SAMTag.RG.name()) = readGroup.getId
        record.header = header
    }
    record
  }

  private def getReadGroupMap(readGroupsFile: FilePath, samFileHeader: SAMFileHeader): Map[String, SAMReadGroupRecord] = {
    val fromReadGroups = samFileHeader.getReadGroups.map { rg => rg.getId -> rg }.toMap

    // get the header
    val reader = SamSource(readGroupsFile)
    val h = reader.header
    reader.safelyClose()

    // map from the "FR" attribute value to the new read group
    val readGroupMap = h.getReadGroups.map { rg =>
      rg.getAttribute(fromReadGroupKey) match {
        case null => fail(s"'$fromReadGroupKey' attribute not found in read group: " + rg.toString)
        case fromId =>
          // remove the "FR" attribute
          rg.setAttribute(fromReadGroupKey, null)
          // verify the `id` is found in `fromReadGroups`
          if (!fromReadGroups.contains(fromId)) fail(s"Read group id '$fromId' not found in input BAM file.")
          // keep the old attributes that are not being overwritten
          if (keepReadGroupAttributes) {
            // find all read group attributes that are not present in the new read group and add those attributes to the new read group
            fromReadGroups(fromId).getAttributes.filter(attr => rg.getAttribute(attr.getKey) == null).foreach { attr => rg.setAttribute(attr.getKey, attr.getValue) }
          }
          fromId -> rg
      }
    }.toMap

    // get the read groups in the original file that are not being updated
    val missingReadGroups = samFileHeader.getReadGroups.filter(rg => !readGroupMap.contains(rg.getId)).map(rg => rg.getId -> rg).toMap
    
    if (ignoreMissingReadGroups) {
      readGroupMap ++ missingReadGroups
    }
    else if (missingReadGroups.nonEmpty) fail("Read groups found that are not being replaced: " + missingReadGroups.keys.mkString(", "))
    else {
      readGroupMap
    }
  }
}
