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

package com.fulcrumgenomics.personal.nhomer

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import htsjdk.samtools.util.CloserUtil
import htsjdk.samtools.{SAMFileWriter, SAMFileWriterFactory, SamReader, SamReaderFactory}

@clp(
  description = "Splits an optional tag in a SAM or BAM into multiple optional tags.",
  group = ClpGroups.Personal
)
class SplitTag
( @arg(doc = "Input SAM or BAM.") val input: PathToBam,
  @arg(doc = "Output SAM or BAM.") val output: PathToBam,
  @arg(doc = "Tag to split.") val tagToSplit: String,
  @arg(doc = "Tag(s) to output.  There should be one per produced token.", minElements = 1) val tagsToOutput: List[String],
  @arg(doc = "The delimiter used to split the string.") val delimiter: String = "-"
) extends FgBioTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  if (tagToSplit.length != 2) throw new ValidationException(s"The tag to split must be of length two (was ${tagToSplit.length}).")
  tagsToOutput.foreach { tagToOutput =>
    if (tagToOutput.length != 2) throw new ValidationException(s"The tag to output '$tagToOutput' must be of length two (was ${tagToOutput.length}).")
  }

  override def execute(): Unit = {
    val reader: SamReader = SamReaderFactory.makeDefault.open(input.toFile)
    val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader, true, output.toFile)
    reader.foreach { record =>
      val value: String = record.getStringAttribute(tagToSplit)
      if (value == null) fail(String.format("Record '%s' was missing the tag '%s'", record.getReadName, tagToSplit))
      val tokens: Array[String] = value.split(delimiter)
      if (tokens.length != tagsToOutput.size) fail(s"Record '${record.getReadName}' did not have '${tagsToOutput.size}' tokens")
      tagsToOutput.zip(tokens).foreach { case (tagToOutput, token) => record.setAttribute(tagToOutput, token) }
      writer.addAlignment(record)
    }
    CloserUtil.close(reader)
    writer.close()
  }
}
