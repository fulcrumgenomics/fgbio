/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.{LazyLogging, StringUtil}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}


@clp(description =
  """
    |Updates then contig names in a GFF.
    |
    |The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
    |new name will be the primary (non-alias) name in the sequence dictionary.
    |
    |Please note: the output GFF will be in the same order as the input GFF.
  """,
  group = ClpGroups.Utilities)
class UpdateGffContigNames
(@arg(flag='i', doc="Input GFF.") val input: FilePath,
 @arg(flag='d', doc="The path to the sequence dictionary with contig aliases.") val dict: PathToSequenceDictionary,
 @arg(flag='o', doc="Output GFF.") val output: FilePath,
 @arg(doc="Skip missing contigs.") val skipMissing: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, dict))
  Io.assertCanWriteFile(output)

  private val sequenceRegionKey = "##sequence-region"

  override def execute(): Unit = {
    val progress = ProgressLogger(logger, noun="features", verb="written", unit=10e5.toInt)
    val dict     = SequenceDictionary(this.dict)
    val in       = Io.readLines(input).bufferBetter
    val out      = Io.toWriter(output)
    val fields   = Array("", "", "", "")

    // writes the header lines prior to the features for a given sequence/contig
    def writeSequenceHeader(headerLines: Seq[String]): Unit = {
      headerLines.foreach { line =>
        if (line.startsWith(sequenceRegionKey)) {
          require(StringUtil.split(line, ' ', fields) == 4, s"Not enough fields: $line")
          val srcName = fields(1)
          dict.get(srcName).map { info =>
            out.append(fields(0))
              .append(' ').append(info.name)
              .append(' ').append(fields(2))
              .append(' ').append(fields(3))
              .append('\n')
          }
        }
        else {
          out.append(line).append('\n')
        }
      }
    }

    // writes the features for the given sequence/contig (src)
    def writeSequenceFeatures(in: BetterBufferedIterator[String],
                              srcName: String,
                              targetName: String): Unit = {
      in.takeWhile(_.startsWith(srcName)).foreach { line =>
        require(StringUtil.split(line, '\t', fields) == fields.length, s"Not enough fields: $line")
        require(fields(0) == srcName)
        val position = fields(3).toInt
        out.append(targetName)
        out.append(line.substring(srcName.length))
        out.append('\n')
        progress.record(s"$srcName => $targetName", position)
      }
    }

    // write out the header
    in.takeWhile(line => line.startsWith("#") && !line.startsWith(sequenceRegionKey))
      .foreach(out.append(_).append('\n'))

    // Process a contig at a time
    while (in.nonEmpty) {

      // read in the header or the first entry to determine the sequence/src name
      val header         = in.takeWhile(_.startsWith("#")).toSeq
      val sequenceRegion = header.find(_.startsWith(sequenceRegionKey))
      val srcName        = sequenceRegion.map { line =>
        require(StringUtil.split(line, ' ', fields) == fields.length, s"Not enough fields: $line")
        fields(1)
      }.orElse {
        in.headOption.map { line =>
          require(StringUtil.split(line, '\t', fields) == fields.length, s"Not enough fields: $line")
          fields.head
        }
      }

      // write the header and the features
      srcName match {
        case None      => header.foreach(line => out.append(line).append('\n')) // just write the header
        case Some(src) =>
          dict.get(src) match {
            case None =>
              val message = s"Did not find contig $src in the list of source names."
              if (skipMissing) {
                logger.warning(message)
                in.dropWhile(_.startsWith(src))
              }
              else throw new IllegalStateException(message)
            case Some(info) =>
              writeSequenceHeader(header)
              writeSequenceFeatures(in, src, info.name)
          }
      }
    }
    progress.logLast()

    out.close()
  }
}