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

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToSequenceDictionary, _}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{LazyLogging, StringUtil}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}

@clp(description =
  """
    |Updates the contig names in columns of a delimited data file (e.g. CSV, TSV).
    |
    |The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
    |new name will be the primary (non-alias) name in the sequence dictionary.
  """,
  group = ClpGroups.Utilities)
class UpdateDelimitedFileContigNames
(@arg(flag='i', doc="Input delimited data file.") val input: FilePath,
 @arg(flag='d', doc="The path to the sequence dictionary with contig aliases.") val dict: PathToSequenceDictionary,
 @arg(flag='c', doc="The column indices for the contig names (0-based).") val columns: Seq[Int],
 @arg(flag='T', doc="The delimiter") val delimiter: Char = '\t',
 @arg(flag='H', doc="Treat lines with this starting string as comments (always printed)") val comment: String = "#",
 @arg(flag='o', doc="Output delimited data file.") val output: FilePath,
 @arg(flag='n', doc="Output the first `N` lines as-is (always printed).") val outputFirstNumLines: Int = 0
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, dict))
  Io.assertCanWriteFile(output)

  this.columns.foreach(column => validate(column >= 0, s"Column index < 0: $column"))

  override def execute(): Unit = {
    val dict     = SequenceDictionary(this.dict)
    val in       = Io.readLines(input).bufferBetter
    val out      = Io.toWriter(output)
    val progress = ProgressLogger(logger, noun="lines", verb="written")

    // Output the requested # of lines at the start of the file
    in.take(outputFirstNumLines).foreach { line =>
      progress.record()
      out.write(line)
      out.write('\n')
    }

    // write the rest of the lines
    val maxColumn = this.columns.max
    val fields = Array.fill(maxColumn + 2)("")  // stuff the remaining in the last array value
    in.foreach { line =>
      progress.record()

      if (line.startsWith(this.comment)) {
        out.write(line)
        out.write('\n')
      }
      else {
        val numFields = StringUtil.split(line=line, delimiter=delimiter, arr=fields, concatenateRemaining=true)
        require(numFields >= maxColumn + 1, f"Too few columns '$numFields' on line ${progress.getCount}%,d")

        // update the fields and output
        columns.foreach { column =>
          dict.get(fields(column)) match {
            case Some(metadata) =>
              fields(column) = metadata.name
            case None           =>
              throw new IllegalStateException(s"Did not find contig ${fields(column)} in the given sequence dictionary.")
          }
        }

        // write it all out
        forloop(from = 0, until = numFields) { i =>
          if (0 < i) out.append(delimiter)
          out.write(fields(i))
        }
        out.write('\n')
      }
    }
    progress.logLast()

    out.close()
  }
}