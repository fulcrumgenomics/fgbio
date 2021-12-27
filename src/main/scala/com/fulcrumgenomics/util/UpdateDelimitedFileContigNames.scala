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
import com.fulcrumgenomics.util.Sorter.Codec
import enumeratum.EnumEntry

import scala.collection.immutable.IndexedSeq

@clp(description =
  """
    |Updates the contig names in columns of a delimited data file (e.g. CSV, TSV).
    |
    |The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
    |new name will be the primary (non-alias) name in the sequence dictionary.  Use `--skip-missing` to ignore lines
    |where a contig name could not be updated (i.e. missing from the sequence dictionary).
  """,
  group = ClpGroups.Utilities)
class UpdateDelimitedFileContigNames
(@arg(flag='i', doc="Input delimited data file.") val input: FilePath,
 @arg(flag='d', doc="The path to the sequence dictionary with contig aliases.") val dict: PathToSequenceDictionary,
 @arg(flag='c', doc="The column indices for the contig names (0-based).") val columns: Seq[Int],
 @arg(flag='T', doc="The delimiter") val delimiter: Char = '\t',
 @arg(flag='H', doc="Treat lines with this starting string as comments (always printed)") val comment: String = "#",
 @arg(flag='o', doc="Output delimited data file.") val output: FilePath,
 @arg(flag='n', doc="Output the first `N` lines as-is (always printed).") val outputFirstNumLines: Int = 0,
 @arg(doc="Skip lines where a contig name could not be updated (i.e. missing from the sequence dictionary).") skipMissing: Boolean = false,
 @arg(flag='s', doc="Sort the output based on the following order.") sortOrder: SortOrder = SortOrder.Unsorted,
 @arg(doc="The column index for the contig (0-based) for sorting. Use the first column if not given.") val contig: Option[Int] = None,
 @arg(doc="The column index for the genomic position (0-based) for sorting by coordinate.") val position: Option[Int] = None,
 @arg(doc="The maximum number of objects to store in memory") val maxObjectsInRam: Int = 1e6.toInt
) extends FgBioTool with LazyLogging {

  import UpdateDelimitedFileContigNames._

  Io.assertReadable(input)
  Io.assertReadable(Seq(input, dict))
  Io.assertCanWriteFile(output)

  this.columns.foreach(column => validate(column >= 0, s"Column index < 0: $column"))

  if (sortOrder == SortOrder.ByCoordinate) {
    validate(position.nonEmpty, "--position is required when sorting by coordinate")
  }

  override def execute(): Unit = {
    val dict              = SequenceDictionary(this.dict)
    val in                = Io.readLines(input).bufferBetter
    val out               = Io.toWriter(output)
    val progress          = ProgressLogger(logger, noun="lines", verb="written")
    val maxColumn         = this.columns.max
    val contigColumnIndex = this.contig.getOrElse(this.columns.head)
    val sorter            = {
      if (SortOrder.Unsorted == this.sortOrder) None
      else {
        val codec: LineInfoCodec = new LineInfoCodec(
          delimiter           = this.delimiter,
          maxColumn           = maxColumn,
          contigColumnIndex   = this.contig.getOrElse(this.columns.head),
          positionColumnIndex = position.getOrElse(-1),
          toContigIndex       = (name: String) => dict(name).index
        )
        Some(new Sorter[LineInfo, LineInfoKey](maxObjectsInRam=maxObjectsInRam, codec=codec, keyfunc=_.key))
      }
    }
    val delimiterString   = this.delimiter.toString

    // Output the requested # of lines at the start of the file
    in.take(outputFirstNumLines).foreach { line =>
      progress.record()
      out.write(line)
      out.write('\n')
    }

    // write the rest of the lines
    val fields: Array[String] = Array.fill(maxColumn + 2)("")  // stuff the remaining in the last array value
    in.foreach { line =>
      progress.record()

      if (line.startsWith(this.comment)) {
        out.write(line)
        out.write('\n')
      }
      else {
        val numFields = StringUtil.split(line=line, delimiter=delimiter, arr=fields, concatenateRemaining=true)
        require(numFields >= maxColumn + 1, f"Too few columns '$numFields' on line ${progress.getCount}%,d")

        // Find the first column that can't be updated, but update the fields with those that can be updated along the way
        val missingContig: Option[Int] = columns.find { column =>
          dict.get(fields(column)) match {
            case Some(metadata) => fields(column) = metadata.name; false
            case None           => true
          }
        }

        // Log or error if there was a column that couldn't be updated, otherwise, write it out
        missingContig match {
          case Some(column) if skipMissing =>
            logger.warning(s"Skipping line, could not update contig with name '${fields(column)}' in column $column on line ${progress.getCount}%,d")
          case Some(column) =>
            throw new IllegalStateException(s"Did not find contig ${fields(column)} in the given sequence dictionary.")
          case None =>
            // sorting
            sorter match {
              case None => // Write out the fields without sorting
                forloop(from = 0, until = numFields) { i =>
                  if (0 < i) out.append(delimiter)
                  out.write(fields(i))
                }
                out.write('\n')
              case Some(_sorter) => // add then to the sorter, to be written later
                _sorter += LineInfo(
                  contig   = dict(fields(contigColumnIndex)).index,
                  position = if (this.sortOrder == SortOrder.ByContigOnly) 0 else position.map(i => fields(i).toInt).getOrElse(0),
                  lineno   = progress.getCount.toInt,
                  line     = fields.take(numFields).mkString(delimiterString)
                )
            }
        }
      }
    }
    progress.logLast()

    // Sort and write the output
    sorter.foreach { _sorter =>
      logger.info("Sorting the output")
      _sorter.foreach { data =>
        out.write(data.line)
        out.write('\n')
      }
      _sorter.close()
    }

    out.close()
  }
}

/** Trait that entries in SortOrder will extend. */
sealed trait SortOrder extends EnumEntry

/** Enum to represent the sort order. */
object SortOrder extends FgBioEnum[SortOrder] {
  case object Unsorted extends SortOrder
  case object ByContigOnly extends SortOrder
  case object ByCoordinate extends SortOrder
  def values: IndexedSeq[SortOrder] = findValues
}

object UpdateDelimitedFileContigNames {

  /** The key used to compare across lines in the input delimited file. */
  private case class LineInfoKey(contig: Int, position: Int, lineno: Int) extends Ordered[LineInfoKey] {
    /** Compares by contig, position, then line number, ascending. */
    override def compare(that: LineInfoKey): Int = {
      var retval = this.contig - that.contig
      if (retval == 0) retval = this.position - that.position
      if (retval == 0) retval = this.lineno - that.lineno
      retval
    }
  }

  private object LineInfo {
    /** Builds a [[LineInfo]]. */
    def apply(contig: Int, position: Int, lineno: Int, line: String): LineInfo = {
      LineInfo(LineInfoKey(contig=contig, position=position, lineno=lineno), line=line)
    }
  }

  /** Container class for a line in the input delimited file, and a key used to sort the lines. */
  private case class LineInfo(key: LineInfoKey, line: String)

  /** [[Codec]] for [[LineInfo]] */
  private class LineInfoCodec(val delimiter: Char,
                              val maxColumn: Int,
                              val contigColumnIndex: Int,
                              val positionColumnIndex: Int = -1,
                              val toContigIndex: String => Int
                             ) extends Codec[LineInfo] {
    // add one for the line number, and stuff any remaining items in the last array value
    private val fields: Array[String] = Array.fill(maxColumn + 3)("")
    private val delimiterString: String = delimiter.toString

    def encode(a: LineInfo): Array[Byte] = f"${a.key.lineno}$delimiter${a.line}".getBytes
    override def decode(bs: Array[Byte], start: Int, length: Int): LineInfo = {
      val line = new String(bs, start, length)
      StringUtil.split(line=line, delimiter=delimiter, arr=fields, concatenateRemaining=true)
      val lineno   = fields(0).toInt
      val contig   = fields(contigColumnIndex + 1)
      val position = if (positionColumnIndex >= 0) fields(positionColumnIndex + 1).toInt else 0
      LineInfo(
        contig   = toContigIndex(contig),
        position = position,
        lineno   = lineno,
        line     = fields.drop(1).mkString(delimiterString)
      )
    }
  }
}