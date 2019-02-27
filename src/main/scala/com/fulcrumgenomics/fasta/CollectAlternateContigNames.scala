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

package com.fulcrumgenomics.fasta


import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Io
import enumeratum.EnumEntry

import scala.collection.immutable.IndexedSeq

sealed trait AssemblyReportColumn extends EnumEntry {
  def key: String
}

object AssemblyReportColumn extends FgBioEnum[AssemblyReportColumn] {
  def values: IndexedSeq[AssemblyReportColumn] = findValues

  override def apply(str: String): AssemblyReportColumn = {
    values.find(_.key == str).getOrElse(super.apply(str))
  }

  case object SequenceName     extends AssemblyReportColumn { val key: String = "Sequence-Name" }
  case object AssignedMolecule extends AssemblyReportColumn { val key: String = "Assigned-Molecule" }
  case object GenBankAccession extends AssemblyReportColumn { val key: String = "GenBank-Accn" }
  case object RefSeqAccession  extends AssemblyReportColumn { val key: String = "RefSeq-Accn" }
  case object UcscName         extends AssemblyReportColumn { val key: String = "UCSC-style-name" }

  val MissingValue: String   = "na"
  val SequenceRole: String   = "Sequence-Role"
  val SequenceLength: String = "Sequence-Length"
}


sealed trait SequenceRole extends EnumEntry {
  def key: String
  def primary: Boolean
}

object SequenceRole extends FgBioEnum[SequenceRole] {
  def values: IndexedSeq[SequenceRole] = findValues

  override def apply(str: String): SequenceRole = {
    values.find(_.key == str).getOrElse(super.apply(str))
  }

  case object AltScaffold         extends SequenceRole { val key: String = "alt-scaffold"; val primary: Boolean = false }
  case object AssembledMolecule   extends SequenceRole { val key: String = "assembled-molecule"; val primary: Boolean = true }
  case object FixPatch            extends SequenceRole { val key: String = "fix-patch"; val primary: Boolean = false }
  case object NovelPatch          extends SequenceRole { val key: String = "novel-patch"; val primary: Boolean = false }
  case object UnlocalizedScaffold extends SequenceRole { val key: String = "unlocalized-scaffold"; val primary: Boolean = false }
  case object UnplacedScaffold    extends SequenceRole { val key: String = "unplaced-scaffold"; val primary: Boolean = false }
}


@clp(description =
  """
    |Gets the alternate contig names from an NCBI assembly report.
    |
    |The input is be the `*.assembly_report.txt` obtained from NCBI.
    |
    |The output will have the first column is the original name and the second column is an alternative name.  If there
    |is more than one alternate name, each alternate name will be on a separate line.
    |
    |First, sequences with the Sequence-Role "assembled-molecule" will be outputted.  Next, the remaining sequences will
    |be sorted by length, with smallest output first.
  """,
  group = ClpGroups.Fasta)
class CollectAlternateContigNames
( @arg(flag='i', doc="Input NCBI assembly report file.") val input: FilePath,
  @arg(flag='o', doc="Output file.")val output: FilePath,
  @arg(flag='p', doc="The assembly report column for the primary contig name.") val primary: AssemblyReportColumn = AssemblyReportColumn.RefSeqAccession,
  @arg(flag='a', doc="The assembly report column(s) for the alternate contig name(s)", minElements=1) val alternates: Seq[AssemblyReportColumn],
  @arg(flag='s', doc="Only output sequences with the given sequence roles.  If none given, all sequences will be output.", minElements=0)
  val sequenceRoles: Seq[SequenceRole] = Seq.empty,
  @arg(doc="Output all alternates on a single line.") val singleLine: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  validate(!alternates.contains(primary), s"Primary is in alternate: $primary in " + alternates.mkString(", "))

  override def execute(): Unit = {
    // Go through each sequence
    val lines = readAssemblyReport().toIterator.flatMap { dict =>
      val primary = dict(this.primary.key)

      val altNames = this.alternates.flatMap { alt =>
        dict(alt.key) match {
          case alternate if alternate == AssemblyReportColumn.MissingValue =>
            logger.warning(s"Missing value for alternate ${alt.key}: $primary")
            None
          case alternate => Some(alternate)
        }
      }
      if (primary == AssemblyReportColumn.MissingValue) {
        logger.warning(s"Primary contig name '${this.primary.key}' had a missing value: $primary")
        Seq.empty
      }
      else if (altNames.isEmpty) {
        logger.warning(s"No alternates found for: $primary")
        Seq.empty
      }
      else if (singleLine) {
        Seq(primary + "\t" + altNames.mkString("\t"))
      }
      else {
        altNames.map { altName => s"$primary\t$altName" }
      }
    }

    Io.writeLines(output, lines)
  }

  private def readAssemblyReport(): Seq[Map[String, String]] = {
    val iter = Io.readLines(input).bufferBetter

    // skip over comments until we reach the header
    iter.dropWhile { line => line.startsWith("#") && !line.startsWith(s"# ${AssemblyReportColumn.SequenceName.key}") }

    // read the header
    require(iter.hasNext, s"Missing header from $input.")
    val header = iter.next().split('\t')

    // Reads the input assembly report into a sequence of key/value maps
    val maps = iter.flatMap { line =>
      val fields = line.split('\t')
      val d      = header.zip(fields).toMap
      if (sequenceRoles.isEmpty || sequenceRoles.exists(_.key == d(AssemblyReportColumn.SequenceRole))) {
        Some(d)
      } else {
        None
      }
    }.toSeq

    // Re-order the entries based on primary vs non primary sequence roles, then based on increasing sequence length
    val primary   = maps.filter { d => SequenceRole(d(AssemblyReportColumn.SequenceRole)).primary }
    val secondary = maps.filter { d => !SequenceRole(d(AssemblyReportColumn.SequenceRole)).primary}.sortBy { d => d(AssemblyReportColumn.SequenceLength)}
    primary ++ secondary
  }
}
