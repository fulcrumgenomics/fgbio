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


import com.fulcrumgenomics.FgBioDef.{FgBioEnum, PathToSequenceDictionary}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Io
import enumeratum.EnumEntry

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable.ListBuffer

/** Trait that entries in AssemblyReportColumn will extend. */
sealed trait AssemblyReportColumn extends EnumEntry {
  /** The column name in the assembly report. */
  def key: String
  /** The tag to store in the SAM header */
  def tag: String
}

/** Enum to represent columns in a NCBI assembly report. */
object AssemblyReportColumn extends FgBioEnum[AssemblyReportColumn] {
  def values: IndexedSeq[AssemblyReportColumn] = findValues

  /** Allows the column to be build from the Enum name or the actual column name. */
  override def apply(str: String): AssemblyReportColumn = {
    values.find(_.key == str).getOrElse(super.apply(str))
  }

  case object SequenceName     extends AssemblyReportColumn { val key: String = "Sequence-Name"; val tag: String = "sn" }
  case object AssignedMolecule extends AssemblyReportColumn { val key: String = "Assigned-Molecule"; val tag: String = "am"  }
  case object GenBankAccession extends AssemblyReportColumn { val key: String = "GenBank-Accn"; val tag: String = "ga"  }
  case object RefSeqAccession  extends AssemblyReportColumn { val key: String = "RefSeq-Accn"; val tag: String = "ra"  }
  case object UcscName         extends AssemblyReportColumn { val key: String = "UCSC-style-name"; val tag: String = "un"  }

  /** The missing value for a column */
  val MissingValue: String = "na"

  /** The column key for Sequence-Role */
  val SequenceRole: String = "Sequence-Role"

  /** The column key for Sequence-Length */
  val SequenceLength: String = "Sequence-Length"

  /** The sam tag to use for Sequence-Role */
  val SequenceRoleTag: String = "sr"

  values.foreach { value =>
    require(value.tag.toLowerCase == value.tag, s"Tag must be lowerase for $value: ${value.tag}")
  }
}

/** Trait that entries in SequenceRole will extend. */
sealed trait SequenceRole extends EnumEntry {
  def key: String
  def primary: Boolean
}

/** Enum to represent the various types of sequence-roles in NCBI assembly report. */
object SequenceRole extends FgBioEnum[SequenceRole] {
  def values: IndexedSeq[SequenceRole] = findValues

  /** Allows the column to be build from the Enum name or the actual sequence-role column name. */
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
    |Collates the alternate contig names from an NCBI assembly report.
    |
    |The input is to be the `*.assembly_report.txt` obtained from NCBI.
    |
    |The output will be a "sequence dictionary", which is a valid SAM file, containing the version header line and one
    |line per contig.  The primary contig name (i.e. `@SQ.SN`) is specified with `--primary` option, while alternate
    |names (i.e. aliases) are specified with the `--alternates` option.
    |
    |First, contig with the Sequence-Role "assembled-molecule" will be outputted.  Next, the remaining contigs will be
    |sorted by descending length.
  """,
  group = ClpGroups.Fasta)
class CollectAlternateContigNames
(@arg(flag='i', doc="Input NCBI assembly report file.") val input: FilePath,
 @arg(flag='o', doc="Output sequence dictionary file.") val output: PathToSequenceDictionary,
 @arg(flag='p', doc="The assembly report column for the primary contig name.") val primary: AssemblyReportColumn = AssemblyReportColumn.RefSeqAccession,
 @arg(flag='a', doc="The assembly report column(s) for the alternate contig name(s)", minElements=1) val alternates: Seq[AssemblyReportColumn],
 @arg(flag='s', doc="Only output sequences with the given sequence roles.  If none given, all sequences will be output.", minElements=0)
  val sequenceRoles: Seq[SequenceRole] = Seq.empty,
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  validate(!alternates.contains(primary), s"Primary is in alternate: $primary in " + alternates.mkString(", "))

  override def execute(): Unit = {
    val iter = Io.readLines(input).bufferBetter

    // skip over comments until we reach the header
    iter.dropWhile { line => line.startsWith("#") && !line.startsWith(s"# ${AssemblyReportColumn.SequenceName.key}") }

    // read the header
    require(iter.hasNext, s"Missing header from $input.")
    val header = iter.next().substring(2).split('\t')

    // store primary and secondaries seperately; the former will be written before the latter
    val primaries   = ListBuffer[SequenceMetadata]()
    val secondaries = ListBuffer[SequenceMetadata]()

    // Reads the input assembly report into a sequence of key/value maps
    iter.foreach { line =>
      val dict   = header.zip(line.split('\t')).toMap
      val name   = dict(this.primary.key)
      val role   = SequenceRole(dict(AssemblyReportColumn.SequenceRole))
      val alts   = this.alternates.flatMap { alt =>
        dict(alt.key) match {
          case alternate if alternate == AssemblyReportColumn.MissingValue =>
            logger.warning(s"Contig '$name' had a missing value for alternate in column '${alt.key}'")
            None
          case alternate => Some(alternate)
        }
      }
      if (sequenceRoles.nonEmpty && !this.sequenceRoles.contains(role)) {
        logger.warning(s"Skipping contig name '$name' with mismatching sequencing role: $role.")
      }
      else if (name == AssemblyReportColumn.MissingValue) {
        logger.warning(s"Skipping contig as it had a missing value for column '${this.primary.key}': $line")
      }
      else if (alts.isEmpty) {
        logger.warning(s"Skipping contig name '$name' with no alternates.")
      }
      else {
       val attributes = ListBuffer[(String, String)]()
        attributes += ((AssemblyReportColumn.SequenceRoleTag, dict(AssemblyReportColumn.SequenceRole)))
        AssemblyReportColumn.values.foreach { column =>
          attributes += ((column.tag, dict(column.key)))
        }
        val metadata = SequenceMetadata(
          name            = name,
          length          = dict(AssemblyReportColumn.SequenceLength).toInt,
          aliases         = alts,
          customAttributes = attributes.toMap
        )
        if (role.primary) {
          primaries += metadata
        }
        else {
          secondaries += metadata
        }
      }
    }
    // Sort the secondary entries based on descending sequence length.
    val infos = (primaries ++ secondaries.sortBy(-_.length))
      .zipWithIndex
      .map { case (metadata, index) => metadata.copy(index=index) }
      .toSeq

    // Write it out!
    SequenceDictionary(infos:_*).write(output)
  }
}
