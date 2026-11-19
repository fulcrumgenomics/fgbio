/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics
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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Alignment, Cigar, CigarElem}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.reference.{FastaSequenceIndex, ReferenceSequenceFile, ReferenceSequenceFileFactory, ReferenceSequenceFileWalker}
import htsjdk.samtools.util.{SequenceUtil, StringUtil}
import htsjdk.samtools.{CigarOperator, SAMException, SamPairUtil}

import java.nio.file.{Files, Path}
import scala.collection.immutable
import scala.collection.mutable
import scala.util.{Failure, Success, Try}

/** What to do with an indel at the start of a read's alignment, for example after it is left-aligned there. */
sealed trait LeadingIndelMode extends EnumEntry

object LeadingIndelMode extends FgBioEnum[LeadingIndelMode] {
  /** Keep the indel at the start of the alignment. */
  case object Keep extends LeadingIndelMode
  /** Remove the indel from the alignment: inserted bases are soft-clipped, and the start of the alignment is moved past
    * deleted bases. */
  case object Trim extends LeadingIndelMode
  /** Remove every record of the read's template from the output. */
  case object Filter extends LeadingIndelMode

  override def values: immutable.IndexedSeq[LeadingIndelMode] = findValues
}

@clp(description =
  """
    |Left-aligns indels in a SAM/BAM file.
    |
    |An indel is shifted left only across repeated bases, such that every aligned read base remains paired with an
    |identical reference base.  Indels are not shifted into clipped bases or across reference skips (`N` CIGAR
    |operators); each block of the alignment between reference skips is left-aligned independently.  Within a run of
    |adjacent insertions and deletions, the insertions are written first (e.g. `2I1D`), so alignments that differ only
    |in the order of adjacent insertions and deletions are written identically.  Shifting may place an insertion directly
    |next to a deletion where the input had aligned bases between them; GATK's `GOOD_CIGAR` read filter (on by default in
    |HaplotypeCaller and Mutect2) removes reads with such adjacent indels.
    |
    |An indel at the start of the alignment (after any clipping), whether shifted there or already there in the input,
    |is not supported by any aligned bases before it, and GATK's `GOOD_CIGAR` read filter removes reads that start with a
    |deletion.  `--leading-indel-mode` controls what is done with such an indel:
    |
    |- `Trim` (the default): the indel is removed.  Inserted bases are soft-clipped, and the alignment start is moved
    |  past deleted bases.  The MD, NM, and UQ tags are re-computed if present.
    |- `Keep`: the indel is kept at the start of the alignment.
    |- `Filter`: all records of the read's template (the read, its mate, and their secondary and supplementary records)
    |  are removed from the output, so that no remaining record refers to a removed one.  A read with an unavailable
    |  query name (`*`) is not grouped with other reads, so only the read itself is removed.
    |
    |Indels at the start of a later block of the alignment (after a reference skip) are always kept.
    |
    |Unmapped reads, reads without stored bases (`SEQ` of `*`), and reads with padding (`P`) CIGAR operators are left
    |unchanged.  Read bases stored as `=` (equal to the reference) are compared using the reference base, and are written
    |as the reference base when the read's CIGAR changes.
    |
    |The reference bases are retrieved from the reference FASTA (`--ref`) when given.  An indexed FASTA (with a `.fai`)
    |may be used with any input sort order, while a FASTA without an index may only be used with coordinate sorted
    |input whose sequence dictionary lists the contigs in the same order as the FASTA.  The FASTA contig lengths must
    |match the input's sequence dictionary.  Without `--ref`, the reference bases are reconstructed from each read's MD
    |tag, which is then required for every mapped read that may be changed, and such reads may not have `=` read bases
    |that the MD tag reports as matching the reference.
    |
    |When a read's alignment changes, its MD, NM, and UQ tags are re-computed if present, counting a base as a mismatch
    |unless the read and reference bases are the same and not `N` (as samtools does).  UQ is removed if the read has no
    |base qualities.  The alignment score (AS) is not changed, even when indels are merged or trimmed.  The CIGAR uses
    |match (`=`) and mismatch (`X`) operators, with the same definition of a mismatch, when `--use-equals-and-x` is
    |given or the read's CIGAR already uses them, otherwise match-or-mismatch (`M`) operators.  A read whose CIGAR
    |already uses `=` and `X` operators keeps its CIGAR when no indel is moved.
    |
    |Mate information (e.g. the MC tag and insert size) is reset for all primary and supplementary reads, including
    |those whose alignment did not change.  The mate start, insert size, and MC tag of secondary reads are updated from
    |their mate's primary alignment, and the positions and CIGARs in SA tags are updated.  An SA tag entry is updated only when
    |exactly one alignment of the template matches it, otherwise it is left unchanged and counted.  To do so, input that
    |is not query name sorted or grouped is sorted into query name order, then sorted back into the input's sort order
    |(if known) when written; input of unknown order is written in query name order.  Query name sorted or grouped input
    |is written in the input order, except that template-coordinate grouped input is sorted again when trimming
    |leading indels may have changed the template's coordinates.  Reads with an unavailable query name (`*`) are not
    |grouped with other reads, so their mate information and SA tags are not updated.  Unpaired reads that share a query
    |name with another unpaired read are written without updating their SA tags.
    |
    |If the tool fails, the partially written output is removed.
  """,
  group = ClpGroups.SamOrBam)
class LeftAlignIndels
( @arg(flag='i', doc="Input BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
  @arg(flag='r', doc="Reference sequence FASTA file, indexed unless the input is coordinate sorted.  If not given, the MD tag is used.")
  val ref: Option[PathToFasta] = None,
  @arg(flag='x', doc="If specified, do not fail when reads marked as paired are missing their mate pairs.")
  val allowMissingMates: Boolean = false,
  @arg(doc="Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).")
  val useEqualsAndX: Boolean = false,
  @arg(flag='l', doc="What to do with an indel at the start of a read's alignment: trim it from the alignment, keep it, or filter the read's template.")
  val leadingIndelMode: LeadingIndelMode = LeadingIndelMode.Trim
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)
  validate(!Files.exists(output) || !Files.isSameFile(input, output), s"The output must not be the same file as the input: $output")

  override def execute(): Unit = {
    val reader       = SamSource(input)
    val header       = reader.header
    val isCoordinate = header.getSortOrder == SortOrder.coordinate
    val isQueryGroup = header.getSortOrder == SortOrder.queryname || header.getGroupOrder == GroupOrder.query
    val reference    = ref.map(path => (path, ReferenceSequenceFileFactory.getReferenceSequenceFile(path, true, true)))
    reference.filter { case (_, file) => !file.isIndexed && !isCoordinate }.foreach { case (path, file) =>
      reader.safelyClose()
      file.safelyClose()
      fail(s"The reference FASTA is not indexed, which requires a coordinate sorted input.  Index the FASTA (e.g. with " +
        s"samtools faidx) or sort the input by coordinate: $path")
    }

    val aligner       = new IndelLeftAligner(useEqualsAndX=useEqualsAndX, trimLeadingIndels=leadingIndelMode == LeadingIndelMode.Trim)
    val refBasesFor   = this.refBasesFor(reference, header, isCoordinate)
    val readProgress  = ProgressLogger(logger=logger, verb="read")
    val writeProgress = ProgressLogger(logger=logger, verb="written")
    var numChanged    = 0L // the number of records whose alignment was changed
    var numTrimmed    = 0L // the number of records whose leading indel was trimmed
    var numFiltered   = 0L // the number of records filtered for having (or being in a template with) a leading indel
    var numAmbiguous  = 0L // the number of SA tag entries not updated as they did not match exactly one alignment
    var numSharedName = 0L // the number of unpaired records that share a query name with another unpaired record

    // If the input is not query grouped, records are re-sorted into query name order to reset the mate information, and
    // so need to be re-sorted back into the input order (if known) when written.  If the input order is not known, the
    // output is written in query name order.  Since records are written in the order they are grouped, query grouped
    // input keeps its order (including any sub-sort, e.g. natural query name order), except that template-coordinate
    // input is re-sorted when trimming leading indels may change the coordinates of its templates.
    val (outHeader, outSort) = {
      if (isQueryGroup) {
        val trims = leadingIndelMode == LeadingIndelMode.Trim
        (header, SamOrder(header).filter(o => trims && o == SamOrder.TemplateCoordinate))
      }
      else SamOrder(header).filter(_ != SamOrder.Unsorted).orElse(if (isCoordinate) Some(SamOrder.Coordinate) else None) match {
        case Some(order) => (header, Some(order))
        case None        => (SamOrder.Queryname.applyTo(header.clone()), None)
      }
    }
    val writer = SamWriter(output, outHeader, sort=outSort)

    // Left-align each record before grouping them by query name, since a FASTA without an index must be read in order
    val records = reader.iterator.map { rec =>
      val result = try aligner.align(rec, refBasesFor(rec)) catch {
        case ex: InvalidRecordException => fail(ex.getMessage)
      }
      if (result.changed) numChanged += 1
      if (result.trimmed) numTrimmed += 1
      readProgress.record(rec)
      rec
    }

    var succeeded = false
    try {
      val grouped = Bams.queryGroupedIterator(records, header, Bams.MaxInMemory, Io.tmpDir)
      try {
        while (grouped.hasNext) {
          // Records without a query name (`*`) are not from the same template, so are not grouped
          val name  = grouped.head.name
          val group = if (name == "*") IndexedSeq(grouped.next()) else grouped.takeWhile(_.name == name).toIndexedSeq
          val keep  = leadingIndelMode != LeadingIndelMode.Filter || !group.exists(r => r.mapped && IndelLeftAligner.hasLeadingIndel(r.cigar))
          if (name != "*") Try(Template(group.iterator)) match {
            case Success(template) =>
              if (!allowMissingMates) {
                template.primaryReads.find(r => r.paired && (template.r1.isEmpty || template.r2.isEmpty)).foreach { r =>
                  fail(s"Record ${r.name} is paired but its mate was not found.  Use --allow-missing-mates to allow this.")
                }
              }
              if (keep) IndelLeftAligner.fixMateAndSupplementaryInfo(template, onUnmatchedSaEntry = _ => numAmbiguous += 1)
            case Failure(ex) =>
              if (group.exists(_.paired)) fail(s"Could not group the records with query name $name: ${ex.getMessage}")
              numSharedName += group.length
          }
          if (keep) group.foreach { rec =>
            writeProgress.record(rec)
            writer.write(rec)
          }
          else numFiltered += group.length
        }
      } finally {
        grouped.safelyClose()
      }
      succeeded = true
    } finally {
      reference.foreach { case (_, file) => file.safelyClose() }
      reader.safelyClose()
      if (succeeded) writer.close() else {
        // Do not leave a truncated, but otherwise valid, output behind
        writer.safelyClose()
        LeftAlignIndels.deleteOutput(output)
      }
    }
    readProgress.logLast()
    writeProgress.logLast()

    val total = readProgress.getCount
    def pct(count: Long): Double = if (total == 0) 0.0 else count * 100.0 / total
    logger.info(f"Changed the alignment of $numChanged%,d out of $total%,d records (${pct(numChanged)}%.2f%%).")
    leadingIndelMode match {
      case LeadingIndelMode.Trim   => logger.info(f"Trimmed a leading indel from $numTrimmed%,d records (${pct(numTrimmed)}%.2f%%).")
      case LeadingIndelMode.Filter => logger.info(f"Filtered $numFiltered%,d records (${pct(numFiltered)}%.2f%%) in templates with a leading indel.")
      case LeadingIndelMode.Keep   => ()
    }
    if (numAmbiguous > 0) logger.warning(f"Did not update $numAmbiguous%,d SA tag entries that did not match exactly one alignment.")
    if (numSharedName > 0) logger.warning(f"Did not update the SA tags of $numSharedName%,d unpaired records that share a query name.")
  }

  /** Returns the method to retrieve the upper-case reference bases spanned by a mapped record, from its start to its
    * end inclusive.
    *
    * @param reference the path to the reference FASTA and the opened FASTA, if given, otherwise the MD tag is used
    * @param header the input's header
    * @param isCoordinate true if the input is coordinate sorted, required for a FASTA without an index
    */
  private def refBasesFor(reference: Option[(PathToFasta, ReferenceSequenceFile)],
                          header: htsjdk.samtools.SAMFileHeader,
                          isCoordinate: Boolean): SamRecord => String = reference match {
    case Some((path, file)) if file.isIndexed =>
      logger.info("Using the indexed reference FASTA to retrieve reference bases.")
      val index   = new FastaSequenceIndex(ReferenceSequenceFileFactory.getFastaIndexFileName(path))
      val checked = mutable.Set[Int]() // the indices of the contigs whose length was checked
      rec => {
        val length = header.getSequence(rec.refIndex).getSequenceLength
        if (!checked.contains(rec.refIndex)) {
          if (!index.hasIndexEntry(rec.refName)) fail(s"Contig ${rec.refName} is not in the reference FASTA: $path")
          val fastaLength = index.getIndexEntry(rec.refName).getSize
          if (fastaLength != length) fail(
            s"Contig ${rec.refName} has length $fastaLength in the reference FASTA but $length in the input's sequence dictionary."
          )
          checked += rec.refIndex
        }
        if (rec.end > length) fail(s"Record ${rec.name} ends at ${rec.end}, past the end of contig ${rec.refName} (length $length).")
        upper(file.getSubsequenceAt(rec.refName, rec.start.toLong, rec.end.toLong).getBases, 0, rec.end - rec.start + 1)
      }
    case Some((_, file)) =>
      require(isCoordinate, "A reference FASTA without an index requires a coordinate sorted input.")
      logger.info("Input is coordinate sorted: using the reference FASTA (without an index) to retrieve reference bases.")
      val walker = new ReferenceSequenceFileWalker(file)
      rec => {
        val contig = try walker.get(rec.refIndex, rec.refName, header.getSequence(rec.refIndex).getSequenceLength) catch {
          case ex: SAMException => fail(
            s"Could not retrieve contig ${rec.refName} from the reference FASTA, whose contigs must match the input's " +
              s"sequence dictionary in name, length, and order (or index the FASTA): ${ex.getMessage}"
          )
        }
        if (rec.end > contig.length) {
          fail(s"Record ${rec.name} ends at ${rec.end}, past the end of contig ${rec.refName} (length ${contig.length}).")
        }
        upper(contig.getBases, rec.start - 1, rec.end - rec.start + 1)
      }
    case None =>
      logger.info("Reference FASTA not provided: using the MD tag to retrieve reference bases.")
      rec => {
        if (!rec.contains("MD")) fail(
          s"Record ${rec.name} has no MD tag.  Provide --ref, or add MD tags to the input (e.g. with samtools calmd)."
        )
        val bases = try rec.refBasesFromMd catch {
          case ex: IllegalArgumentException => fail(s"${ex.getMessage}.  Provide --ref, or fix the MD tags (e.g. with samtools calmd).")
        }
        upper(bases, 0, bases.length)
      }
  }

  /** Returns `length` bases starting at `offset` as an upper-case string. */
  private def upper(bases: Array[Byte], offset: Int, length: Int): String = {
    val copy = java.util.Arrays.copyOfRange(bases, offset, offset + length)
    StringUtil.toUpperCase(copy)
    StringUtil.bytesToString(copy)
  }
}

object LeftAlignIndels {
  /** Deletes the output and its index, if they are regular files, ignoring any failure. */
  private[bam] def deleteOutput(output: Path): Unit = {
    val name    = output.getFileName.toString
    val indices = Seq(output.resolveSibling(name + ".bai")) ++
      (if (name.endsWith(".bam")) Seq(output.resolveSibling(name.dropRight(4) + ".bai")) else Seq.empty)
    (output +: indices).filter(p => Files.isRegularFile(p)).foreach(p => Try(Files.delete(p)))
  }
}

/** Thrown when a record cannot be left-aligned because its data is inconsistent (e.g. with the reference). */
private[bam] class InvalidRecordException(message: String) extends RuntimeException(message)

/** The result of left-aligning a record.
  *
  * @param changed true if the record's alignment (CIGAR or start) was changed
  * @param trimmed true if an indel was trimmed from the start of the record's alignment
  */
private[bam] case class LeftAlignResult(changed: Boolean, trimmed: Boolean)

/** Left-aligns indels in [[com.fulcrumgenomics.bam.api.SamRecord]]s.
  *
  * Indels are shifted left only across repeat-equivalent bases, so that every aligned base keeps its pairing (and hence
  * match or mismatch status).  Indels are left-aligned independently within each block of the alignment separated by
  * reference skips (`N` CIGAR operators), and are never shifted into clipped bases.  When the alignment changes, the
  * MD, NM, and UQ tags are re-computed if present (see [[IndelLeftAligner.isMatch]]).
  *
  * @param useEqualsAndX use match (=) and mismatch (X) CIGAR operators instead of match-or-mismatch (M).  Records
  *                      whose CIGAR already uses = or X operators always keep using them.
  * @param trimLeadingIndels true to remove an indel at the start of the alignment (after any clipping), by
  *                          soft-clipping inserted bases and moving the start past deleted bases
  */
private[bam] class IndelLeftAligner(useEqualsAndX: Boolean, trimLeadingIndels: Boolean) {
  import IndelLeftAligner._

  /** Left-aligns the given [[com.fulcrumgenomics.bam.api.SamRecord]] given the reference bases spanned by the record
    * (from its start to its end, inclusive, including any deleted or skipped bases).  The reference bases are only
    * evaluated if the record may change.  Unmapped records, records without bases, and records with padding (P)
    * operators are not changed.  Read bases stored as `=` are compared using their reference base, and are written as
    * their reference base if the alignment changes.  A record whose CIGAR uses `=` and `X` operators is not changed
    * unless an indel moves or `useEqualsAndX` is true and its CIGAR also uses `M` operators.
    *
    * @throws InvalidRecordException if the record's CIGAR, bases, and reference bases are inconsistent
    */
  def align(rec: SamRecord, refBases: => String): LeftAlignResult = {
    val originalCigar = rec.cigar
    val usesEqualsAndX = originalCigar.exists(e => e.operator == CigarOperator.EQ || e.operator == CigarOperator.X)
    val emitEqAndX     = useEqualsAndX || usesEqualsAndX
    val mayChange      = originalCigar.exists(_.operator.isIndel) || (useEqualsAndX && originalCigar.exists(_.operator == CigarOperator.M))

    if (rec.unmapped || !mayChange || rec.length == 0 || originalCigar.exists(_.operator == CigarOperator.P)) Unchanged else {
      val target        = refBases
      val originalBases = rec.basesString
      if (originalCigar.lengthOnQuery != originalBases.length) throw new InvalidRecordException(
        s"Query length mismatch for record ${rec.name}: ${originalCigar.lengthOnQuery} (CIGAR) vs. ${originalBases.length} (bases)")
      if (originalCigar.lengthOnTarget != target.length) throw new InvalidRecordException(
        s"Reference length mismatch for record ${rec.name}: ${originalCigar.lengthOnTarget} (CIGAR) vs. ${target.length} " +
          "(reference bases).  Check that the reference FASTA or MD tag matches the alignment.")
      if (target.contains(IndelPaddedAlignment.Indel)) throw new InvalidRecordException(
        s"The reference bases for record ${rec.name} contain the gap character '${IndelPaddedAlignment.Indel}': $target")

      // Clipping may only be the first and/or last elements, and indels are not shifted into clipped bases
      val (leadingClips, inner, trailingClips) = splitClipping(originalCigar)

      // Left-align each block of the alignment separated by reference skips (N), keeping the reference skips in place
      val bases       = resolveMatchMarkers(originalCigar, originalBases, target)
      val newElems    = IndexedSeq.newBuilder[CigarElem]
      var queryOffset = leadingClips.map(_.lengthOnQuery).sum
      var refOffset   = 0
      var blockStart  = 0
      forloop(from=0, until=inner.length + 1) { i =>
        if (i == inner.length || inner(i).operator == CigarOperator.N) {
          if (blockStart < i) {
            val block = Cigar(inner.slice(blockStart, i))
            val query = bases.substring(queryOffset, queryOffset + block.lengthOnQuery)
            val ref   = target.substring(refOffset, refOffset + block.lengthOnTarget)
            newElems ++= IndelPaddedAlignment(block, query, ref).leftAligned.cigar(useEqualsAndX=emitEqAndX).elems
            queryOffset += block.lengthOnQuery
            refOffset   += block.lengthOnTarget
          }
          if (i < inner.length) {
            newElems  += inner(i)
            refOffset += inner(i).length
          }
          blockStart = i + 1
        }
      }

      // Optionally trim the leading indels of the first block, provided aligned bases remain in that block
      val alignedElems = newElems.result()
      val numLeading   = alignedElems.indexWhere(!_.operator.isIndel)
      val trim         = trimLeadingIndels && numLeading > 0 && alignedElems(numLeading).operator.isAlignment
      val (newLeadingClips, newInner, numTrimmedDeleted) = if (!trim) (leadingClips, alignedElems, 0) else {
        val leading  = alignedElems.take(numLeading)
        val inserted = leading.filter(_.operator == CigarOperator.I).map(_.length).sum
        val deleted  = leading.filter(_.operator == CigarOperator.D).map(_.length).sum
        val clips    = if (inserted == 0) leadingClips else leadingClips.lastOption match {
          case Some(CigarElem(CigarOperator.S, len)) => leadingClips.dropRight(1) :+ CigarElem(CigarOperator.S, len + inserted)
          case _                                     => leadingClips :+ CigarElem(CigarOperator.S, inserted)
        }
        (clips, alignedElems.drop(numLeading), deleted)
      }
      val newCigar = Cigar(newLeadingClips ++ newInner ++ trailingClips)

      // A CIGAR whose indels did not move is kept, unless its M operators are converted to = and X
      val keep = newCigar == originalCigar || (!trim && withoutEqualsAndX(newCigar) == withoutEqualsAndX(originalCigar) &&
        !(useEqualsAndX && originalCigar.exists(_.operator == CigarOperator.M)))
      if (keep) Unchanged else {
        rec.cigar  = newCigar
        rec.start += numTrimmedDeleted
        // An `=` read base may no longer be aligned to its reference base, so write the resolved bases
        if (bases != originalBases) rec.bases = bases
        if (rec.contains("MD") || rec.contains("NM") || rec.contains("UQ")) {
          val tags = AlignmentTags(newCigar, bases, rec.quals, target, targetOffset=numTrimmedDeleted)
          if (rec.contains("MD")) rec("MD") = tags.md
          if (rec.contains("NM")) rec("NM") = tags.nm
          if (rec.contains("UQ")) rec("UQ") = tags.uq.map(Int.box).orNull
        }
        LeftAlignResult(changed=true, trimmed=trim)
      }
    }
  }
}

private[bam] object IndelLeftAligner {
  /** The result for a record that was not changed. */
  private val Unchanged = LeftAlignResult(changed=false, trimmed=false)

  /** True if the read and reference bases match: they are the same base (ignoring case) and not `N`, as in samtools. */
  def isMatch(readBase: Char, refBase: Char): Boolean = {
    val base = Character.toUpperCase(readBase)
    base == Character.toUpperCase(refBase) && base != 'N'
  }

  /** True if the alignment starts with an insertion or deletion (after any clipping). */
  def hasLeadingIndel(cigar: Cigar): Boolean = cigar.elems.find(!_.operator.isClipping).exists(_.operator.isIndel)

  /** Splits the CIGAR into its leading clipping elements, the elements between the clipping, and its trailing clipping
    * elements. */
  def splitClipping(cigar: Cigar): (IndexedSeq[CigarElem], IndexedSeq[CigarElem], IndexedSeq[CigarElem]) = {
    val (leading, rest) = cigar.elems.span(_.operator.isClipping)
    val numTrailing     = rest.reverseIterator.takeWhile(_.operator.isClipping).length
    (leading, rest.dropRight(numTrailing), rest.takeRight(numTrailing))
  }

  /** Returns the CIGAR with its = and X operators replaced with M operators. */
  private def withoutEqualsAndX(cigar: Cigar): Cigar = Cigar(cigar.elems.map { e =>
    if (e.operator == CigarOperator.EQ || e.operator == CigarOperator.X) e.copy(operator=CigarOperator.M) else e
  }).coalesce

  /** Returns the read bases with each `=` (equal to the reference) in an aligned (M, =, or X) position replaced by the
    * reference base it is aligned to.
    *
    * @param cigar the read's CIGAR
    * @param bases the read's bases
    * @param target the reference bases spanned by the read, including any deleted or skipped bases
    */
  def resolveMatchMarkers(cigar: Cigar, bases: String, target: String): String = {
    if (!bases.contains('=')) bases else {
      val chars       = bases.toCharArray
      var queryOffset = 0
      var refOffset   = 0
      cigar.foreach { elem =>
        if (elem.operator.isAlignment) forloop(from=0, until=elem.length) { i =>
          if (chars(queryOffset + i) == '=') chars(queryOffset + i) = target.charAt(refOffset + i)
        }
        queryOffset += elem.lengthOnQuery
        refOffset   += elem.lengthOnTarget
      }
      new String(chars)
    }
  }

  /** The alignment-dependent tags of a record: the MD tag, the edit distance (NM), and the sum of the base qualities of
    * mismatched bases (UQ, if the record has base qualities). */
  private[bam] case class AlignmentTags(md: String, nm: Int, uq: Option[Int])

  private[bam] object AlignmentTags {
    /** Computes the MD, NM, and UQ tags, counting a base as a mismatch as given by [[isMatch]].
      *
      * @param cigar the record's CIGAR
      * @param bases the record's bases
      * @param quals the record's base qualities, or empty if not stored
      * @param target reference bases that include those spanned by the alignment
      * @param targetOffset the offset in `target` of the first reference base of the alignment
      */
    def apply(cigar: Cigar, bases: String, quals: Array[Byte], target: String, targetOffset: Int): AlignmentTags = {
      val md          = new StringBuilder
      var matches     = 0
      var nm          = 0
      var uq          = 0
      var queryOffset = 0
      var refOffset   = targetOffset
      cigar.foreach { elem =>
        elem.operator match {
          case op if op.isAlignment =>
            forloop(from=0, until=elem.length) { i =>
              val refBase = target.charAt(refOffset + i)
              if (isMatch(bases.charAt(queryOffset + i), refBase)) matches += 1 else {
                md.append(matches).append(refBase)
                matches = 0
                nm     += 1
                if (quals.nonEmpty) uq += quals(queryOffset + i)
              }
            }
          case CigarOperator.D =>
            md.append(matches).append('^').append(target.substring(refOffset, refOffset + elem.length))
            matches = 0
            nm     += elem.length
          case CigarOperator.I => nm += elem.length
          case _               => () // S, H, N, and P
        }
        queryOffset += elem.lengthOnQuery
        refOffset   += elem.lengthOnTarget
      }
      md.append(matches)
      AlignmentTags(md=md.toString, nm=nm, uq=if (quals.isEmpty) None else Some(uq))
    }
  }

  /** The properties of a CIGAR and its end position that left-aligning (including trimming leading indels) does not
    * change, used to match an SA tag entry (which has the alignment before left-aligning) to the alignment it refers to.
    * Hard and soft clipping are not distinguished, since the SA tag may list a hard-clipped alignment with
    * soft-clipping. */
  private case class AlignmentInvariants(end: Int, leadingClipOrInserted: Int, trailingClip: Int, aligned: Int, skipped: Int)

  private object AlignmentInvariants {
    def from(start: Int, cigar: Cigar): AlignmentInvariants = {
      val (leading, inner, trailing) = splitClipping(cigar)
      def sum(elems: Seq[CigarElem], p: CigarElem => Boolean): Int = elems.filter(p).map(_.length).sum
      AlignmentInvariants(
        end                   = start + cigar.lengthOnTarget - 1,
        leadingClipOrInserted = sum(leading, _ => true) + sum(inner, _.operator == CigarOperator.I),
        trailingClip          = sum(trailing, _ => true),
        aligned               = sum(inner, _.operator.isAlignment),
        skipped               = sum(inner, _.operator == CigarOperator.N)
      )
    }
  }

  private def isR1(rec: SamRecord): Boolean = !rec.paired || rec.firstOfPair

  private def strand(rec: SamRecord): String = if (rec.positiveStrand) "+" else "-"

  /** Fixes the mate information on the primary and supplementary records of a template, sets the mate start and mate
    * CIGAR (MC) on secondary records from their mate's primary record, and updates the positions and CIGARs in the SA
    * tags of all records from the primary and supplementary records they refer to.
    *
    * @param template the template to fix
    * @param onUnmatchedSaEntry called with each SA tag entry that was not updated, as it did not match exactly one
    *                           alignment
    */
  def fixMateAndSupplementaryInfo(template: Template, onUnmatchedSaEntry: String => Unit = _ => ()): Unit = {
    template.fixMateInfo()

    // Update the mate start, insert size, and mate cigar on secondary records from the mate's primary record.  The secondary's mate
    // start may be before the mate's start was moved by trimming a leading deletion, which does not change its end.
    def fixSecondaries(secondaries: Seq[SamRecord], matePrimary: Option[SamRecord]): Unit = matePrimary.filter(_.mapped).foreach { mate =>
      val mateCigar = mate.cigar.toString
      secondaries.filter { r =>
        r.paired && r.mateMapped && r.mateRefIndex == mate.refIndex &&
          (r.mateStart == mate.start || r.mateCigar.exists(c => r.mateStart + c.lengthOnTarget - 1 == mate.end))
      }.foreach { r =>
        r.mateStart  = mate.start
        r.insertSize = SamPairUtil.computeInsertSize(r.asSam, mate.asSam)
        r("MC")      = mateCigar
      }
    }
    fixSecondaries(template.r1Secondaries, template.r2)
    fixSecondaries(template.r2Secondaries, template.r1)

    // Update the SA tags, which list the other primary or supplementary alignments for the same read
    val records = template.allReads.filter(_.contains("SA")).toIndexedSeq
    if (records.nonEmpty) {
      val alignments = (template.primaryReads ++ template.r1Supplementals ++ template.r2Supplementals).filter(_.mapped).toIndexedSeq
      records.foreach { rec =>
        val entries = rec[String]("SA").split(';').filter(_.nonEmpty).map { entry =>
          updateSaEntry(entry, rec, alignments).getOrElse { onUnmatchedSaEntry(entry); entry }
        }
        rec("SA") = entries.mkString("", ";", ";")
      }
    }
  }

  /** Returns the SA tag entry with its position and CIGAR replaced by those of the alignment it refers to, or None if
    * the entry cannot be parsed or does not match exactly one alignment.  The alignment must be for the same read as
    * `rec` (but not `rec` itself), on the same reference and strand, and match the entry in all properties that
    * left-aligning does not change (see [[AlignmentInvariants]]).  The entry's clipping operators are kept (e.g. soft
    * clipping for a hard-clipped alignment) where its clipping at that end has no hard clipping.
    *
    * @param entry a single entry from the SA tag of `rec`
    * @param rec the record with the SA tag
    * @param alignments the primary and supplementary alignments of the template
    */
  private def updateSaEntry(entry: String, rec: SamRecord, alignments: Seq[SamRecord]): Option[String] = entry.split(',') match {
    case Array(refName, pos, strandString, cigarString, mapq, nm) =>
      Try((pos.toInt, Cigar(cigarString))).toOption.flatMap { case (entryStart, entryCigar) =>
        val invariants = AlignmentInvariants.from(entryStart, entryCigar)
        val matching   = alignments.filter { r =>
          !(r eq rec) && isR1(r) == isR1(rec) && r.refName == refName && strand(r) == strandString &&
            AlignmentInvariants.from(r.start, r.cigar) == invariants
        }
        matching match {
          case Seq(alignment) =>
            val (entryLeading, _, entryTrailing) = splitClipping(entryCigar)
            val (leading, inner, trailing)       = splitClipping(alignment.cigar)
            // Use soft clipping for an end whose clipping in the entry has no hard clipping
            def clipping(entryClips: Seq[CigarElem], clips: IndexedSeq[CigarElem]): IndexedSeq[CigarElem] = {
              val length = clips.map(_.length).sum
              if (length == 0 || entryClips.exists(_.operator == CigarOperator.H)) clips
              else IndexedSeq(CigarElem(CigarOperator.S, length))
            }
            val cigar = Cigar(clipping(entryLeading, leading) ++ inner ++ clipping(entryTrailing, trailing))
            Some(Seq(refName, alignment.start.toString, strandString, cigar.toString, mapq, nm).mkString(","))
          case _ => None
        }
      }
    case _ => None
  }
}

private[bam] object IndelPaddedAlignment {
  /** The character in the read or reference for insertions or deletions respectively. */
  val Indel: Char = '-'

  /** True if the two bases are the same (ignoring case), so that shifting an indel across them keeps every aligned base
    * pair the same. */
  def sameBase(a: Char, b: Char): Boolean = SequenceUtil.basesEqual(a.toByte, b.toByte)

  /** Builds a padded alignment for the given [[com.fulcrumgenomics.alignment.Cigar]], query and target bases.  The CIGAR
    * may contain only alignment (M, =, X), insertion (I), and deletion (D) operators: clipping must be removed, and each
    * block between reference skips (N) aligned separately.  The query and target may not contain the gap character. */
  def apply(cigar: Cigar, query: String, target: String): IndelPaddedAlignment = {
    require(cigar.forall(e => e.operator.isAlignment || e.operator.isIndel),
      s"only alignment (M, =, X), insertion (I), and deletion (D) operators are supported: $cigar")
    require(cigar.lengthOnQuery == query.length,
      s"query length mismatch: ${cigar.lengthOnQuery} vs. ${query.length}")
    require(cigar.lengthOnTarget == target.length,
      s"target length mismatch: ${cigar.lengthOnTarget} vs. ${target.length}")
    require(!query.contains(Indel) && !target.contains(Indel),
      s"the query and target may not contain the gap character '$Indel': $query vs. $target")

    val Seq(paddedQuery, _, paddedTarget) =
      Alignment(query=query, target=target, queryStart=1, targetStart=1, cigar=cigar, score=0).paddedString(padChar=Indel)
    IndelPaddedAlignment(paddedQuery, paddedTarget)
  }
}

/** Stores an alignment padded for indels.
  *
  * The query and target strings are padded as follows:
  * - an inserted base (in the query) has a corresponding gap character ('-') in the target
  * - a deleted base (in the target) has a corresponding gap character ('-') in the query
  *
  * Clipped bases are not part of the padded alignment.  For example:
  *
  * query:  GATTACAGATT--ACA
  * target: GA--ACAGATTAAACA
  *
  * In this example, the query bases in columns 3-4 of the padded alignment ("TT") are inserted, and the target bases in
  * columns 12-13 ("AA") are deleted.
  */
private[bam] case class IndelPaddedAlignment(query: String, target: String) {
  import IndelPaddedAlignment.{Indel, sameBase}
  require(query.length == target.length,
    f"Lengths differ ${query.length} != ${target.length}\n$query\n$target")
  require(length > 0)
  require((0 until length).forall(i => query(i) != Indel || target(i) != Indel), s"a column may not be a gap in both:\n$query\n$target")

  def length: Int = query.length

  /** Returns the CIGAR for this alignment, using match (=) and mismatch (X) operators (as given by
    * [[IndelLeftAligner.isMatch]]) when `useEqualsAndX` is true, otherwise match-or-mismatch (M) operators. */
  def cigar(useEqualsAndX: Boolean = false): Cigar = {
    def operatorFrom(index: Int): CigarOperator = {
      if (target(index) == Indel) CigarOperator.I
      else if (query(index) == Indel) CigarOperator.D
      else if (!useEqualsAndX) CigarOperator.M
      else if (IndelLeftAligner.isMatch(query(index), target(index))) CigarOperator.EQ
      else CigarOperator.X
    }

    val elems = IndexedSeq.newBuilder[CigarElem]
    var op    = operatorFrom(0)
    var len   = 1
    forloop(from=1, until=length) { index =>
      val nextOp = operatorFrom(index)
      if (nextOp == op) len += 1 else {
        elems += CigarElem(op, len)
        op  = nextOp
        len = 1
      }
    }
    elems += CigarElem(op, len)
    Cigar(elems.result())
  }

  /** Builds a left-aligned padded alignment from this alignment.
    *
    * Within each run of adjacent insertion and deletion columns, the insertions are first placed before the deletions,
    * which does not change the order of the query or target bases, nor any aligned base pair.  Then, the insertions of
    * a run are shifted left by one column if the aligned column before the run has a query base equal to the last
    * inserted base; otherwise the deletions of the run are shifted left (before the insertions) by one column if the
    * aligned column before the run has a target base equal to the last deleted base.  Either shift keeps every aligned
    * base paired with an identical base, so the matches and mismatches are unchanged.  This repeats until no indel can
    * be shifted, merging indels that become adjacent.  Since every run is first put in the same order, alignments that
    * differ only in the order of adjacent insertions and deletions give the same result.
    */
  def leftAligned: IndelPaddedAlignment = {
    val queryArr  = this.query.toCharArray
    val targetArr = this.target.toCharArray

    def isIndel(index: Int): Boolean = queryArr(index) == Indel || targetArr(index) == Indel

    // Returns the index after the end of the run of indel columns that starts at `start`
    def runEnd(start: Int): Int = {
      var end = start
      while (end < length && isIndel(end)) end += 1
      end
    }

    // Within each run of adjacent indel columns, places the insertions before the deletions
    def reorderAll(): Unit = {
      var index = 0
      while (index < length) {
        if (!isIndel(index)) index += 1 else {
          val end      = runEnd(index)
          val inserted = (index until end).filter(i => targetArr(i) == Indel).map(queryArr)
          val deleted  = (index until end).filter(i => queryArr(i) == Indel).map(targetArr)
          inserted.indices.foreach { i => queryArr(index + i) = inserted(i); targetArr(index + i) = Indel }
          deleted.indices.foreach  { i => queryArr(index + inserted.length + i) = Indel; targetArr(index + inserted.length + i) = deleted(i) }
          index = end
        }
      }
    }

    // Shifts the insertions or deletions of each run left by one column (the runs must be ordered with insertions
    // first), returning true if any were shifted
    def shiftAll(): Boolean = {
      var shifted = false
      var start   = 0
      while (start < length) {
        if (!isIndel(start)) start += 1 else {
          val end         = runEnd(start)
          val numInserted = (start until end).count(i => targetArr(i) == Indel)
          val numDeleted  = end - start - numInserted
          val prev        = start - 1 // the aligned column before the run, if any
          val inserted    = queryArr.slice(start, start + numInserted)
          val deleted     = targetArr.slice(start + numInserted, end)
          if (prev >= 0 && numInserted > 0 && sameBase(queryArr(prev), inserted.last)) {
            // [prev] [inserted] [deleted] => [prev query base +: inserted.init] [inserted.last / prev target base] [deleted]
            val prevQuery = queryArr(prev)
            val prevTarget = targetArr(prev)
            forloop(from=0, until=numInserted) { i =>
              queryArr(prev + i)  = if (i == 0) prevQuery else inserted(i - 1)
              targetArr(prev + i) = Indel
            }
            queryArr(prev + numInserted)  = inserted.last
            targetArr(prev + numInserted) = prevTarget
            shifted = true
          }
          else if (prev >= 0 && numDeleted > 0 && sameBase(targetArr(prev), deleted.last)) {
            // [prev] [inserted] [deleted] => [prev target base +: deleted.init] [prev query base / deleted.last] [inserted]
            val prevQuery  = queryArr(prev)
            val prevTarget = targetArr(prev)
            forloop(from=0, until=numDeleted) { i =>
              queryArr(prev + i)  = Indel
              targetArr(prev + i) = if (i == 0) prevTarget else deleted(i - 1)
            }
            queryArr(prev + numDeleted)  = prevQuery
            targetArr(prev + numDeleted) = deleted.last
            forloop(from=0, until=numInserted) { i =>
              queryArr(prev + numDeleted + 1 + i)  = inserted(i)
              targetArr(prev + numDeleted + 1 + i) = Indel
            }
            shifted = true
          }
          start = end
        }
      }
      shifted
    }

    // Each shift moves indel columns left (and placing insertions first does not move any indel column), so this ends
    var shifted = true
    while (shifted) {
      reorderAll()
      shifted = shiftAll()
    }

    IndelPaddedAlignment(new String(queryArr), new String(targetArr))
  }
}
