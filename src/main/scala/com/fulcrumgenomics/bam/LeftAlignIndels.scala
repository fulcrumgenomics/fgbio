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
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioCommonArgs, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.io.HtsPath
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.util.{FileExtensions, IOUtil, SequenceUtil, StringUtil}
import htsjdk.samtools.{CigarOperator, GenomicIndexUtil, SAMException, SAMFileHeader, SamPairUtil, SamReader}

import java.nio.file.Files
import scala.collection.immutable
import scala.collection.mutable
import scala.jdk.CollectionConverters._
import scala.util.{Failure, Success, Try, Using}

/** What to do with an indel at the start of a read's alignment, for example after it is left-aligned there. */
sealed trait LeadingIndelMode extends EnumEntry

object LeadingIndelMode extends FgBioEnum[LeadingIndelMode] {
  /** Keep the indel at the start of the alignment. */
  case object Keep extends LeadingIndelMode
  /** Remove the indel from the alignment: inserted bases are soft-clipped, and the start of the alignment is moved past
    * deleted bases.  The indel is kept if no aligned bases follow it before the first reference skip (`N`). */
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
    |  past deleted bases.  The MD, NM, and UQ tags are re-computed if present.  The indel is kept if no aligned bases
    |  follow it before the first reference skip (e.g. `2I10N2M`).
    |- `Keep`: the indel is kept at the start of the alignment.
    |- `Filter`: all records of the read's template (the read, its mate, and their secondary and supplementary records)
    |  are removed from the output, so that no remaining record refers to a removed one.  A read with an unavailable
    |  query name (`*`) is not grouped with other reads, so only the read itself is removed.
    |
    |Indels at the start of a later block of the alignment (after a reference skip) are always kept.
    |
    |Unmapped reads, reads without stored bases (`SEQ` of `*`), and reads with padding (`P`) CIGAR operators are left
    |unchanged.  Read bases stored as `=` (equal to the reference) are compared using the reference base, and are written
    |as the reference base when the read's CIGAR changes.  If a resolved base cannot be stored unchanged in SAM/BAM
    |(e.g. `U`, which is normalized to `T`), the tool fails with the read's name.  The tool also fails if materializing
    |`=` changes SEQ on a record with base-modification tags (MM, Mm, ML, Ml, or MN), rather than shifting modification
    |calls.  These tags are preserved when SEQ is unchanged.  On contigs marked `TP:circular`, reference lookup and
    |trimmed alignment starts wrap at the contig boundary.
    |
    |When `--ref` is given, the entire reference FASTA is loaded into memory before processing reads.  Plain, gzip, and
    |BGZF FASTAs are supported; no index is needed, and any `.fai`, `.gzi`, or `.dict` sidecars are ignored.  Every contig
    |in the input's sequence dictionary must be present with the same length, even if it has no reads.  Additional FASTA
    |contigs (which are read but not kept in memory) and a different contig order are allowed.  A gzip or BGZF FASTA must
    |have a `.gz` or `.gzip` extension (e.g. not `.bgz`); a plain FASTA need not use a standard filename extension.  The
    |reference is loaded and checked even for an empty input or one with no alignments that need changing.
    |
    |A human reference needs about 3.1 GB for its bases alone, plus memory for parsing, reads, and sorting buffers.
    |The launcher's default heap may be too small.  Use an explicit larger heap (e.g. `fgbio -Xmx8g ...` or
    |`java -Xmx8g -jar fgbio.jar ...`) as a starting point; larger references, reads, or tags may need more.
    |
    |Without `--ref`, the reference bases are reconstructed from each read's MD tag, which is required for every mapped
    |read that may be changed.  Such reads may not have `=` read bases that the MD tag reports as matching the reference.
    |
    |When a read's alignment changes, its MD, NM, and UQ tags are re-computed if present, counting a base as a mismatch
    |unless the read and reference bases are the same, ignoring case (as HTSJDK and other fgbio tools do).  An `N` (or
    |another IUPAC code) in both the read and the reference is therefore a match.  This differs from the SAM
    |specification's NM definition, under which only `A`, `C`, `G`, and `T` bases can match.  With
    |`--samtools-mismatches`, a read and reference base that are both `N` are a mismatch instead (as samtools does).
    |Repaired SA tag entries use the same edit distance.  UQ is removed if the read has no base qualities.  The BQ,
    |minimap2 cs, ds, and de tags are removed when the alignment edit path changes, but kept when only M versus =/X
    |notation changes.  The de tag is also removed when adjacent indel operators are coalesced.
    |The alignment score (AS) is not changed, even when indels are merged or trimmed.  The CIGAR uses match (`=`) and mismatch (`X`) operators, with the same definition of a mismatch, when
    |`--use-equals-and-x` is given or the read's CIGAR already uses them, otherwise match-or-mismatch (`M`) operators.
    |Existing `=` and `X` notation is kept when no indel is moved or coalesced, unless remaining `M` operators are converted.
    |
    |Mate information (e.g. the MC tag and insert size) is reset for all primary and supplementary reads, including
    |those whose alignment did not change; existing ms scores are preserved and absent ms tags are not added.  The mate
    |start, insert size, and MC tag of secondary reads are updated from their mate's primary alignment only when its
    |original reference, start, and strand agree, and its MC (if present) exactly matches the original primary CIGAR.
    |The positions, CIGARs, and edit distances in SA tags are updated.  An SA tag entry is matched against the original start and CIGAR, accepting exact CIGARs and minimap2-style summaries
    |(end clips, the smaller of the aligned query and reference spans as M, then any net I or D).  Hard and soft clips,
    |and M versus =/X, are treated equivalently.  An entry is updated only when exactly one alignment of the template
    |matches either form, otherwise it is left unchanged and counted.
    |
    |Records are grouped by query name before any alignments are changed.  Coordinate and template-coordinate input is
    |sorted into query name order, then sorted into coordinate or template-coordinate order after repairs, in every
    |leading indel mode.  These two sorts require additional memory and temporary disk space.  The template-coordinate
    |order is fgbio's, so may order ties differently than the input (e.g. if sorted by samtools).  Template-coordinate
    |`GO:query` is not trusted because secondary or supplementary alignments may be separated from their primary pair.
    |Other query name sorted input, or `SO:unsorted GO:query` input, keeps its input order and ordering metadata.  Random
    |input is re-sorted into random order.  Unknown or ordinary unsorted input is written in query name order, including
    |`SO:unknown GO:query` (or `GO:query` without `SO`): the unknown sort order takes precedence over the grouping
    |declaration.
    |
    |With `--leading-indel-mode Trim`, trimming changes the coordinates of a read's template, so in template-coordinate
    |output (e.g. from `GroupReadsByUmi`) the reads of a molecule (`MI` tag) may no longer be adjacent, as required by
    |the consensus callers.  Left-align indels before grouping reads by UMI, or use `Keep` or `Filter`, which do not
    |change the coordinates of the reads that are written.
    |
    |For template-coordinate output, every retained paired record with a mapped mate must have an MC tag after repairs,
    |including in `Keep` and `Filter` modes.  If one does not, the tool fails with a diagnostic instead of changing the
    |output order.  Add mate CIGARs (e.g. with `samtools fixmate -m` or fgbio `SetMateInformation` on query name grouped
    |input), or query name sort the input to obtain query name output.  `--allow-missing-mates` does not waive this
    |requirement.
    |
    |Reads with an unavailable query name (`*`) are not grouped with other reads, so their mate information and SA tags
    |are not updated.  Unpaired reads that share a query name with another unpaired read are written without updating
    |their SA tags, and with `--leading-indel-mode Filter` only those reads with a leading indel are removed.
  """,
  group = ClpGroups.SamOrBam)
class LeftAlignIndels
( @arg(flag='i', doc="Input BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
  @arg(flag='r', doc="Reference sequence FASTA file, loaded entirely into memory without requiring an index.  If not given, the MD tag is used.")
  val ref: Option[PathToFasta] = None,
  @arg(flag='x', doc="If specified, do not fail when reads marked as paired are missing their mate pairs.")
  val allowMissingMates: Boolean = false,
  @arg(doc="Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).")
  val useEqualsAndX: Boolean = false,
  @arg(flag='l', doc="What to do with an indel at the start of a read's alignment: trim it from the alignment, keep it, or filter the read's template.")
  val leadingIndelMode: LeadingIndelMode = LeadingIndelMode.Trim,
  @arg(doc="Count a read base and reference base that are both N as a mismatch (as samtools does) in the MD, NM, and UQ tags and the =/X operators, instead of as a match (as HTSJDK and other fgbio tools do).")
  val samtoolsMismatches: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)
  validate(!Files.exists(output) || !Files.isSameFile(input, output), s"The output must not be the same file as the input: $output")
  validate(!Files.exists(output) || ref.forall(r => !Files.isSameFile(r, output)),
    s"The output must not be the same file as the reference FASTA: $output")

  override def execute(): Unit = {
    // Common arguments are installed after tool construction; this reference may be used by either CRAM reader or writer.
    val references = (ref.toIndexedSeq ++ FgBioCommonArgs.args.cramRefFasta).filter(Files.exists(_))
    if (Files.exists(output) && references.exists(r => Files.isSameFile(r, output))) fail(
      s"The output must not be the same file as a reference FASTA: $output")

    val readProgress  = ProgressLogger(logger=logger, verb="read")
    val writeProgress = ProgressLogger(logger=logger, verb="written")
    var numChanged    = 0L // the number of records whose alignment was changed
    var numTrimmed    = 0L // the number of records whose leading indel was trimmed
    var numFiltered   = 0L // the number of records filtered for having (or being in a template with) a leading indel
    var numAmbiguous  = 0L // the number of SA tag entries not updated as they did not match exactly one alignment
    var numSharedName = 0L // the number of unpaired records that share a query name with another unpaired record

    Using.resource(SamSource(input)) { reader =>
      val header       = reader.header
      // A missing SO means the sort order is unknown (per the SAM spec), although htsjdk reports it as unsorted
      val sortOrder    = if (header.getAttribute("SO") == null) SortOrder.unknown else header.getSortOrder
      val groupOrder   = Option(header.getGroupOrder).getOrElse(GroupOrder.none)
      val isCoordinate = sortOrder == SortOrder.coordinate
      // The sub-sort may or may not be prefixed by the sort order (e.g. `unsorted:template-coordinate`), and may be
      // refined by further colon-separated values (e.g. `unsorted:template-coordinate:<more>`)
      def hasSubSort(order: SamOrder): Boolean = order.subSort.exists { o =>
        Option(header.getAttribute("SS")).exists { ss =>
          Seq(ss, ss.substring(ss.indexOf(':') + 1)).exists(s => s == o || s.startsWith(o + ":"))
        }
      }
      val isTemplateCoordinate = sortOrder == SortOrder.unsorted && groupOrder == GroupOrder.query &&
        hasSubSort(SamOrder.TemplateCoordinate)
      val isQueryGroup = sortOrder == SortOrder.queryname ||
        (sortOrder == SortOrder.unsorted && groupOrder == GroupOrder.query && !isTemplateCoordinate)

      val aligner     = new IndelLeftAligner(useEqualsAndX=useEqualsAndX, trimLeadingIndels=leadingIndelMode == LeadingIndelMode.Trim,
        samtoolsMismatches=samtoolsMismatches)
      val refBasesFor = this.refBasesFor(header)

      // Restore positional orders after repairs, even when the input also declares GO:query.  Only genuinely
      // query-grouped input keeps its order and header; the writer must not mutate the header used for grouping.
      val outSort: Option[SamOrder] = {
        if (isCoordinate) Some(SamOrder.Coordinate)
        else if (isTemplateCoordinate) Some(SamOrder.TemplateCoordinate)
        else if (sortOrder == SortOrder.unsorted && groupOrder == GroupOrder.none && hasSubSort(SamOrder.Random)) Some(SamOrder.Random)
        else None
      }
      val outHeader = if (!isQueryGroup && outSort.isEmpty) SamOrder.Queryname.applyTo(header.clone()) else header.clone()

      // Match SamWriter's indexing gate and HTSJDK's selected index name, not other possible sidecars.
      if (outHeader.getSortOrder == SortOrder.coordinate && IOUtil.isRegularPath(output) &&
        outHeader.getSequenceDictionary.getSequences.asScala.forall(_.getSequenceLength <= GenomicIndexUtil.BIN_GENOMIC_SPAN)) {
        val name = output.getFileName.toString
        val lowerName = name.toLowerCase(java.util.Locale.ROOT)
        val isSam = SamReader.Type.SAM_TYPE.hasValidFileExtension(name) ||
          FileExtensions.BLOCK_COMPRESSED.asScala.exists(ext => lowerName.endsWith(".sam" + ext))
        val index = if (SamReader.Type.CRAM_TYPE.hasValidFileExtension(name)) Some(IOUtil.addExtension(output, FileExtensions.CRAM_INDEX))
        else if (isSam) {
          if (IOUtil.hasBlockCompressedExtension(output)) Some(IOUtil.addExtension(output, FileExtensions.CSI)) else None
        }
        else if (SamReader.Type.BAM_TYPE.hasValidFileExtension(name)) {
          Some(output.resolveSibling(name.dropRight(FileExtensions.BAM.length) + FileExtensions.BAI_INDEX))
        }
        else Some(IOUtil.addExtension(output, FileExtensions.BAI_INDEX))
        index.filter(Files.exists(_)).foreach { path =>
          if (Files.isSameFile(input, path) || references.exists(r => Files.isSameFile(r, path))) fail(
            s"The output index must not be the same file as the input or reference FASTA: $path")
        }
      }
      Using.resource(SamWriter(output, outHeader, sort=outSort)) { out =>

        val records = reader.iterator.map { rec =>
          readProgress.record(rec)
          rec
        }

        // Original-alignment snapshots are transient, so align only after the grouping sort is complete.
        def alignRecords(group: IndexedSeq[SamRecord]): Unit = group.foreach { rec =>
          val (start, cigar) = (rec.start, rec.cigar)
          val result = try aligner.align(rec, refBasesFor(rec)) catch {
            case ex: InvalidRecordException => fail(ex.getMessage)
          }
          if (result.changed) {
            numChanged += 1
            IndelLeftAligner.setOriginalAlignment(rec, start, cigar, result.nm)
          }
          if (result.trimmed) numTrimmed += 1
        }

        val filter = leadingIndelMode == LeadingIndelMode.Filter
        def hasLeadingIndel(rec: SamRecord): Boolean = rec.mapped && IndelLeftAligner.hasLeadingIndel(rec.cigar)

        val grouped = if (isQueryGroup) Bams.queryGroupedIterator(records, header, Bams.MaxInMemory, Io.tmpDir)
                      else Bams.querySortedIterator(records, header, Bams.MaxInMemory, Io.tmpDir)
        try {
          while (grouped.hasNext) {
            // Records without a query name (`*`) are not from the same template, so are not grouped
            val name  = grouped.head.name
            val group = if (name == "*") IndexedSeq(grouped.next()) else grouped.takeWhile(_.name == name).toIndexedSeq
            // Unpaired records that share a query name with paired records are not part of the pair (e.g. they must not
            // be used as the mate of a paired record whose own mate is missing), so are aligned and filtered on their own
            val hasPaired = group.exists(_.paired)
            def isStray(r: SamRecord): Boolean = hasPaired && !r.paired
            val template = if (name == "*") None else Try(Template(group.iterator.filterNot(isStray))) match {
              case Success(t) =>
                if (!allowMissingMates) {
                  // Check every record, as the group may have secondary or supplementary records but no primary records
                  t.allReads.find(r => r.paired && (if (r.firstOfPair) t.r2 else t.r1).isEmpty).foreach { r =>
                    fail(s"Record ${r.name} is paired but its mate was not found.  Use --allow-missing-mates to allow this.")
                  }
                }
                Some(t)
              case Failure(ex) =>
                if (group.exists(_.paired)) fail(s"Could not group the records with query name $name: ${ex.getMessage}")
                numSharedName += group.length
                None
            }
            alignRecords(group)
            val kept = template match {
              case Some(t) =>
                val keepTemplate = !(filter && group.exists(r => !isStray(r) && hasLeadingIndel(r)))
                if (keepTemplate) numAmbiguous += IndelLeftAligner.fixMateAndSupplementaryInfo(t)
                group.filter(r => if (isStray(r)) !(filter && hasLeadingIndel(r)) else keepTemplate)
              // Unavailable names and unrelated unpaired records are aligned and filtered individually, without repair.
              case None => group.filterNot(r => filter && hasLeadingIndel(r))
            }
            numFiltered += group.length - kept.length
            kept.foreach { rec =>
              if (isTemplateCoordinate && rec.paired && rec.mateMapped && !rec.contains("MC")) fail(
                s"Record ${rec.id} has a mapped mate but no MC tag, required for template-coordinate output.  " +
                  "Add mate CIGARs (e.g. with samtools fixmate -m or fgbio SetMateInformation on query-name grouped input), " +
                  "or query-name sort the input."
              )
              writeProgress.record(rec)
              out.write(rec)
            }
          }
        } finally {
          grouped.safelyClose()
        }
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

  /** Returns a lookup function for the upper-case reference bases spanned by a mapped record, omitting N skips.
    * The returned string concatenates the alignment blocks, including deleted bases but not skipped bases.
    * A supplied FASTA is loaded and checked against every contig in the input's dictionary before returning the
    * lookup function.  Otherwise each lookup reconstructs the bases from the record's MD tag.
    *
    * @param header the input's header
    */
  private[bam] def refBasesFor(header: SAMFileHeader): SamRecord => String = ref match {
    case Some(path) =>
      logger.info("Loading the entire reference FASTA into memory.")
      val dict   = header.getSequenceDictionary
      val seen   = mutable.Set.empty[String]
      val loaded = mutable.Map.empty[String, Array[Byte]] // only the contigs in the input's sequence dictionary
      // An explicit null dictionary avoids discovery of sidecars, whose stale lengths can change FASTA parsing.
      try {
        Using.resource(new FastaSequenceFile(new HtsPath(path.toUri.toString), null, true)) { file =>
          ReferenceSequenceIterator(file).foreach { sequence =>
            val name = sequence.getName
            if (!seen.add(name)) fail(
              s"Contig $name occurs more than once in the reference FASTA (names are compared up to the first " +
                s"whitespace).  Remove the duplicate, or use the reference the input was aligned to: $path"
            )
            // Only the dictionary's own names are used for lookups, not any alternative names (AN)
            if (Option(dict.getSequence(name)).exists(_.getSequenceName == name)) {
              val bases = sequence.getBases
              StringUtil.toUpperCase(bases)
              loaded(name) = bases
            }
          }
        }
      } catch {
        case ex: SAMException => fail(s"Could not read reference FASTA $path: ${ex.getMessage}")
      }
      dict.getSequences.asScala.foreach { sequence =>
        val name   = sequence.getSequenceName
        val length = sequence.getSequenceLength
        val bases  = loaded.getOrElse(name, fail(
          s"Contig $name is not in the reference FASTA.  Use the reference the input was aligned to: $path"
        ))
        if (bases.length != length) fail(
          s"Contig $name has length ${bases.length} in the reference FASTA but $length in the input's sequence " +
            s"dictionary.  Use the reference the input was aligned to: $path"
        )
      }
      rec => {
        val bases = loaded(rec.refName)
        val circular = header.getSequence(rec.refIndex).getAttribute("TP") == "circular"
        if (rec.start > bases.length || (!circular && rec.end > bases.length)) fail(
          s"Record ${rec.name} ends at ${rec.end}, past the end of contig ${rec.refName} (length ${bases.length})."
        )
        val target = new StringBuilder
        var offset = rec.start - 1
        rec.cigar.foreach { elem =>
          if (elem.operator != CigarOperator.N) {
            var consumed = 0
            while (consumed < elem.lengthOnTarget) {
              val pos = if (circular) ((offset.toLong + consumed) % bases.length).toInt else offset + consumed
              val count = math.min(elem.lengthOnTarget - consumed, bases.length - pos)
              target.append(StringUtil.bytesToString(bases, pos, count))
              consumed += count
            }
          }
          offset = if (circular) ((offset.toLong + elem.lengthOnTarget) % bases.length).toInt else offset + elem.lengthOnTarget
        }
        target.result()
      }
    case None =>
      logger.info("Reference FASTA not provided: using the MD tag to retrieve reference bases.")
      rec => {
        if (!rec.contains("MD")) fail(
          s"Record ${rec.name} has no MD tag.  Provide --ref, or add MD tags to the input (e.g. with samtools calmd)."
        )
        val bases = try rec.refBasesFromMd(includeSkippedBases=false) catch {
          case ex: IllegalArgumentException => fail(s"${ex.getMessage}.  Provide --ref, or fix the MD tags (e.g. with samtools calmd).")
        }
        StringUtil.toUpperCase(bases) // the bases are newly allocated, so may be changed in place
        StringUtil.bytesToString(bases)
      }
  }
}

/** Thrown when a record cannot be left-aligned because its data is inconsistent (e.g. with the reference). */
private[bam] class InvalidRecordException(message: String) extends RuntimeException(message)

/** The result of left-aligning a record.
  *
  * @param changed true if the record's alignment (CIGAR or start) was changed
  * @param trimmed true if an indel was trimmed from the start of the record's alignment
  * @param nm the edit distance of the changed alignment (see [[IndelLeftAligner.isMatch]]), or -1 if not changed
  */
private[bam] case class LeftAlignResult(changed: Boolean, trimmed: Boolean, nm: Int = -1)

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
  *                          soft-clipping inserted bases and moving the start past deleted bases, provided aligned
  *                          bases follow it before the first reference skip
  * @param samtoolsMismatches true to count a read and reference base that are both `N` as a mismatch, as samtools
  *                           does, otherwise as a match, as HTSJDK does (see [[IndelLeftAligner.isMatch]])
  */
private[bam] class IndelLeftAligner(useEqualsAndX: Boolean, trimLeadingIndels: Boolean, samtoolsMismatches: Boolean = false) {
  import IndelLeftAligner._

  /** Left-aligns the given [[com.fulcrumgenomics.bam.api.SamRecord]] given the reference bases spanned by the record
    * (including deleted bases, but omitting skipped bases: the blocks between N operators are concatenated).
    * The reference bases are only evaluated if the record may change.  Unmapped records, records without bases, and records with padding (P)
    * operators are not changed.  Read bases stored as `=` are compared using their reference base, and are written as
    * their reference base if the alignment changes.  A record whose CIGAR uses `=` and `X` operators is not changed
    * unless an indel moves, adjacent indel operators need coalescing, or `useEqualsAndX` is true and its CIGAR also
    * uses `M` operators.
    *
    * @throws InvalidRecordException if the record's CIGAR, bases, and reference bases are inconsistent
    */
  def align(rec: SamRecord, refBases: => String): LeftAlignResult = {
    val originalCigar = rec.cigar
    val usesEqualsAndX = originalCigar.exists(e => e.operator == CigarOperator.EQ || e.operator == CigarOperator.X)
    val emitEqAndX     = useEqualsAndX || usesEqualsAndX
    val convertM      = useEqualsAndX && originalCigar.exists(_.operator == CigarOperator.M)
    val mayChange      = originalCigar.exists(_.operator.isIndel) || convertM

    if (rec.unmapped || !mayChange || rec.length == 0 || originalCigar.exists(_.operator == CigarOperator.P)) Unchanged else {
      val target        = refBases
      val originalBases = rec.basesString
      val targetLength  = originalCigar.iterator.filterNot(_.operator == CigarOperator.N).map(_.lengthOnTarget).sum
      if (originalCigar.lengthOnQuery != originalBases.length) throw new InvalidRecordException(
        s"Query length mismatch for record ${rec.name}: ${originalCigar.lengthOnQuery} (CIGAR) vs. ${originalBases.length} (bases)")
      if (targetLength != target.length) throw new InvalidRecordException(
        s"Reference length mismatch for record ${rec.name}: $targetLength (CIGAR excluding skips) vs. ${target.length} " +
          "(reference bases).  Check that the reference FASTA or MD tag matches the alignment.")
      if (target.contains(IndelPaddedAlignment.Indel)) throw new InvalidRecordException(
        s"The reference bases for record ${rec.name} contain the gap character '${IndelPaddedAlignment.Indel}': $target")

      // Clipping may only be the first and/or last elements, and indels are not shifted into clipped bases
      val (leadingClips, inner, trailingClips) = splitClipping(originalCigar)

      // Left-align each block of the alignment separated by reference skips (N), keeping the reference skips in place
      val bases       = resolveMatchMarkers(rec, originalBases, target)
      val newElems    = IndexedSeq.newBuilder[CigarElem]
      var queryOffset = leadingClips.map(_.lengthOnQuery).sum
      var refOffset   = 0
      var blockStart  = 0
      forloop(from=0, until=inner.length + 1) { i =>
        if (i == inner.length || inner(i).operator == CigarOperator.N) {
          if (blockStart < i) {
            val block   = Cigar(inner.slice(blockStart, i))
            val query   = bases.substring(queryOffset, queryOffset + block.lengthOnQuery)
            val ref     = target.substring(refOffset, refOffset + block.lengthOnTarget)
            val aligned = IndelPaddedAlignment(block, query, ref).leftAligned
            newElems ++= aligned.cigar(useEqualsAndX=emitEqAndX, samtoolsMismatches=samtoolsMismatches).elems
            queryOffset += block.lengthOnQuery
            refOffset   += block.lengthOnTarget
          }
          if (i < inner.length) newElems += inner(i)
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

      // Keep an unchanged path's notation unless indels need coalescing or M operators are converted to = and X.
      val samePath = newCigar == originalCigar || withoutEqualsAndX(newCigar) == withoutEqualsAndX(originalCigar)
      val coalesceIndels = originalCigar.iterator.zip(originalCigar.iterator.drop(1)).exists { case (a, b) =>
        a.operator.isIndel && a.operator == b.operator
      }
      val keep = newCigar == originalCigar || (!trim && samePath && !convertM && !coalesceIndels)
      if (keep) Unchanged else {
        if (bases != originalBases && Seq("MM", "Mm", "ML", "Ml", "MN").exists(rec.contains)) throw new InvalidRecordException(
          s"Cannot materialize '=' read bases with base-modification tags for record ${rec.name}.")
        // Resolving '=' must not introduce a base that SAM/BAM storage normalizes (e.g. U to T).
        if (bases != originalBases && !bases.forall("=ACMGRSVTWYHKDBN".contains(_))) throw new InvalidRecordException(
          s"Cannot store the resolved '=' read bases unchanged for record ${rec.name}: $bases")
        rec.cigar = newCigar
        val sequence = rec.header.getSequence(rec.refIndex)
        if (numTrimmedDeleted > 0 && sequence.getAttribute("TP") == "circular") {
          rec.start = ((rec.start.toLong - 1 + numTrimmedDeleted) % sequence.getSequenceLength + 1).toInt
        }
        else rec.start += numTrimmedDeleted
        // An `=` read base may no longer be aligned to its reference base, so write the resolved bases
        if (bases != originalBases) rec.bases = bases
        val tags = AlignmentTags(newCigar, bases, rec.quals, target, targetOffset=numTrimmedDeleted, samtoolsMismatches=samtoolsMismatches)
        if (rec.contains("MD")) rec("MD") = tags.md
        if (rec.contains("NM")) rec("NM") = tags.nm
        if (rec.contains("UQ")) rec("UQ") = tags.uq.map(Int.box).orNull
        if (numTrimmedDeleted > 0 || !samePath) Seq("cs", "ds", "de", "BQ").foreach(rec.remove)
        else if (coalesceIndels) rec.remove("de")
        LeftAlignResult(changed=true, trimmed=trim, nm=tags.nm)
      }
    }
  }
}

private[bam] object IndelLeftAligner {
  /** The result for a record that was not changed. */
  private val Unchanged = LeftAlignResult(changed=false, trimmed=false)

  /** True if the read and reference bases match.  By default this is HTSJDK's `SequenceUtil.basesEqual` (the same
    * base ignoring case, including `N` and other IUPAC codes), as used by other fgbio tools.  If `samtoolsMismatches`
    * is true, they must be the same base ignoring case and not `N`, as in samtools.
    */
  def isMatch(readBase: Char, refBase: Char, samtoolsMismatches: Boolean): Boolean = {
    if (samtoolsMismatches) {
      val base = Character.toUpperCase(readBase)
      base == Character.toUpperCase(refBase) && base != 'N'
    }
    else SequenceUtil.basesEqual(readBase.toByte, refBase.toByte)
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
    * @param rec the read with its original alignment
    * @param bases the read's bases
    * @param target the reference bases spanned by the read, including deleted bases but omitting N skips
    */
  def resolveMatchMarkers(rec: SamRecord, bases: String, target: String): String = {
    if (!bases.contains('=')) bases else {
      val chars = bases.toCharArray
      var queryOffset = 0
      var refOffset   = 0
      rec.cigar.foreach { elem =>
        if (elem.operator.isAlignment) forloop(from=0, until=elem.length) { i =>
          if (chars(queryOffset + i) == '=') chars(queryOffset + i) = target.charAt(refOffset + i)
        }
        queryOffset += elem.lengthOnQuery
        if (elem.operator != CigarOperator.N) refOffset += elem.lengthOnTarget
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
      * @param target reference bases that include those spanned by the alignment, omitting N skips
      * @param targetOffset the offset in `target` of the first reference base of the alignment
      * @param samtoolsMismatches true to count bases that are both `N` as a mismatch (see [[isMatch]])
      */
    def apply(cigar: Cigar, bases: String, quals: Array[Byte], target: String, targetOffset: Int,
              samtoolsMismatches: Boolean): AlignmentTags = {
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
              if (isMatch(bases.charAt(queryOffset + i), refBase, samtoolsMismatches)) matches += 1 else {
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
        if (elem.operator != CigarOperator.N) refOffset += elem.lengthOnTarget
      }
      md.append(matches)
      AlignmentTags(md=md.toString, nm=nm, uq=if (quals.isEmpty) None else Some(uq))
    }
  }

  /** The original alignment and the changed alignment's edit distance, used when repairing the template. */
  private case class OriginalAlignment(start: Int, cigar: Cigar, nm: Int)

  private object OriginalAlignmentKey

  /** Retains the original start/CIGAR and the changed alignment's edit distance without writing a SAM tag.
    * This must be called after query-name grouping, since transient attributes do not survive a disk sort (nor are
    * they written to the output). */
  def setOriginalAlignment(rec: SamRecord, start: Int, cigar: Cigar, nm: Int): Unit = {
    val _ = rec.transientAttrs.update(OriginalAlignmentKey, OriginalAlignment(start, cigar, nm))
  }

  /** Returns the original start/CIGAR, or the current alignment if it was not changed. */
  private def originalAlignment(rec: SamRecord): (Int, Cigar) = rec.transientAttrs.get[OriginalAlignment](OriginalAlignmentKey) match {
    case Some(alignment) => (alignment.start, alignment.cigar)
    case None            => (rec.start, rec.cigar)
  }

  /** Returns the changed alignment's edit distance, even if the record has no NM tag, or None if unchanged. */
  private def changedNm(rec: SamRecord): Option[Int] = rec.transientAttrs.get[OriginalAlignment](OriginalAlignmentKey).map(_.nm)

  /** Returns the CIGAR with hard clipping replaced by soft clipping, and = and X operators replaced by M operators, as
    * SA tags may list an alignment with either. */
  private def normalized(cigar: Cigar): Cigar = withoutEqualsAndX(Cigar(cigar.elems.map { e =>
    if (e.operator == CigarOperator.H) e.copy(operator=CigarOperator.S) else e
  }))

  /** Minimap2 format.c writes end clips, min(query span, reference span) M, then the net I or D.  Reference skips
    * contribute to its reference span; padding and internal clipping do not describe a minimap2 alignment. */
  private def minimapSummary(cigar: Cigar): Option[Cigar] = {
    val (leading, inner, trailing) = splitClipping(normalized(cigar))
    if (inner.isEmpty || inner.exists(e => !e.operator.isAlignment && !e.operator.isIndel && e.operator != CigarOperator.N)) None
    else {
      val query   = inner.map(_.lengthOnQuery).sum
      val target  = inner.map(_.lengthOnTarget).sum
      val matches = math.min(query, target)
      val elems   = IndexedSeq.newBuilder[CigarElem]
      elems ++= leading
      if (matches > 0) elems += CigarElem(CigarOperator.M, matches)
      if (query > target) elems += CigarElem(CigarOperator.I, query - target)
      if (target > query) elems += CigarElem(CigarOperator.D, target - query)
      elems ++= trailing
      Some(Cigar(elems.result()))
    }
  }

  private def isR1(rec: SamRecord): Boolean = !rec.paired || rec.firstOfPair

  private def strand(rec: SamRecord): String = if (rec.positiveStrand) "+" else "-"

  /** Fixes the mate information on the primary and supplementary records of a template, sets the mate start and mate
    * CIGAR (MC) on secondary records from their mate's primary record, and updates the positions, CIGARs, and edit
    * distances in the SA tags of all records from the primary and supplementary records they refer to.
    *
    * @param template the template to fix
    * @return the number of SA entries that did not match exactly one alignment and were left unchanged
    */
  def fixMateAndSupplementaryInfo(template: Template): Int = {
    // fixMateInfo sets ms from the mate's AS, but existing ms tags may be fixmate or minimap2 scores.
    val scores = template.allReads.map(r => (r, r.get[Any]("ms"))).toIndexedSeq
    template.fixMateInfo()
    scores.foreach { case (rec, score) => rec("ms") = score.orNull }

    // Update the mate start, insert size, and mate cigar on secondary records from the mate's primary record.  The
    // secondary's mate start (and mate cigar, if present) must match the mate's primary alignment before it was changed.
    def fixSecondaries(secondaries: Seq[SamRecord], matePrimary: Option[SamRecord]): Unit = matePrimary.filter(_.mapped).foreach { mate =>
      val mateCigar                      = mate.cigar.toString
      val (originalStart, originalCigar) = originalAlignment(mate)
      secondaries.iterator.filter { r =>
        r.paired && r.mateMapped && r.mateRefIndex == mate.refIndex && r.mateStart == originalStart &&
          r.mateNegativeStrand == mate.negativeStrand &&
          r.mateCigar.forall(_ == originalCigar)
      }.foreach { r =>
        r.mateStart  = mate.start
        r.insertSize = SamPairUtil.computeInsertSize(r.asSam, mate.asSam)
        r("MC")      = mateCigar
      }
    }
    fixSecondaries(template.r1Secondaries, template.r2)
    fixSecondaries(template.r2Secondaries, template.r1)

    // Update the SA tags, which list the other primary or supplementary alignments for the same read
    var unmatched = 0
    val records = template.allReads.filter(_.contains("SA")).toIndexedSeq
    if (records.nonEmpty) {
      val alignments = (template.primaryReads ++ template.r1Supplementals ++ template.r2Supplementals).filter(_.mapped).toIndexedSeq
      records.foreach { rec =>
        val entries = rec[String]("SA").split(';').filter(_.nonEmpty).map { entry =>
          updateSaEntry(entry, rec, alignments).getOrElse { unmatched += 1; entry }
        }
        rec("SA") = entries.mkString("", ";", ";")
      }
    }
    unmatched
  }

  /** Returns the SA tag entry with its position, CIGAR, and edit distance (NM) replaced by those of the alignment it
    * refers to, or None if the entry cannot be parsed or does not match exactly one alignment.  The alignment must be for
    * the same read as `rec` (but not `rec` itself), on the same reference and strand, and its start and CIGAR before it
    * was changed (see [[setOriginalAlignment]]) must match the entry's, ignoring the difference between hard and soft
    * clipping, and between M and = or X operators.  A minimap2-style summary of the original CIGAR is also accepted
    * (see [[minimapSummary]]), but exact and summarized candidates together must identify exactly one alignment.
    * The entry's clipping operators are kept (e.g. soft clipping for a
    * hard-clipped alignment) where its clipping at that end has no hard clipping.  The edit distance is that computed for
    * the alignment when it was changed (whether or not the alignment has an NM tag), otherwise the entry's is kept.
    *
    * @param entry a single entry from the SA tag of `rec`
    * @param rec the record with the SA tag
    * @param alignments the primary and supplementary alignments of the template
    */
  private def updateSaEntry(entry: String, rec: SamRecord, alignments: Seq[SamRecord]): Option[String] = entry.split(',') match {
    case Array(refName, pos, strandString, cigarString, mapq, nm) =>
      Try((pos.toInt, Cigar(cigarString))).toOption.flatMap { case (entryStart, entryCigar) =>
        val normalizedEntry = normalized(entryCigar)
        val matching = alignments.filter { r =>
          !(r eq rec) && isR1(r) == isR1(rec) && r.refName == refName && strand(r) == strandString && {
            val (start, cigar) = originalAlignment(r)
            start == entryStart && (normalized(cigar) == normalizedEntry || minimapSummary(cigar).contains(normalizedEntry))
          }
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
            val newNm = changedNm(alignment).map(_.toString).getOrElse(nm)
            Some(Seq(refName, alignment.start.toString, strandString, cigar.toString, mapq, newNm).mkString(","))
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
  def sameBase(a: Char, b: Char): Boolean = Character.toUpperCase(a) == Character.toUpperCase(b)

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
  def cigar(useEqualsAndX: Boolean = false, samtoolsMismatches: Boolean = false): Cigar = {
    def operatorFrom(index: Int): CigarOperator = {
      if (target(index) == Indel) CigarOperator.I
      else if (query(index) == Indel) CigarOperator.D
      else if (!useEqualsAndX) CigarOperator.M
      else if (IndelLeftAligner.isMatch(query(index), target(index), samtoolsMismatches)) CigarOperator.EQ
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

    // Returns the index after the last column of the run of indel columns that contains `index`
    def runEnd(index: Int): Int = {
      var end = index
      while (end < length && isIndel(end)) end += 1
      end
    }

    // Returns the index of the first column of the run of indel columns that contains `index`
    def runStart(index: Int): Int = {
      var start = index
      while (start > 0 && isIndel(start - 1)) start -= 1
      start
    }

    // Places the insertions of the run of indel columns [start, end) before its deletions
    def reorder(start: Int, end: Int): Unit = {
      val inserted = (start until end).filter(i => targetArr(i) == Indel).map(queryArr)
      val deleted  = (start until end).filter(i => queryArr(i) == Indel).map(targetArr)
      inserted.indices.foreach { i => queryArr(start + i) = inserted(i); targetArr(start + i) = Indel }
      deleted.indices.foreach  { i => queryArr(start + inserted.length + i) = Indel; targetArr(start + inserted.length + i) = deleted(i) }
    }

    // Shifts the insertions or deletions of the run of indel columns [start, end) left by one column (the run must be
    // ordered with insertions first), returning true if any were shifted.  Only the columns [start - 1, end) change.
    def shift(start: Int, end: Int): Boolean = {
      val numInserted = (start until end).count(i => targetArr(i) == Indel)
      val numDeleted  = end - start - numInserted
      val prev        = start - 1 // the aligned column before the run, if any
      val inserted    = queryArr.slice(start, start + numInserted)
      val deleted     = targetArr.slice(start + numInserted, end)
      if (prev >= 0 && numInserted > 0 && sameBase(queryArr(prev), inserted.last)) {
        // [prev] [inserted] [deleted] => [prev query base +: inserted.init] [inserted.last / prev target base] [deleted]
        val prevQuery  = queryArr(prev)
        val prevTarget = targetArr(prev)
        forloop(from=0, until=numInserted) { i =>
          queryArr(prev + i)  = if (i == 0) prevQuery else inserted(i - 1)
          targetArr(prev + i) = Indel
        }
        queryArr(prev + numInserted)  = inserted.last
        targetArr(prev + numInserted) = prevTarget
        true
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
        true
      }
      else false
    }

    // In a pure run the non-gap sequence does not change; only the gap sequence's two boundary columns move.
    def shiftPure(start: Int, end: Int, insertion: Boolean): Boolean = {
      val bases = if (insertion) queryArr else targetArr
      val gaps  = if (insertion) targetArr else queryArr
      if (start == 0 || !sameBase(bases(start - 1), bases(end - 1))) false else {
        gaps(end - 1) = gaps(start - 1)
        gaps(start - 1) = Indel
        true
      }
    }

    // Each pass orders, then shifts, every run by at most one column, until no run can be shifted.  Since each shift
    // moves indel columns left (and ordering a run does not move any indel column), this ends.  Whether a run can be
    // shifted depends only on its columns and the aligned column before it, and shifting it changes only those columns,
    // so a run that was not shifted can only be shifted in a later pass once one of those columns is changed by shifting
    // another run (e.g. one that joins it).  Each pass therefore only visits the runs that contain a column changed in
    // the previous pass, which gives the same result as visiting every run, without re-scanning the whole alignment.
    var changed  = IndexedSeq((0, length)) // the ranges of columns changed in the previous pass, in order
    var pureRuns = Map.empty[Int, (Int, Boolean)] // start -> (end, isInsertion), avoiding rescans after pure shifts
    while (changed.nonEmpty) {
      val nextChanged  = IndexedSeq.newBuilder[(Int, Int)]
      val nextPureRuns = Map.newBuilder[Int, (Int, Boolean)]
      var visitedEnd   = 0 // the end of the last run visited in this pass
      changed.foreach { case (from, until) =>
        var index = math.max(from, visitedEnd)
        while (index < until) {
          if (!isIndel(index)) index += 1 else {
            // A neighboring run may have joined this one since the previous pass; then rediscover the whole run.
            val cached = pureRuns.get(index).filter { case (end, _) =>
              (index == 0 || !isIndel(index - 1)) && (end == length || !isIndel(end))
            }
            val start = if (cached.isDefined) index else runStart(index)
            val end   = cached.map(_._1).getOrElse(runEnd(index))
            val pureInsertion = cached.map(_._2).orElse {
              if ((start until end).forall(i => targetArr(i) == Indel)) Some(true)
              else if ((start until end).forall(i => queryArr(i) == Indel)) Some(false)
              else None
            }
            val shifted = pureInsertion match {
              case Some(insertion) =>
                val moved = shiftPure(start, end, insertion)
                if (moved) nextPureRuns += ((start - 1, (end - 1, insertion)))
                moved
              case None =>
                reorder(start, end)
                shift(start, end)
            }
            if (shifted) nextChanged += ((start - 1, end))
            visitedEnd = end
            index      = end
          }
        }
      }
      changed  = nextChanged.result()
      pureRuns = nextPureRuns.result()
    }

    IndelPaddedAlignment(new String(queryArr), new String(targetArr))
  }
}
