/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics
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

import com.fulcrumgenomics.FgBioDef.{PathToFasta, SafelyClosable, _}
import com.fulcrumgenomics.alignment.{Cigar, CigarElem}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef.PathToBam
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sequences}
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.SamPairUtil.SetMateInfoIterator
import htsjdk.samtools.reference.ReferenceSequenceFileWalker

import scala.annotation.tailrec
import scala.collection.mutable.ArrayBuffer

@clp(description =
  """
    |Left-aligns indels in a SAM/BAM file.
    |
    |This tool is preferred over GATK's LeftAlignIndels, as the latter can only left-align reads with at most one indel.
    |
    |Mate information will be reset.
  """,
  group = ClpGroups.SamOrBam)
class LeftAlignIndels
( @arg(flag='i', doc="Input BAM file.") input: PathToBam,
  @arg(flag='o', doc="Output BAM file.") output: PathToBam,
  @arg(flag='r', doc="Reference sequence FASTA file, otherwise the MD tag is required.") val ref: Option[PathToFasta] = None,
  @arg(flag='x', doc="If specified, do not fail when reads marked as paired are missing their mate pairs.")
  val allowMissingMates: Boolean = false,
  @arg(flag='X', doc="Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).") useEqualsAndX: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val reader       = SamSource(input)
    val header       = reader.header
    val walker       = ref.map(r => new ReferenceSequenceFileWalker(r))
    val aligner      = new LeftIndelAligner(useEqualsAndX=useEqualsAndX)
    val readProgress = ProgressLogger(logger=logger, verb="read")
    var leftAligned  = 0L
    val inSamOrder   = SamOrder(header)
    // Check if we need to resort the read into query-name/query-grouped order to be able to reset mate information
    val resortInput  = !inSamOrder.exists(order => order.sortOrder == SortOrder.queryname || order.groupOrder == GroupOrder.query)
    // We can only use the reference when the input is coordinate order
    val useRef       = ref.isDefined && inSamOrder.contains(SamOrder.Coordinate)
    // If the input is not re-sorted, then we do not need to re-sort again into the original order when writing.
    val writer       = SamWriter(output, header, sort = if (resortInput) SamOrder(header) else None)

    // The method to left-align a record, based on if a walker is available and the input is coordinate sorted
    val alignFunc: SamRecord => Boolean = walker match {
      case Some(_walker) if useRef =>
        logger.info("Input is coordinate sorted and reference FASTA given: using the input reference FASTA to retrieve reference bases.")
        rec => aligner.align(rec, _walker)
      case _walker                 =>
        if (walker.isEmpty) logger.info("Reference FASTA not provided: using the MD tag retrieve reference bases.")
        else logger.warning("Reference FASTA provided but the input is not coordinate sorted: using the MD tag retrieve reference bases.")
        rec => aligner.align(rec, new String(rec.refBasesFromMd(includeDeletedBases=true)))
    }

    // The method to left-align a record, and update counters.
    def process(rec: SamRecord): SamRecord = {
      if (alignFunc(rec)) leftAligned += 1
      readProgress.record(rec)
      rec
    }

    // Left-align the records and return a query-grouped iterator over the SamRecords
    val queryGroupedIterator = {
      // If the input is already sorted into query-name or query-grouped order, just map, otherwise we will need to re-sort.
      if (!resortInput) reader.map(process) else {
        val sorter = Bams.sorter(SamOrder.Queryname, reader.header)
        reader.foreach { rec => sorter += process(rec) }
        sorter
      }
    }

    // Write out the records after re-setting mate information.
    val writeProgress = ProgressLogger(logger=logger, verb="Written")
    val iterator = new SetMateInfoIterator(queryGroupedIterator.iterator.map(_.asSam), true, allowMissingMates).map(_.asInstanceOf[SamRecord])
    iterator.foreach { rec =>
      writeProgress.record(rec)
      writer.write(rec)
    }
    writeProgress.logLast()

    // Clean up
    walker.foreach(_.close())
    reader.safelyClose()
    writer.close()

    val pct = leftAligned * 100.0 / writeProgress.getCount
    logger.info(f"Left aligned $leftAligned%,d out of ${writeProgress.getCount}%,d records ($pct%.2f%%)")
  }
}

class LeftIndelAligner(useEqualsAndX: Boolean) {
  /** Left-aligns the given [[SamRecord]]. */
  def align(rec: SamRecord, walker: ReferenceSequenceFileWalker): Boolean = if (rec.unmapped) false else {
    val refBases = walker.get(rec.refIndex).getBases.slice(from = rec.start-1, until = rec.end).map(_.toChar.toUpper).mkString("")
    this.align(rec, refBases)
  }

  /** Left-aligns the given [[SamRecord]]. */
  def align(rec: SamRecord, refBases: String): Boolean = if (rec.unmapped) false else {
    val original    = IndelPaddedAlignment(rec, refBases)
    val leftAligned = original.leftAligned
    val cigar       = leftAligned.cigar(useEqualsAndX=useEqualsAndX)

    // Set the cigar
    rec.cigar = {
      val leadingHardClip  = rec.cigar.takeWhile(_.operator == CigarOperator.H).map(_.length).sum
      val trailingHardClip = rec.cigar.elems.drop(1).reverse.takeWhile(_.operator == CigarOperator.H).map(_.length).sum
      (leadingHardClip, trailingHardClip) match {
        case (0, 0) => cigar
        case (n, 0) => new Cigar(CigarElem(CigarOperator.H, n) +: cigar.elems)
        case (0, m) => new Cigar(cigar.elems :+ CigarElem(CigarOperator.H, m))
        case (n, m) => new Cigar(CigarElem(CigarOperator.H, n) +: cigar.elems :+ CigarElem(CigarOperator.H, m))
      }
    }

    original == leftAligned
  }
}


object IndelPaddedAlignment {
  /** The character in the read or reference for insertions or deletions respectively. */
  val Indel: Char = '-'
  /** The character in the reference to indicate the read was soft-clipped. */
  val SoftClip: Char = 'X'

  // Convenience values
  private val IndelString: String = Indel.toString
  private val SoftClipString: String = SoftClip.toString

  /** Builds a padded alignment for the given [[SamRecord]] and its corresponding target bases.  The record must be
    * mapped. NB: hard-clipped bases are lost. */
  def apply(rec: SamRecord, target: String): IndelPaddedAlignment = {
    require(rec.mapped)
    this.apply(rec.cigar, rec.basesString, target)
  }

  /** Builds a padded alignment for the given [[Cigar]], query and target target bases. NB: hard-clipped bases are lost. */
  def apply(cigar: Cigar, query: String, target: String): IndelPaddedAlignment = {
    require(cigar.lengthOnQuery == query.length,
      s"query length mismatch: ${cigar.lengthOnQuery} vs. ${query.length}")
    require(cigar.lengthOnTarget == target.length,
      s"target length mismatch: ${cigar.lengthOnTarget} vs. ${target.length}")

    val queryBuilder  = new ArrayBuffer[Char]()
    val targetBuilder = new ArrayBuffer[Char]()

    // Build the padded query
    cigar.foldLeft(0) { case (index, elem) =>
      if (elem.lengthOnQuery > 0) { // M, X, EQ, I, S
        queryBuilder.appendAll(query.slice(index, index + elem.length))
        index + elem.length
      }
      else { // H or D
        if (elem.operator == CigarOperator.D) {
          require(elem.operator == CigarOperator.D)
          queryBuilder.appendAll(IndelString * elem.length)
        }
        index
      }
    }

    // Build the padded target
    cigar.foldLeft(0) { case (index, elem) =>
      // target
      if (elem.lengthOnTarget == 0) { // S, H, I
        if (elem.operator == CigarOperator.I) targetBuilder.appendAll(IndelString * elem.length)
        else if (elem.operator == CigarOperator.S) targetBuilder.appendAll(SoftClipString * elem.length)
      }
      else { // M, X, EQ, D
        targetBuilder.appendAll(target.slice(index, index + elem.length))
      }
      index + elem.lengthOnTarget
    }

    new IndelPaddedAlignment(queryBuilder.mkString, targetBuilder.mkString)
  }}

/** Stores an alignment padded for indels.
  *
  * The query and target strings are padded as follows:
  * - an inserted base (in the query) has a corresponding gap character ('-') in the target
  * - a deleted base (in the target) has a corresponding gap character ('-') in the query
  * - a soft-clipped base (in the query) has a corresponding soft-clip character ('X') in the target
  *
  * For example:
  *
  * query:  AAAGATTACAGATT--ACA
  * target: XXXGA--ACAGATTAAACA
  *
  * In this example, query bases 1-3 ("AAA") are soft-clipped, query bases 6-7 ("TT") are inserted, and target bases
  * 10-11 ("AA") are deleted.
  *
  * */
case class IndelPaddedAlignment private(query: String, target: String){
  import IndelPaddedAlignment.{Indel, SoftClip}
  require(query.length == target.length,
    f"Lengths differ ${query.length} != ${target.length}\n${query.mkString("")}\n${target.mkString("")}")
  require(length > 0)

  def length: Int = query.length

  def cigar(useEqualsAndX: Boolean = false): Cigar = {
    val firstOperator = operatorAt(0, useEqualsAndX=useEqualsAndX)
    _cigar(index = 1, elems = Seq.empty, lastElem = CigarElem(firstOperator, 1), useEqualsAndX=useEqualsAndX)
  }

  /** Builds a left-aligned padded alignment from this alignment. */
  def leftAligned: IndelPaddedAlignment = {
    val queryArr  = this.query.toArray
    val targetArr = this.target.toArray
    var index     = 0

    // skip over leading indels and soft-clipping
    while (query(index) == Indel || targetArr(index) == Indel || targetArr(index) == SoftClip) {
      index += 1
    }

    while (index < queryArr.length) {
      if (queryArr(index) == Indel && index > 0) { // start of a deletion
        // find the span of the deletion
        if (index > 0) require(queryArr(index-1) != Indel)
        var end = index
        while (end < length && queryArr(end) == Indel) {
          end += 1
        }
        val offset = end - index // to skip over the shifted deletion
        // shift over the deletion, if necessary
        if (index < end) {
          end -= 1 // move to the last deleted base
          while (0 < index && Sequences.compatible(queryArr(index-1), targetArr(end)) && targetArr(index-1) != Indel && targetArr(index-1) != SoftClip) {
            queryArr(end) = queryArr(index-1) // or ref(end)
            queryArr(index-1) = Indel
            index -= 1
            end -= 1
          }
        }
        index += offset
      }
      else if (targetArr(index) == Indel) { // start of an insertion
        // find the span of the insertion
        if (index > 0) require(targetArr(index-1) != Indel)
        var end = index
        while (end < length && targetArr(end) == Indel) {
          end += 1
        }
        // shift over the insertion, if necessary
        val offset = end - index // to skip over the shifted insertion
        if (index < end) {
          end -= 1 // move to the last inserted base
          while (0 < index && Sequences.compatible(targetArr(index-1), queryArr(end)) && queryArr(index-1) != Indel) {
            targetArr(end) = targetArr(index-1) // or read(end)
            targetArr(index-1) = Indel
            index -= 1
            end -= 1
          }
        }
        index += offset
      }
      else {
        index += 1
      }
    }

    IndelPaddedAlignment(queryArr.mkString(""), targetArr.mkString(""))
  }

  /** Gets the operator in the alignment at the given index */
  private def operatorAt(index: Int, useEqualsAndX: Boolean): CigarOperator = {
    (target(index), query(index)) match {
      case (Indel, _)          => CigarOperator.I
      case (_, Indel)          => CigarOperator.D
      case (SoftClip, _)       => CigarOperator.S
      case _ if !useEqualsAndX => CigarOperator.M
      case (qBase, tBase)      => if (Sequences.compatible(qBase, tBase)) CigarOperator.EQ else CigarOperator.X
    }
  }

  /** Recursively builds a cigar starting at a given index.  */
  @tailrec
  private def _cigar(index: Int, elems: Seq[CigarElem], lastElem: CigarElem, useEqualsAndX: Boolean): Cigar = {
    if (index == length) Cigar(elems.toIndexedSeq :+ lastElem)
    else {
      val nextOp = operatorAt(index=index, useEqualsAndX=useEqualsAndX)
      if (nextOp == lastElem.operator) {
        _cigar(index + 1, elems, lastElem.copy(length=lastElem.length + 1), useEqualsAndX=useEqualsAndX)
      }
      else {
        _cigar(index + 1, elems :+ lastElem, CigarElem(nextOp, length=1), useEqualsAndX=useEqualsAndX)
      }
    }
  }
}
