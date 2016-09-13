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
 *
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Metric, NumericCounter, ProgressLogger}
import dagr.commons.CommonsDef._
import dagr.commons.io.{Io, PathUtil}
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.SequenceUtil
import htsjdk.samtools.util.SequenceUtil.{isValidBase => validBase}
import htsjdk.samtools.{CigarOperator, SAMRecord, SamReaderFactory}

@clp(description =
  """
    |Calculates the error rate by read position on mapped BAMs.
    |
    |If a reference is given, the sort order must be coordinate.  If no reference is given, the MD tag must be present
    |for every mapped read and computed correctly, and no sort order is required.
    |
    |The output metrics text file will have the extension 'error_rate_by_read_position.txt'.  The output file will
    |contain a row per logical read and per position in that read.  The 'read_number' will be 0 for fragments, 1 for the
    |first end of a pair, and 2 for the second end of a pair. Both the error rate and number of bases observed will be
    |provided.  Soft-clipped bases, and unmapped, secondary, supplementary, and vendor failing reads are ignored.
    |Optionally, duplicates are ignored.
    |
    |Mismatches and optionally insertions are considered errors, while deletions and non-ACGT bases are ignored.
  """,
  group = ClpGroups.SamOrBam)
class ErrorRateByReadPosition
( @arg(flag="i", doc="Input BAM file.") val input: PathToBam,
  @arg(flag="o", doc="Output metrics prefix. If not given, will use the input BAM basename as a prefix.")
  val output: Option[PathPrefix] = None,
  @arg(flag="r", doc="Reference fasta. If not given, will use the MD tag.") val ref: Option[PathToFasta] = None,
  @arg(flag="d", doc="Include duplicate reads, otherwise ignore.") val includeDuplicates: Boolean = false,
  @arg(flag="I", doc="Count insertions as errors, otherwise as matches.") val includeInsertions: Boolean = false,
  @arg(flag="m", doc="The minimum mapping quality for a read to be included.") val minMappingQuality: Int = 20,
  @arg(flag="q", doc="The minimum base quality for a base to be included.") val minBaseQuality: Int = 0
) extends FgBioTool with LazyLogging {
  import scala.collection.JavaConversions.asScalaIterator

  Io.assertReadable(input)
  ref.foreach(Io.assertReadable)
  output.foreach(out => Io.assertCanWriteFile(out))

  def execute(): Unit = {
    val progress = new ProgressLogger(logger, verb="read", unit=1e6.toInt)
    val in       = SamReaderFactory.make.open(input)
    val iterator = in.iterator()
    val counter  = ErrorRateByReadPositionCounter(ref=ref, includeInsertions=includeInsertions, minBaseQuality=minBaseQuality)

    if (ref.exists(r => in.getFileHeader.getSortOrder != SortOrder.coordinate)) {
      fail("Input BAM must be coordinate sorted when the reference is given.")
    }

    iterator.filterNot { rec =>
      rec.getReadUnmappedFlag || rec.isSecondaryOrSupplementary ||
        (!includeDuplicates && rec.getDuplicateReadFlag) || rec.getReadFailsVendorQualityCheckFlag ||
        (rec.getMappingQuality < minMappingQuality)
    }
    .foreach { rec =>
      counter.accept(rec)
      progress.record(rec)
    }

    val path = PathUtil.pathTo(output.getOrElse(PathUtil.removeExtension(input)) + ErrorRateByReadPositionMetric.FileExtension)
    Metric.write(metrics=counter.metrics, path=path)

    in.safelyClose()
  }
}

object ErrorRateByReadPositionMetric {
  val FileExtension = ".error_rate_by_read_position.txt"
}

/** Gives the error rate and total number of bases observed for each position within the read.
  *
  * @param read_number the read number (ex. 0 for fragments, 1 for read one of a pair, 2 for read two of a pair)
  * @param position the read position (1-based).
  * @param error_rate the error rate at the given position n the read
  * @param bases_total the total number of bases observed at this position.
  */
case class ErrorRateByReadPositionMetric
( read_number: Int,
  position: Int,
  error_rate: Double,
  bases_total: Long
) extends Metric

private object ErrorRateByReadPositionCounter {
  def apply(ref: Option[PathToFasta], includeInsertions: Boolean, minBaseQuality: Int): ErrorRateByReadPositionCounter = {
    ref match {
      case Some(r) => new ErrorRateByReadPositionWithReferenceCounter(ref=r, includeInsertions=includeInsertions, minBaseQuality=minBaseQuality)
      case None    => new ErrorRateByReadPositionWithMdCounter(includeInsertions=includeInsertions, minBaseQuality=minBaseQuality)
    }
  }
}

/** Base class for computing the error rate by position.
  *
  * Sub-classes should implement the [updateErrors] method to extract the read bases and qualities, and reference bases.
  * They should then call [incrementIfValid] or [increment] to update the counters for each read offset.
  *
  * @param minBaseQuality the minimum base quality to allow.
  */
private abstract class ErrorRateByReadPositionCounter(includeInsertions: Boolean, minBaseQuality: Int) {
  private val errorsByReadOffset = List.tabulate(3)(x => new NumericCounter[Int]())
  private val totalByReadOffset  = List.tabulate(3)(x => new NumericCounter[Int]())

  /** Gets the list of metrics sorted by read number and then position. */
  def metrics: Seq[ErrorRateByReadPositionMetric] = {
    if (totalByReadOffset.forall(_.isEmpty)) return Seq.empty
    // create the metrics
    errorsByReadOffset.zip(totalByReadOffset).zipWithIndex.filter {
      case ((errCounter, totalCounter), idx) => totalCounter.nonEmpty
    }.flatMap { case ((errCounter, totalCounter), idx) =>
      val maxOffset = totalCounter.map(_._1).max
      Range.inclusive(0, maxOffset).map { offset =>
        val error_bases = errCounter.countOf(offset)
        val error_total = totalCounter.countOf(offset)
        val error_rate  = if (0 == error_total) 0 else error_bases / error_total.toDouble
        new ErrorRateByReadPositionMetric(
          read_number = idx,
          position    = offset + 1,
          error_rate  = error_rate,
          bases_total = error_total
        )
      }
    }.sortBy(m => (m.read_number, m.position))
  }

  /** Checks that the given bases are valid, and increments the # of errors and total bases. The read index is
    * 0 for fragment reads, 1 for the first of a pair, and 2 for second of a pair.  The read position is zero-based. */
  protected def incrementIfValidBase(readIndex: Int, readOffset: Int, readBase: Byte, qual: Byte, refBase: Byte): Unit = {
    // Assert the bases must be valid
    if (validBase(refBase) && validBase(readBase)) {
      increment(readIndex=readIndex, readOffset=readOffset, qual=qual, isError = !SequenceUtil.basesEqual(readBase, refBase))
    }
  }

  /** Increments the error count if `isError` is true. Also increments the total bases observed. The read index is
    * 0 for fragment reads, 1 for the first of a pair, and 2 for second of a pair.  The read offset is zero-based. */
  protected def increment(readIndex: Int, readOffset: Int, qual: Byte, isError: Boolean = true): Unit = {
    if (minBaseQuality <= qual) {
      if (isError) errorsByReadOffset(readIndex).count(readOffset)
      totalByReadOffset(readIndex).count(readOffset)
    }
  }

  /** Implement me! */
  protected def updateErrors(rec: SAMRecord, readIndex: Int): Unit

  def accept(rec: SAMRecord): Unit = {
    val readIndex: Int = if (!rec.getReadPairedFlag) 0 else if (rec.getFirstOfPairFlag) 1 else 2
    updateErrors(rec=rec, readIndex=readIndex)
  }

  protected def maybeReverseComplement(bases: Array[Byte], rc: Boolean): Array[Byte] = {
    if (rc) SequenceUtil.reverseComplement(bases)
    bases
  }

  protected def maybeReverse(qualities: Array[Byte], rc: Boolean): Array[Byte] = {
    if (rc) SequenceUtil.reverse(qualities, 0, qualities.length)
    qualities
  }
}

/** Computes the error rate using a reference FASTA.  Does not rely on the MD tag. */
private class ErrorRateByReadPositionWithReferenceCounter
( ref: PathToFasta,
  includeInsertions: Boolean,
  minBaseQuality:Int
) extends ErrorRateByReadPositionCounter(includeInsertions=includeInsertions, minBaseQuality=minBaseQuality) {
  import scala.collection.JavaConversions.iterableAsScalaIterable

  val referenceFile     = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
  var referenceSequence = referenceFile.getSequence(referenceFile.getSequenceDictionary.getSequence(0).getSequenceName)

  protected def updateErrors(rec: SAMRecord, readIndex: Int): Unit = {
    // Move to the correct reference.
    if (rec.getReferenceIndex != referenceSequence.getContigIndex) {
      referenceSequence = referenceFile.getSequence(referenceFile.getSequenceDictionary.getSequence(rec.getReferenceIndex).getSequenceName)
    }

    // Get the reference bases
    val referenceBases = referenceSequence.getBases

    // We proceed through the read in the order it was aligned, but if the read is on the reverse strand, we will
    // need to invert the read position since the last base in the read will be the first base sequenced
    val readOffsetIncrement = if (rec.getReadNegativeStrandFlag) -1 else 1
    val readLength = rec.getReadLength

    // Get the read bases and qualities.
    val read  = rec.getReadBases
    val quals = rec.getBaseQualities

    // For each alignment block
    for (block <- rec.getAlignmentBlocks) {
      var refPosition  = block.getReferenceStart - 1
      var readPosition = block.getReadStart - 1
      val readEnd      = readPosition + block.getLength
      var readOffset   = if (rec.getReadNegativeStrandFlag) readLength - readPosition - 1 else readPosition

      // Matches/mismatches
      while (readPosition < readEnd) {
        val readBase = read(readPosition)
        val qual     = quals(readPosition)
        val refBase  = referenceBases(refPosition)
        incrementIfValidBase(readIndex=readIndex, readOffset=readOffset, readBase=readBase, qual=qual, refBase=refBase)
        refPosition  += 1
        readPosition += 1
        readOffset   += readOffsetIncrement
      }
    }

    // Handle insertions
    if (includeInsertions) {
      var readOffset = if (rec.getReadNegativeStrandFlag) rec.getReadLength - 1 else 0
      rec.getCigar.getCigarElements.foreach { elem =>
        val op: CigarOperator = elem.getOperator
        if (op == CigarOperator.INSERTION) {
          var i = 0
          while (i < elem.getLength) {
            val qual = quals(readOffset)
            increment(readIndex = readIndex, readOffset = readOffset, qual=qual)
            readOffset += readOffsetIncrement
            i += 1
          }
        }
        else if (op.consumesReadBases) {
          readOffset += readOffsetIncrement * elem.getLength
        }
      }
    }
  }
}

/** Computes the error rate using the MD tag. Does not rely on a reference FASTA. */
private class ErrorRateByReadPositionWithMdCounter(includeInsertions: Boolean, minBaseQuality: Int)
  extends ErrorRateByReadPositionCounter(includeInsertions=includeInsertions, minBaseQuality=minBaseQuality) {

  protected def updateErrors(rec: SAMRecord, readIndex: Int): Unit  = {
    // Use the MD tag!
    val read              = maybeReverseComplement(rec.getReadBases.clone(), rc=rec.getReadNegativeStrandFlag)
    val quals             = maybeReverse(rec.getBaseQualities.clone(), rc=rec.getReadNegativeStrandFlag)
    val reference         = maybeReverseComplement(SequenceUtil.makeReferenceFromAlignment(rec, false), rc=rec.getReadNegativeStrandFlag)
    var readPosition      = 0

    while (readPosition < read.length) {
      val readBase = read(readPosition)
      val qual     = quals(readPosition)
      reference(readPosition) match {
        case '-'                       => // insertion
          if (includeInsertions) this.increment(readIndex=readIndex, readOffset=readPosition, qual=qual) // insertion
        case refBase if '0' != refBase => // not a soft-clip
          this.incrementIfValidBase(readIndex=readIndex, readOffset=readPosition, readBase=readBase, qual=qual, refBase=refBase)
        case refBase                   => // ignore
      }
      readPosition += 1
    }
  }
}
