/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam.identifyprimers

import com.fulcrumgenomics.alignment.Alignable
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.util.{Metric, Sequences}
import htsjdk.samtools.util.{Locatable, OverlapDetector, SequenceUtil}

private[identifyprimers] object Primer {

  /** Writes the given primers to file. */
  def write(path: FilePath, primers: TraversableOnce[Primer]): Unit = Metric.write(path, primers)

  /** Reads primers from a given path and validates the primers (see [[validatePrimers()]]).
    *
    * NB: the reverse strand primers stored at the path are assumed to have sequence in the primer order (i.e. not the
    * genomic top strand).
    * */
  def read(path: FilePath, multiPrimerPairs: Boolean = false): Seq[Primer] = {
    val parser  = DelimitedDataParser(path, '\t')
    val primers = parser.map { row =>
      val positive_strand = isPositiveStrand(row.get[String]("positive_strand", true).getOrElse(row.apply[String]("strand")))
      Primer(
        pair_id        = row.apply[String]("pair_id"),
        primer_id      = row.apply[String]("primer_id"),
        sequence       = row.apply[String]("sequence").toUpperCase,
        ref_name       = row.apply[String]("ref_name"),
        start          = row.apply[Int]("start"),
        end            = row.apply[Int]("end"),
        positive_strand = positive_strand
      )
    }.toSeq
    validatePrimers(primers, multiPrimerPairs)
    primers
  }

  /** Converts strand to a True if forward, false if reverse. */
  def isPositiveStrand(strand: String): Boolean = strand.toLowerCase match {
    case "+" | "f" | "fwd" | "for" | "forward" | "positive" | "true" => true
    case "-" | "r" | "rev" | "reverse" | "negative" | "false"        => false
    case _ => throw new IllegalArgumentException(s"Could not parse strand: $strand")
  }

  /** Validates a set of primers.
    *
    * If `multiPrimerPairs` is true, validates that there is at least one forward and one reverse primer for each pair.
    * Otherwise, validates that two primers exist for each pair, and that one is forward and the other is reverse.
    *
    * Also validates that no forward primers that have the same refName/start, or reverse primers with the same
    * refName/end.
    * */
  def validatePrimers(primers: Seq[Primer], multiPrimerPairs: Boolean): Unit = {
    import com.fulcrumgenomics.FgBioDef.javaIteratorAsScalaIterator

    if (multiPrimerPairs) {
      // Validate that there is at least one forward and one reverse primer for each pair
      primers.groupBy(_.pair_id).foreach {
        case (pairId, primersForPair) if primersForPair.length >= 2 =>
          if (!primersForPair.exists(p => p.positive_strand) || primersForPair.forall(p => p.positive_strand)) {
            val tpe = if (!primersForPair.exists(p => p.positive_strand)) "forward" else "reverse"
            throw new IllegalArgumentException(s"Found two $tpe primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
          }
        case (pairId, primersForPair) =>
          throw new IllegalArgumentException(s"Found ${primersForPair.length} primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
      }
    }
    else {
      // Validate that two primers exist for each pair, and that one is forward and the other is reverse
      primers.groupBy(_.pair_id).foreach {
        case (pairId, Seq(first, second)) =>
          if (first.positive_strand == second.positive_strand) {
            val tpe = if (first.positive_strand) "forward" else "reverse"
            throw new IllegalArgumentException(s"Found two $tpe primers for pair with id '$pairId': " + Seq(first, second).map(_.primer_id).mkString(", "))
          }
        case (pairId, primersForPair) =>
          throw new IllegalArgumentException(s"Found ${primersForPair.length} primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
      }
    }

    // Validate we do not have the same refName/start/end for primers
    val fwdDetector = {
      val detector = new OverlapDetector[Primer](0, 0)
      primers.filter(_.positive_strand).filter(_.mapped).foreach { p => detector.addLhs(p, p) }
      detector
    }
    val revDetector = {
      val detector = new OverlapDetector[Primer](0, 0)
      primers.filterNot(_.positive_strand).filter(_.mapped).foreach { p => detector.addLhs(p, p) }
      detector
    }
    primers.filter(_.mapped).foreach { primer =>
      val overlaps = (fwdDetector.getOverlaps(primer).iterator() ++ revDetector.getOverlaps(primer).iterator())
        .filter { overlap => primer != overlap && overlap.positive_strand == primer.positive_strand && overlap.ref_name == primer.ref_name }
        .filter { overlap => overlap.start == primer.start &&  overlap.end == primer.end }
      if (overlaps.nonEmpty) {
        throw new IllegalArgumentException(s"Found primers had the same location:\n$primer\n" + overlaps.map(_.toString).mkString("\n"))
      }
    }
  }
}

/** A locatable class for a single primer of a primer pair.
  *
  * The bases should be 5'->3' on the genomic forward strand, which facilitates matching against [[SamRecord]]s.
  *
  * @param pair_id the canonical primer pair identifier, unique across all primer pairs.
  * @param primer_id the canonical primer identifier, unique across all primers.
  * @param sequence the primer sequence, in sequencing order.
  * @param ref_name the reference name to which this primer targets.
  * @param start the reference start position at which this primer starts.
  * @param end the reference end position at which this primer starts.
  * @param positive_strand true if the primer maps to the forward genomic strand, false otherwise.
  * Primers without a mapping to the reference should have an empty `ref_name`.
  */
private[identifyprimers] case class Primer(pair_id: String,
                                           primer_id: String,
                                           sequence: String,
                                           ref_name: String,
                                           start: Int,
                                           end: Int,
                                           positive_strand: Boolean) extends Locatable with Alignable with Metric {
  override def getContig: String = ref_name
  override def getStart: Int = start
  override def getEnd: Int = end
  override def length: Int = end - start + 1
  def positiveStrand: Boolean = this.positive_strand
  def negativeStrand: Boolean = !this.positive_strand

  /** Returns true if the primer has a mapping to the reference, false otherwise. */
  def mapped: Boolean = this.ref_name.nonEmpty

  /** The bases as bytes in sequencing order */
  val bases: Array[Byte] = sequence.getBytes

  /** The bases in sequencing order. */
  def basesInSequencingOrder: Array[Byte] = this.bases

  /** The bases in genomic order. */
  val basesInGenomicOrder: Array[Byte] = if (positiveStrand) this.bases else SequenceUtil.reverseComplement(sequence).getBytes()

  // Validate the primer bases
  this.sequence.zipWithIndex.foreach { case (base, index) =>
    if (!SequenceUtil.isValidBase(base.toByte) && !SequenceUtil.isIUPAC(base.toByte)) {
      val prefix = sequence.substring(0, index)
      val base   = sequence.charAt(index)
      val suffix = sequence.substring(index+1)
      throw new IllegalArgumentException(s"Found an invalid base for primer $pair_id/$primer_id: $prefix[$base]$suffix")
    }
  }
}