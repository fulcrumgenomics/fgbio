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
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.util.{Locatable, OverlapDetector, SequenceUtil}

private[identifyprimers] object Primer {

  /** Writes the given primers to file.  Reverse strand primers are written with their sequence reverse complemented. */
  def write(path: FilePath, primers: TraversableOnce[Primer]): Unit = {
    val newPrimers = primers.map { primer =>
      if (primer.positiveStrand) primer
      else primer.copy(sequence = SequenceUtil.reverseComplement(primer.sequence))
    }
    Metric.write(path, newPrimers)
  }

  /** Reads primers from a given path and validates the primers (see [[validatePrimers()]]). */
  def read(path: FilePath, multiPrimerPairs: Boolean = false): Seq[Primer] = {
    val parser  = DelimitedDataParser(path, '\t')
    val primers = parser.map { row =>
      val forward  = strandToForward(row.apply[String]("strand"))
      val sequence = row.apply[String]("sequence").toUpperCase
      Primer(
        pair_id   = row.apply[String]("pair_id"),
        primer_id = row.apply[String]("primer_id"),
        sequence  = if (forward) sequence else SequenceUtil.reverseComplement(sequence),
        ref_name  = row.apply[String]("ref_name"),
        start     = row.apply[Int]("start"),
        end       = row.apply[Int]("end"),
        positiveStrand   = forward
      )
    }.toSeq
    validatePrimers(primers, multiPrimerPairs)
    primers
  }

  /** Converts strand to a True if forward, false if reverse. */
  private def strandToForward(strand: String): Boolean = strand.toLowerCase match {
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

    val fwdDetector = {
      val detector = new OverlapDetector[Primer](0, 0)
      primers.filter(_.positiveStrand).foreach { p => detector.addLhs(p, p) }
      detector
    }
    val revDetector = {
      val detector = new OverlapDetector[Primer](0, 0)
      primers.filterNot(_.positiveStrand).foreach { p => detector.addLhs(p, p) }
      detector
    }

    if (multiPrimerPairs) {
      // Validate that there is at least one forward and one reverse primer for each pair
      primers.groupBy(_.pair_id).foreach {
        case (pairId, primersForPair) if primersForPair.length >= 2 =>
          if (!primersForPair.exists(p => p.positiveStrand) || primersForPair.forall(p => p.positiveStrand)) {
            val tpe = if (!primersForPair.exists(p => p.positiveStrand)) "forward" else "reverse"
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
          if (first.positiveStrand == second.positiveStrand) {
            val tpe = if (first.positiveStrand) "forward" else "reverse"
            throw new IllegalArgumentException(s"Found two $tpe primers for pair with id '$pairId': " + Seq(first, second).map(_.primer_id).mkString(", "))
          }
        case (pairId, primersForPair) =>
          throw new IllegalArgumentException(s"Found ${primersForPair.length} primers for pair with id '$pairId': " + primersForPair.map(_.primer_id).mkString(", "))
      }
    }

    // Validate we do not have forward primers that have the same refName/start or reverse primers with the same refName/end
    primers.foreach { primer =>
      val overlaps = (fwdDetector.getOverlaps(primer).iterator() ++ revDetector.getOverlaps(primer).iterator())
        .filter { overlap => primer != overlap && overlap.positiveStrand == primer.positiveStrand && overlap.ref_name == primer.ref_name }
        .filter { overlap => if (primer.positiveStrand) overlap.start == primer.start else overlap.end == primer.end }
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
  * @param sequence the primer sequence, in 5'->3' on the genomic forward strand.
  * @param ref_name the reference name to which this primer targets.
  * @param start the reference start position at which this primer starts.
  * @param end the reference end position at which this primer starts.
  * @param positiveStrand true if the primer maps to the forward genomic strand, false otherwise.
  *
  * Primers without a mapping to the reference should have an empty `ref_name`.
  */
private[identifyprimers] case class Primer(pair_id: String,
                                           primer_id: String,
                                           sequence: String,
                                           ref_name: String,
                                           start: Int,
                                           end: Int,
                                           positiveStrand: Boolean,
                                           private val reverseComplementBases: Boolean = false) extends Locatable with Alignable with Metric {
  override def getContig: String = ref_name
  override def getStart: Int = start
  override def getEnd: Int = end
  override def length: Int = end - start + 1
  def negativeStrand: Boolean = !this.positiveStrand

  /** Returns true if the primer has a mapping to the reference, false otherwise. */
  def mapped: Boolean = this.ref_name.nonEmpty

  val bases: Array[Byte] = if (reverseComplementBases) SequenceUtil.reverseComplement(sequence).getBytes else sequence.getBytes

  this.sequence.zipWithIndex.foreach { case (base, index) =>
    if (!SequenceUtil.isValidBase(base.toByte) && !SequenceUtil.isIUPAC(base.toByte)) {
      val prefix = sequence.substring(0, index)
      val base   = sequence.charAt(index)
      val suffix = sequence.substring(index+1)
      throw new IllegalArgumentException(s"Found an invalid base for primer $pair_id/$primer_id: $prefix[$base]$suffix")
    }
  }
}