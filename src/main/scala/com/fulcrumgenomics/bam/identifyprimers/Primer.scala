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

import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.{LazyLogging, StringUtil}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.util.{Locatable, SequenceUtil}

private[identifyprimers] object Primer extends LazyLogging {
  
  private val Headers = Seq(
    "pair_id",
     "primer_id",
     "sequence",
     "ref_name",
     "start",
     "end",
     "strand"
  )

  /** Writes the given primers to file. */
  def write(path: FilePath, primers: TraversableOnce[Primer]): Unit = Metric.write(path, primers)

  /** Reads primers from a given path and validates the primers (see [[validatePrimers()]]).
    *
    * NB: the reverse strand primers stored at the path are assumed to have sequence in the primer order (i.e. not the
    * genomic top strand).
    * */
  def read(path: FilePath, multiPrimerPairs: Boolean = false): Seq[Primer] = {
    val lines = Io.toSource(path).getLines()
    val headers = if (lines.hasNext) lines.next().split('\t').toIndexedSeq else IndexedSeq.empty

    // A temporary array used for parsing each line into fields
    val tmp = new Array[String](headers.size * 2)

    def toInt(str: String): Int = if (str == "") 0 else str.toInt

    val primers = lines.zipWithIndex.map { case (line, lineIndex) =>
      // Get the fields
      val fields = {
        val count  = StringUtil.split(line=line, delimiter='\t', arr=tmp)
        val f = new Array[String](count)
        System.arraycopy(tmp, 0, f, 0, count)
        f
      }

      // Parse the fields
      val row = headers.zip(fields).toMap
      fields.length match {
        case 3 => // assume no mapping information
          Headers.take(3).find(h => !row.contains(h)).foreach { h =>
            throw new IllegalStateException(s"Missing value for $h on row ${lineIndex+1}")
          }
          Primer(
            pair_id         = row("pair_id"),
            primer_id       = row("primer_id"),
            sequence        = row("sequence").toUpperCase,
            ref_name        = "",
            start           = 0,
            end             = 0,
            positive_strand = true
          )
        case headers.length => // all columns must exist
          Headers.find(h => !row.contains(h)).foreach { h =>
            throw new IllegalStateException(s"Missing value for $h on row ${lineIndex+1}")
          }
          val positive_strand = row("strand") == "" || isPositiveStrand(row("strand"))
          Primer(
            pair_id         = row("pair_id"),
            primer_id       = row("primer_id"),
            sequence        = row("sequence").toUpperCase,
            ref_name        = row("ref_name"),
            start           = toInt(row("start")),
            end             = toInt(row("end")),
            positive_strand = positive_strand
          )
        case _ =>
          logger.error(s"Line has incorrect number of fields. Header length is ${headers.length}, row length is ${fields.length}")
            logger.error(s"Headers: ${headers.mkString("\t")}")
            logger.error(s"Line: $line")
            throw new IllegalStateException(s"Incorrect number of values in line: ${lineIndex+1}.")
      }
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
    * Also validates that no primers have the same ref_name/start/stop/strand
    * */
  def validatePrimers(primers: Seq[Primer], multiPrimerPairs: Boolean): Unit = {
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

    // Validate we do not have the same refName/start/end/strand for primers
    primers.sortBy { primer => (primer.positiveStrand, primer.ref_name, primer.start, primer.end) }
      .sliding(2)
      .filter { case Seq(primer, other) =>
        primer.ref_name == other.ref_name && primer.start == other.start &&
          primer.end == other.end && primer.positiveStrand == other.positiveStrand
      }
      .foreach { case Seq(primer, other) =>
        throw new IllegalArgumentException(s"Found multiple primers with the same coordinates: $primer and $other")
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
                                           positive_strand: Boolean) extends Locatable with Metric {
  require(!pair_id.contains(','), s"pair_id contained a comma: '$pair_id'")
  require(!primer_id.contains(','), s"primer_id contained a comma: '$pair_id'")

  // Validate the primer bases
  this.sequence.zipWithIndex.foreach { case (base, index) =>
    if (!SequenceUtil.isValidBase(base.toByte) && !SequenceUtil.isIUPAC(base.toByte)) {
      val prefix = sequence.substring(0, index)
      val base   = sequence.charAt(index)
      val suffix = sequence.substring(index+1)
      throw new IllegalArgumentException(s"Found an invalid base for primer $pair_id/$primer_id: $prefix[$base]$suffix")
    }
  }

  override def getContig: String = ref_name
  override def getStart: Int = start
  override def getEnd: Int = end
  def length: Int = if (ref_name.nonEmpty) end - start + 1 else sequence.length
  def positiveStrand: Boolean = this.positive_strand
  def negativeStrand: Boolean = !this.positive_strand

  /** Returns true if the primer has a mapping to the reference, false otherwise. */
  def mapped: Boolean = this.ref_name.nonEmpty

  /** The bases in sequencing order. */
  val basesInSequencingOrder: Array[Byte] = sequence.getBytes

  /** The bases as a [[String]] in sequencing order. */
  def baseStringInSequencingOrder: String = sequence

  /** The bases in genomic order. */
  val basesInGenomicOrder: Array[Byte] = if (positiveStrand) basesInSequencingOrder else SequenceUtil.reverseComplement(sequence).getBytes()
}