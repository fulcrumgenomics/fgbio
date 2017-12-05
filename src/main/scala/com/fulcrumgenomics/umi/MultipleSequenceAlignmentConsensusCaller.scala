/*
 * The MIT License
 *
 * Copyright (c) $year Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.UmiConsensusCaller.SourceRead
import com.fulcrumgenomics.util.NumericTypes.PhredScore

import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag
import scala.util.{Failure, Random, Success, Try}
import sys.process._

class MultipleSequenceAlignmentConsensusCaller(val msaCommand: String = DuplexConsensusCaller.MsaCommand,
                                               val fallbackToVanillaCaller: Boolean = false,
                                               val options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions()) extends ConsensusCallerTrait {

  private val NoCall: Byte                   = 'N'.toByte
  private val Gap: Byte                      = '-'.toByte
  private val GapQual: PhredScore            = 30.toByte
  private val NotEnoughReadsQual: PhredScore = 0.toByte // Score output when masking to N due to insufficient input reads
  private val TooLowQualityQual: PhredScore  = 2.toByte  // Score output when masking to N due to too low consensus quality
  private val vanillaCaller = new VanillaUmiConsensusCaller(readNamePrefix="x", options=options)

  private val random = new Random(42)
  private val caller = new ConsensusCaller(errorRatePreLabeling=options.errorRatePreUmi, errorRatePostLabeling=options.errorRatePostUmi)

  private var numMsaFailures: Long = 0

  private def updateReads(reads: Seq[SourceRead], result: Seq[String]): Seq[SourceRead] = {
    // create one string builder for each read
    val basesBuilders = reads.map { _ => new ArrayBuffer[Byte] }.toArray

    // aggregate bases for each read
    result
      .drop(3)
      .filterNot { line => line.startsWith(" ") || line.startsWith("*") || line.isEmpty }
      .zipWithIndex.foreach { case (line, i) =>
        val Seq(id, bases) = line.split("[ ]+").toSeq
        basesBuilders(id.toInt) ++= bases.map(_.toByte)
      }

    val expectedLength = basesBuilders.head.length
    require(basesBuilders.forall(_.length == expectedLength))

    basesBuilders.zipWithIndex.map { case (bases, i) =>
      // pad the qualities
      val read      = reads(i)
      var qualIndex = 0
      val quals = bases.map { b =>
        if (b != Gap) read.quals(qualIndex)
        else {
          qualIndex += 1
          Gap
        }
      }.toArray
      read.copy(bases=bases.toArray, quals=quals)
    }.toSeq
  }

  /**
    * Calculates the length of the consensus read that should be produced. The length is calculated
    * as the maximum length at which minReads reads still have bases.
    *
    * @param reads the set of reads being fed into the consensus
    * @param minReads the minimum number of reads required
    * @return the length of consensus read that should be created
    */
  private def consensusReadLength(reads: Seq[SourceRead], minReads: Int): Int = {
    require(reads.size >= minReads, "Too few reads to create a consensus.")
    reads.map(_.length).sortBy(len => -len).drop(minReads-1).head
  }

  private def align(reads: Seq[SourceRead], fastaPath: FilePath): Try[Seq[SourceRead]] = {
    // First limit to max reads if necessary
    //val capped = if (reads.size <= 32) reads else this.random.shuffle(reads).take(32)
    val capped = if (reads.size <= this.options.maxReads) reads else this.random.shuffle(reads).take(this.options.maxReads)

    // Create a FASTA
    val lines = capped.zipWithIndex.flatMap { case (read, id) => Seq(s">$id", read.baseString) }
    Io.writeLines(fastaPath, lines)

    // Run the MSA tool
    // NB: can run pretty much anything that produces valid CLUSTAL output.
    Try { s"$msaCommand $fastaPath" !! ProcessLogger(_ => ()) split '\n' }.map { result => updateReads(capped, result) }
  }

  private def call(aligned: Seq[SourceRead]): VanillaConsensusRead = {
    // get the most likely consensus bases and qualities
    val consensusLength = consensusReadLength(aligned, this.options.minReads)
    val consensusBases  = new Array[Base](consensusLength)
    val consensusQuals  = new Array[PhredScore](consensusLength)
    val consensusDepths = new Array[Short](consensusLength)
    val consensusErrors = new Array[Short](consensusLength)

    var positionInRead = 0
    val builder = this.caller.builder()
    while (positionInRead < consensusLength) {
      val reads = aligned.filter(_.length > positionInRead)

      // Find the base that has the fewest # of observations, and ignore it, and use it for gap in the model
      val baseCounter = new SimpleCounter[Byte]()
      Seq('A', 'C', 'T', 'G').foreach { base => baseCounter.count(base.toByte, 0) }
      reads.map { _.bases(positionInRead) }.filterNot(_ == NoCall).foreach { base => baseCounter.count(base) }
      val zeroObservationsBase = baseCounter.minBy(_._2)._1

      // Add the evidence from all reads that are long enough to cover this base
      reads.foreach { read =>
        val base = read.bases(positionInRead)
        val qual = read.quals(positionInRead)
        if (base != NoCall && base != zeroObservationsBase) {
          if (base == Gap) builder.add(base=zeroObservationsBase, qual=GapQual)
          else builder.add(base=base, qual=qual)
        }
      }

      // Call the consensus and do any additional filtering
      val (rawBase, rawQual) = builder.call()
      val (base, qual) = {
        if (builder.contributions < this.options.minReads)       (NoCall, NotEnoughReadsQual)
        else if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual)
        else (rawBase, rawQual)
      }

      consensusBases(positionInRead) = if (base == zeroObservationsBase) Gap else base
      consensusQuals(positionInRead) = qual

      // Generate the values for depth and count of errors
      val depth  = builder.contributions
      val errors = if (rawBase == NoCall) depth else depth - builder.observations(rawBase)
      consensusDepths(positionInRead) = if (depth  > Short.MaxValue) Short.MaxValue else depth.toShort
      consensusErrors(positionInRead) = if (errors > Short.MaxValue) Short.MaxValue else errors.toShort

      // Get ready for the next pass
      builder.reset()
      positionInRead += 1
    }

    val gapMask = consensusBases.map { base => base == Gap }
    def filterGap[T](things: Array[T])(implicit ev1: ClassTag[T]): Array[T] = things.toList.zip(gapMask).filterNot(_._2).map(_._1).toArray

    VanillaConsensusRead(
      id     = aligned.head.id,
      bases  = filterGap(consensusBases),
      quals  = filterGap(consensusQuals),
      depths = filterGap(consensusDepths),
      errors = filterGap(consensusErrors)
    )
  }

  /** Creates a consensus read from the given read and qualities sequences.
    * If no consensus read was created, None is returned.
    *
    * The same number of base sequences and quality sequences should be given.
    * */
  def consensusCall(reads: Seq[SourceRead]): Option[VanillaConsensusRead] = {
    // check to see if we have enough reads.
    if (reads.size < this.options.minReads) {
      None
    }
    else if (reads.size <= 2) {
      vanillaCaller.consensusCall(reads)
    }
    else {
      val fastaPath  = Io.makeTempFile(s"batch.${this.random.nextInt}.", ".fasta")
      align(reads, fastaPath) match {
        case Failure(thr)       =>
          fastaPath.toFile.delete()
          numMsaFailures += 1
          if (this.fallbackToVanillaCaller) {
            vanillaCaller.consensusCall(reads)
          } else {
            throw new IllegalStateException(s"MSA command failed: '${this.msaCommand}'", thr)
          }
        case Success(aligned) =>
          fastaPath.toFile.delete()
          Some(call(aligned))

      }
    }
  }
}
