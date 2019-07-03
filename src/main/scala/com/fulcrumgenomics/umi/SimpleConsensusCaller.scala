/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef.forloop
import com.fulcrumgenomics.util.NumericTypes.{LogProbability, PhredScore}

import scala.collection.immutable.BitSet

/**
  * A class that can be mixed in to call sequences represented as strings of the same length.
  *
  * If a non-ACGTN bases is encountered at any position in the strings, that position is not consensus called, and the
  * character from the first sequenced is used.  All strings should have the same non-ACGTN character at this position.
  *
  * @param errorRatePreLabeling The error probability prior to adding the UMI.  By default, set to almost surely zero.
  * @param errorRatePostLabeling The error probability post to adding the UMI.  By default, set to almost surely zero
  * @param qError The adjusted error probability.  By default set to Q20.
  */
private[umi] class SimpleConsensusCaller(val errorRatePreLabeling: Byte = 90.toByte,
                                         val errorRatePostLabeling: Byte = 90.toByte,
                                         qError: PhredScore = 20.toByte
                                        ) {
  val pError: LogProbability = LogProbability.fromPhredScore(this.qError)
  val pTruth: LogProbability = LogProbability.not(this.pError)

  /** A consensus caller used to generate consensus UMI sequences */
  private val consensusBuilder = new ConsensusCaller(
    errorRatePreLabeling  = this.errorRatePreLabeling,
    errorRatePostLabeling = this.errorRatePostLabeling
  ).builder()

  private val DnaBases = Set('A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n')
  private val DnaBasesBitSet = BitSet(DnaBases.map(_.toInt).toSeq:_*)

  /** Calls a simple consensus sequences from a set of sequences all the same length. */
  def callConsensus(sequences: Seq[String]): String = {
    require(sequences.nonEmpty, "Can't call consensus on an empty set of sequences!")
    require(sequences.forall(_.length == sequences.head.length), "Sequences must all have the same length")
    val buffer = new StringBuilder
    val firstRead  = sequences.head
    val readLength = firstRead.length
    val sequencesLength = sequences.length

    forloop (from=0, until=readLength) { i =>
      this.consensusBuilder.reset()
      var nonDna = 0
      sequences.foreach { sequence =>
        val char = sequence.charAt(i)
        if (!this.DnaBasesBitSet.contains(char.toInt)) {
          nonDna += 1
          // verify that all non-DNA bases are the same character
          require(firstRead.charAt(i) == char,
            s"Sequences must have character '${firstRead.charAt(i)}' at position $i, found '$char'")
        }
        else this.consensusBuilder.add(char.toByte, pError=this.pError, pTruth=this.pTruth)
      }

      if (nonDna == 0) buffer.append(this.consensusBuilder.call()._1.toChar)
      else if (nonDna == sequencesLength) buffer.append(firstRead.charAt(i)) // NB: we have previously verified they are all the same character
      else throw new IllegalStateException(s"Sequences contained a mix of DNA and non-DNA characters at offset $i: $sequences")
    }

    buffer.toString()
  }
}
