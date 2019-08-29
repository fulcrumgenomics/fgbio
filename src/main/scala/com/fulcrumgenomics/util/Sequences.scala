/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools.util.SequenceUtil

import scala.annotation.switch


/**
  * Utility methods for working with DNA or RNA sequences
  */
object Sequences {
  /** Common contig/chrom names for non-autosomal sequences in mammals. */
  val CommonNonAutosomalContigNames: Seq[String] = Seq("M", "chrM", "MT", "X", "chrX", "Y", "chrY")

  /** Array for determining what bases are compatible. Indices are all DNA/RNA bases. Contents
    * are Int masks with a bit set for each possible base.  U is treated as identical to T.
    */
  private val IupacMasks: Array[Int] = {
    val masks = new Array[Int](Byte.MaxValue)
    val (a, c, g, t) = (1, 2, 4, 8)

    def up(c: Char): Int = c.toUpper.toInt
    def down(c: Char): Int = c.toLower.toInt

    Seq(up(_), down(_)).foreach { f =>
      masks(f('A')) = a
      masks(f('C')) = c
      masks(f('G')) = g
      masks(f('T')) = t
      masks(f('U')) = t
      masks(f('M')) = a | c
      masks(f('R')) = a | g
      masks(f('W')) = a | t
      masks(f('S')) = c | g
      masks(f('Y')) = c | t
      masks(f('K')) = g | t
      masks(f('V')) = a | c | g
      masks(f('H')) = a | c | t
      masks(f('D')) = a | g | t
      masks(f('B')) = c | g | t
      masks(f('N')) = a | c | g | t
    }

    masks
  }

  /** Returns true if two bases are compatible.  Compatibility is defined as having any overlap in the
    * set of acceptable bases.  E.g.:
    *   - `compatible('T', 't') == true`
    *   - `compatible('T', 'U') == true`
    *   - `compatible('A', 'R') == true`
    *   - `compatible('M', 'S') == true`
    *   - `compatible('A', 'C') == false`
    *   - `compatible('M', 'K') == false`
    *
    * @param base1 the first base to be compared
    * @param base2 the second base to be compared
    * @return true if the bases share at least one concrete base in common, false otherwise
    */
  final def compatible(base1: Byte, base2: Byte): Boolean = base1 == base2 || (IupacMasks(base1) & IupacMasks(base2)) > 0


  /** Counts the number of mismatches between two sequences of the same length. */
  def countMismatches(s1: String, s2: String): Int = {
    require(s1.length == s2.length, s"Cannot count mismatches in strings of differing lengths: $s1 $s2")

    var count = 0
    forloop (from=0, until=s1.length) { i =>
      val a = Character.toUpperCase(s1.charAt(i))
      val b = Character.toUpperCase(s2.charAt(i))
      if (a != b) count += 1
    }

    count
  }

  /** Class to store the zero-based offset and length of match for various properties of a sequence.  For example, see
    * [[longestHomopolymer]] and [[longestDinuc]]. */
  case class OffsetAndLength(offset: Int, length: Int)

  /**
    * Returns the offset (0-based) and length of the longest homopolymer in the sequence. In the case
    * that there are multiple homopolymers of the same length, the earliest one is returned.
    *
    * @param s a DNA or RNA sequence
    * @return the offset and length of the longest homopolymer
    */
  def longestHomopolymer(s: String) : OffsetAndLength = {
    var (bestStart, bestLength) = (0,0)
    forloop(0, s.length) { start =>
      val firstBase = s.charAt(start).toByte
      var i = start
      while (i < s.length && basesEqual(firstBase, s.charAt(i))) i += 1
      val length = i - start
      if (length > bestLength) {
        bestStart = start
        bestLength = length
      }
    }

    OffsetAndLength(bestStart, bestLength)
  }

  /**
    * Returns the offset (0-based) and length of the longest dinucleotide run in the sequence. In the case
    * that there are multiple dinucleotide runs of the same length, the earliest one is returned.
    *
    * @param s a DNA or RNA sequence
    * @return the offset and length of the longest dinucleotide sequence
    */
  def longestDinuc(s: String) : OffsetAndLength = {
    var (bestStart, bestLength) = (0,0)
    forloop(0, s.length-1) { start =>
      val (b1, b2) = (s.charAt(start).toByte, s.charAt(start+1).toByte)
      var i = start
      while (i < s.length-1 && basesEqual(b1, s.charAt(i)) && basesEqual(b2, s.charAt(i+1))) i += 2
      val length = i - start
      if (length > bestLength) {
        bestStart = start
        bestLength = length
      }
    }

    OffsetAndLength(bestStart, bestLength)
  }

  /** Reverse complements a string of bases. */
  def revcomp(s: String): String = {
    val bs = s.getBytes
    revcomp(bs)
    new String(bs)
  }

  /** Returns the sequence that is the complement of the provided sequence. */
  def complement(s: String): String = {
    val bs = s.getBytes
    complement(bs)
    new String(bs)
  }

  /** Reverse complements an array of bases in place. See [[complement()]] for how bases are complemented. */
  def revcomp(bs: Array[Byte]): Unit = {
    var (i, j) = (0, bs.length - 1)
    while (i < j) {
      val tmp = bs(i)
      bs(i) = complement(bs(j))
      bs(j) = complement(tmp)
      i += 1
      j -= 1
    }

    if (i == j) bs(i) = complement(bs(i))
  }

  /**
    * Complements the bases in the array in place. See [[complement()]] for how complementing is performed.
    *
    * @param bs an array of bases as bytes to be complemented in place
    */
  def complement(bs: Array[Byte]): Unit = forloop (from=0, until=bs.length) { i => bs(i) = complement(bs(i)) }

  /** Complements a single base. `U` and `u` are always complemented to `A` and `a` respectively.
    * IUPAC ambiguity codes are handled appropriately - they are complemented to the IUPAC code that represents
    * the set of bases that can pair with the input IUPAC code.
    *
    * Any byte that is not a valid upper- or lower-case DNA base, RNA base or IUPAC code will cause
    * an exception to be thrown.
    *
    * @param b the base to be complemented
    * @return the complement of that base
    */
  @inline final def complement(b: Byte): Byte = (b: @switch) match {
    // Discrete bases
    case 'A' => 'T'
    case 'C' => 'G'
    case 'G' => 'C'
    case 'T' => 'A'
    case 'U' => 'A'
    // IUPAC codes that represent two bases
    case 'M' => 'K'
    case 'K' => 'M'
    case 'R' => 'Y'
    case 'Y' => 'R'
    case 'W' => 'S'
    case 'S' => 'W'
    // IUPAC codes that represent three bases
    case 'B' => 'V'
    case 'V' => 'B'
    case 'H' => 'D'
    case 'D' => 'H'
    // IUPAC universal code
    case 'N' => 'N'

    // Discrete bases
    case 'a' => 't'
    case 'c' => 'g'
    case 'g' => 'c'
    case 't' => 'a'
    case 'u' => 'a'
    // IUPAC codes that represent two bases
    case 'm' => 'k'
    case 'k' => 'm'
    case 'r' => 'y'
    case 'y' => 'r'
    case 'w' => 's'
    case 's' => 'w'
    // IUPAC codes that represent three bases
    case 'b' => 'v'
    case 'v' => 'b'
    case 'h' => 'd'
    case 'd' => 'h'
    // IUPAC universal code
    case 'n' => 'n'

    case ch  => throw new IllegalArgumentException(s"Don't know how to complement '$ch'.")
  }

  /** Compares if two bases are equal ignoring case.  The bases may be IUPAC codes, but their relationships are not
    * considered. */
  @inline private def basesEqual(b1: Byte, b2: Char): Boolean = SequenceUtil.basesEqual(b1, b2.toByte)

  /** Returns the IUPAC code for the set of bases given. */
  def iupacCode(bases: Iterable[Byte]): Byte = {
    val mask: Int = bases.foldLeft(0) { case (code, next) => code | IupacMasks(next) }
    IupacMasks.indexOf(mask).toByte
  }
}
