/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef._
import htsjdk.samtools.util.SequenceUtil

/**
  * Represents an allele in a VCF or [[Variant]].  The only requirement of the trait is that implementing
  * classes override [[toString]] to provide the correct string representation of the allele.
  */
sealed trait Allele {
  /** Returns the displayed value of the allele as a String. */
  def value: String

  /** Ensures that the toString is always the same as the value. */
  override final def toString: String = value
}

object Allele {
  // Constant alleles for all the expected single-base alleles
  private val A = SimpleAllele("A")
  private val C = SimpleAllele("C")
  private val G = SimpleAllele("G")
  private val T = SimpleAllele("T")
  private val N = SimpleAllele("N")
  private val a = SimpleAllele("a")
  private val c = SimpleAllele("c")
  private val g = SimpleAllele("g")
  private val t = SimpleAllele("t")
  private val n = SimpleAllele("n")

  /**
    * Returns an allele for the string in question.  The string may be a string of bases containing
    * `[ACGTNacgtn]`, the no-call string (`.`), a spanned allele (`*`) or a symbolic or break-end allele.
    *
    * The allele that is returned may be newly created or cached.

    * @param value the String value of the allele
    * @return an instance of Allele that represents the given String
    */
  def apply(value: String): Allele = if (value.length == 1) singleCharAllele(value) else value match {
    case v if v.forall(ch => SequenceUtil.isUpperACGTN(ch.toUpper.toByte)) => SimpleAllele(v)
    case v if v.startsWith("<") && v.endsWith(">") => SymbolicAllele(v)
    case v if v.exists(c => c == '[' || c == ']')  => BreakendAllele(v)
    case v => throw new IllegalArgumentException(s"Oops, should probably handle allele: $v")
  }

  /** Optimized match-based lookup for expected single-base allele constants. */
  private def singleCharAllele(s: String): Allele = s.charAt(0) match {
    case 'A' => A
    case 'C' => C
    case 'G' => G
    case 'T' => T
    case 'N' => N
    case 'a' => a
    case 'c' => c
    case 'g' => g
    case 't' => t
    case 'n' => n
    case '.' => NoCallAllele
    case '*' => SpannedAllele
    case ch  => throw new IllegalArgumentException(s"Can't create an allele for '$ch'.")
  }

  /** Singleton representing an allele that is uncalled. */
  case object NoCallAllele extends Allele {
    override val value: String = "."
  }

  /** Singleton representing an allele that doesn't exist because it is spanned by an upstream allele. */
  case object SpannedAllele extends Allele  {
    override val value: String = "*"
  }

  /** Class for alleles composed of a simple string of bases. */
  case class SimpleAllele private(bases: String) extends Allele {
    def length: Int = bases.length
    override def value: String = bases

    /** Provides an implementation of equality based on case-insensitive matching of the bases. */
    override def equals(other: Any): Boolean = other.isInstanceOf[SimpleAllele] && {
      val that = other.asInstanceOf[SimpleAllele]
      (this eq that) || this.bases.equalsIgnoreCase(that.bases)
    }

    /** Generate the hash based on the upper-case version of the bases. */
    override def hashCode(): Int = this.bases.toUpperCase.hashCode
  }

  /** Class for symbolic alleles. */
  case class SymbolicAllele private (id: String) extends Allele {
    override def value: String = id
  }

  /** Class for break-end alleles. */
  case class BreakendAllele private (override val value: String) extends Allele // TODO: what does this need to have
}
