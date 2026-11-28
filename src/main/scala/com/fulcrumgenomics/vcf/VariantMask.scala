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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}

import scala.collection.mutable

object VariantMask {
  /** Generate a variant mask from the provided VCF.  The mask reads the VCF lazily, one contig at a time, as it is
    * queried, so the VCF is left open; use the overload that takes a `VcfSource` to close it once done with the mask.
    */
  def apply(path: PathToVcf): VariantMask = {
    val source = VcfSource(path, allowMissingFieldsInHeader=true)
    try apply(source) catch { case ex: Throwable => source.safelyClose(); throw ex }
  }

  /** Generate a variant mask from the provided VCF source.  The mask reads the source lazily, one contig at a time, as
    * it is queried, so the source must not be closed until the mask is no longer queried. */
  def apply(source: VcfSource): VariantMask = {
    val dict = source.header.dict
    require(dict.length > 0, "Generating a VariantMask requires VCFs contain contig lines.")
    apply(source.iterator, dict)
  }

  /** Generates a VariantMask from the variants in the provided iterator. */
  def apply(variants: Iterator[Variant], dict: SequenceDictionary) = new VariantMask(variants, dict)
}

/**
  * Simple mask that loads variants one reference sequence at a time and creates a compact
  * representation allowing for rapid querying of whether or not positions are overlapped
  * by one or more variants.
  *
  * @param variants a coordinate sorted iterator of variants
  * @param dict the sequence dictionary for the reference
  */
class VariantMask(variants: Iterator[Variant], val dict: SequenceDictionary) {
  private val iterator                    = variants.bufferBetter
  private var currentMask: mutable.BitSet = new mutable.BitSet(0)
  private var currentIndex: Int           = -1
  private var currentRefLength: Int       = -1
  advanceTo(0)

  /**
    * Causes the variant mask to advance to the chosen reference sequence and load up the mask.
    */
  private def advanceTo(refIndex: Int): Unit = {
    require(refIndex > this.currentIndex, "Requested to navigate to an earlier contig. Is your VCF sorted?")
    require(refIndex >= 0 && refIndex < dict.length, s"Requested navigation to contig #$refIndex which is invalid.")

    val ref  = dict(refIndex)
    val bits = new mutable.BitSet(ref.length)
    iterator.dropWhile(v => dict(v.chrom).index < refIndex)
      .takeWhile(v => dict(v.chrom).index == refIndex)
      .filterNot(_.isFiltered)
      .foreach { v =>
        forloop (from=v.pos, until=v.end+1) { i => bits(i-1) = true }
      }

    this.currentMask   = bits
    this.currentIndex  = refIndex
    this.currentRefLength = dict(refIndex).length
  }

  /** Returns true if the position is overlapped by one or more variants, false otherwise. */
  def isVariant(refIndex: Int, position: Int): Boolean = {
    if (refIndex != this.currentIndex) advanceTo(refIndex)
    require(position > 0 && position <= this.currentRefLength, "Position outside of 1-contig length.")
    this.currentMask(position-1)
  }

  /** Returns true if the position is overlapped by one or more variants, false otherwise. */
  def isVariant(refName: String, position: Int): Boolean = isVariant(dict(refName).index, position)
}
