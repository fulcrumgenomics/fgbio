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
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.mutable

object VariantMask {
  /** Generate a variant mask from the provided VCF. */
  def apply(path: PathToVcf): VariantMask = apply(new VCFFileReader(path.toFile))

  /** Generate a variant mask from the provided VCF reader. */
  def apply(reader: VCFFileReader): VariantMask = {
    val dict = reader.getFileHeader.getSequenceDictionary
    require(dict != null && dict.getSequences.size() > 0, "Generating a VariantMask requires VCFs contain contig lines.")
    apply(reader.iterator(), dict)
  }

  /** Generates a VariantMask from the variants in the provided iterator. */
  def apply(variants: Iterator[VariantContext], dict: SAMSequenceDictionary) = new VariantMask(variants, dict)
}

/**
  * Simple mask that loads variants one reference sequence at a time and creates a compact
  * representation allowing for rapid querying of whether or not positions are overlapped
  * by one or more variants.
  *
  * @param variants a coordinate sorted iterator of variants
  * @param dict the sequence dictionary for the reference
  */
class VariantMask(variants: Iterator[VariantContext], val dict: SAMSequenceDictionary) {
  private val iterator                    = variants.bufferBetter
  private var currentMask: mutable.BitSet = new mutable.BitSet(0)
  private var currentIndex: Int           = -1
  private var currentContig: String       = "n/a"
  private var currentRefLength: Int       = -1
  advanceTo(0)

  /**
    * Causes the variant mask to advance to the chosen reference sequence and load up the mask.
    */
  private def advanceTo(refIndex: Int): Unit = {
    require(refIndex > this.currentIndex, "Requested to navigate to an earlier contig. Is your VCF sorted?")
    require(refIndex >= 0 && refIndex < dict.size(), s"Requested navigation to contig #${refIndex} which is invalid.")

    val ref  = dict.getSequence(refIndex)
    val bits = new mutable.BitSet(ref.getSequenceLength)
    iterator.dropWhile(v => dict.getSequenceIndex(v.getContig) < refIndex)
      .takeWhile(v => dict.getSequenceIndex(v.getContig) == refIndex)
      .filterNot(v => v.isFiltered)
      .foreach { v =>
        forloop (from=v.getStart, until=v.getEnd+1) { i => bits(i-1) = true }
      }

    this.currentMask   = bits
    this.currentIndex  = refIndex
    this.currentContig = dict.getSequence(refIndex).getSequenceName
    this.currentRefLength = dict.getSequence(refIndex).getSequenceLength
  }

  /** Returns true if the position is overlapped by one or more variants, false otherwise. */
  def isVariant(refIndex: Int, position: Int): Boolean = {
    if (refIndex != this.currentIndex) advanceTo(refIndex)
    require(position > 0 && position <= this.currentRefLength, "Position outside of 1-contig length.")
    this.currentMask(position-1)
  }

  /** Returns true if the position is overlapped by one or more variants, false otherwise. */
  def isVariant(refName: String, position: Int): Boolean = isVariant(dict.getSequenceIndex(refName), position)
}
