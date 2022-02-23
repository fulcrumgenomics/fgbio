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

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.{SAMFileHeader, SAMProgramRecord, SAMReadGroupRecord, SAMTextHeaderCodec}

/**
  * Trait that can be mixed into any class that provides access to a SAMFileHeader in
  * order to provide some useful methods for more direct access.
  */
private[api] trait HeaderHelper {
  /** The associated [[SAMFileHeader]]. */
  def header: SAMFileHeader

  /** The sequence dictionary. */
  lazy val dict: SequenceDictionary = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    header.getSequenceDictionary.fromSam
  }

  /** The list of read groups from the associated header. */
  def readGroups: Seq[SAMReadGroupRecord] = header.getReadGroups.toIndexedSeq

  /** The list of program groups from the associated header. */
  def programGroups: Seq[SAMProgramRecord] = header.getProgramRecords.toIndexedSeq

  /** Gets the comments from the header, removing the stupid leading @CO that HTSJDK returns. */
  def comments: Seq[String] =
    header.getComments.iterator().map(_.replace(SAMTextHeaderCodec.COMMENT_PREFIX, "")).toIndexedSeq
}
