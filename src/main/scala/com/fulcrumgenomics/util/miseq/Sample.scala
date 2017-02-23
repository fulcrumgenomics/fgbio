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
package com.fulcrumgenomics.util.miseq

import com.fulcrumgenomics.util.ReadStructure
import dagr.commons.CommonsDef._

/**
  * Represents information about a sample within an Illumina Experiment Manager sample sheet.
  *
  * Optionally contains information about dual indexes: i7 and i5.
  *
  * @author Nils Homer
  * @param sampleOrdinal the sample ordinal if this sample belongs to a sample set.
  * @param sampleId the unique sample identifier.
  * @param sampleName the sample name.
  * @param libraryId the library identifier.
  * @param project the project identifier.
  * @param description the sample description.
  * @param lane the lane number for the sample.
  * @param i7IndexBases the sample barcode bases in the i7 read.
  * @param i5IndexBases the sample barcode bases in the i5 read.
  * @param sampleBarcode the full sample barcode for optimized access.
  * @param extendedAttributes a map of non-standard or site-specific extended attributes. Keys should be upper-case.
  */
class Sample(val sampleOrdinal: Int,
             val sampleId: String,
             val sampleName: String,
             val libraryId: Option[String] = None,
             val project: Option[String] = None,
             val description: Option[String] = None,
             val lane: Option[Int] = None,
             val i7IndexBases: Option[String] = None,
             val i5IndexBases: Option[String] = None,
             var sampleBarcode: Option[SampleBarcode] = None,
             val extendedAttributes: Map[String, String] = Map.empty) {
  extendedAttributes.keysIterator.foreach { key => require(key == key.toUpperCase) }

  /** Returns the sample barcodes in order of sequencing. */
  def sampleBarcodeBases: Seq[Option[String]] = {
    Seq(i7IndexBases, i5IndexBases)
  }

  /** Gets an extended attribute for the given column name.  The column name can be any case. */
  def extendedAttribute(name: String): Option[String] = {
    extendedAttributes.get(name.toUpperCase) match {
      case None => None
      case Some(str) if str.isEmpty => None
      case Some(str) => Some(str.trim)
    }
  }

  /** Sets the sample barcode based on the read structure.
    *
    * @param sampleBarcodeReadStructures a read structure per read (ex. i7 and i5) that contains only sample barcode segments.
    */
  def setSampleBarcode(sampleBarcodeReadStructures: Seq[ReadStructure]): this.type = {
    import com.fulcrumgenomics.util.{SampleBarcode => ReadSegmentSampleBarcode}
    if (sampleBarcodeReadStructures.exists(rs => !rs.forall { case s: ReadSegmentSampleBarcode => true; case _ => false } )) {
      throw new IllegalArgumentException("Read structures contained segments that were not sample barcodes.")
    }
    val sampleBarcodeBases = this.sampleBarcodeBases.flatten
    if (sampleBarcodeReadStructures.length != sampleBarcodeBases.length) unreachable("Different # of read structures than sample barcode reads.")
    val sampleBarcodes = sampleBarcodeBases.zip(sampleBarcodeReadStructures).map { case (bases, rs) =>
      if (rs.map(_.length).sum != bases.length) {
        throw new IllegalArgumentException(s"Different number of sample barcode bases than in the read structure")
      }
      rs.structureRead(bases).map(_.bases).mkString
    }
    if (sampleBarcodes.nonEmpty) this.sampleBarcode = Some(new SampleBarcode(sampleBarcodes))
    this
  }
}

