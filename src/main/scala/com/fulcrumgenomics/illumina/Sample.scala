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

package com.fulcrumgenomics.illumina

object Sample {
  val SampleBarcodeDelimiter: String = "-"
}

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
  * @param extendedAttributes a map of non-standard or site-specific extended attributes. Keys should be upper-case.
  */
case class Sample(sampleOrdinal: Int,
                  sampleId: String,
                  sampleName: String,
                  libraryId: String,
                  project: Option[String] = None,
                  description: Option[String] = None,
                  lane: Option[Int] = None,
                  i7IndexBases: Option[String] = None,
                  i5IndexBases: Option[String] = None,
                  extendedAttributes: Map[String, String] = Map.empty) {
  extendedAttributes.keysIterator.foreach {
    key => require(key == key.toUpperCase, s"Extended attribute key is not all upper case: $key")
  }

  import Sample.SampleBarcodeDelimiter

  def sampleBarcodes: Seq[Option[String]] = Seq(i7IndexBases, i5IndexBases)

  lazy val sampleBarcodeBytes: Array[Byte]  = sampleBarcodes.flatten.filter(_.nonEmpty).flatMap(_.getBytes).toArray
  lazy val sampleBarcodeString: String      = sampleBarcodes.flatten.filter(_.nonEmpty).mkString(SampleBarcodeDelimiter)

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
}

