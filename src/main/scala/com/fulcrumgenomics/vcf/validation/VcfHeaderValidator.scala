/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf.validation

import com.fulcrumgenomics.vcf.api._
import com.fulcrumgenomics.vcf.validation.ValidationResult._

import scala.collection.mutable.ListBuffer

sealed trait VcfHeaderValidator extends VcfValidator {
  def validate(header: VcfHeader): Seq[ValidationResult]
}

object VcfHeaderValidator {
  val Validators: Seq[VcfHeaderValidator] = Seq(
    new VcfHeaderHasContigsValidator,
    new VcfHeaderUniqueContigsValidator,
    new VcfHeaderUniqueInfosValidator,
    new VcfHeaderUniqueFormatsValidator
  )

  class VcfHeaderHasContigsValidator extends VcfHeaderValidator {
    def validate(header: VcfHeader): Seq[ValidationResult] = {
      if (header.contigs.nonEmpty) Seq.empty
      else Seq(warning(message="No contig lines in the header."))
    }
  }

  sealed trait UniqueValuesValidator[T<:VcfHeaderEntry] extends VcfHeaderValidator {
    def toEntries(header: VcfHeader): Seq[T]
    def toKey(entry: T): String
    def name: String
    def validate(header: VcfHeader): Seq[ValidationResult] = {
      toEntries(header=header)
        .groupBy(toKey)
        .filter(_._2.length > 2)
        .map { case (key, values) =>
          error(s"Found ${this.name} name ${values.length} times: `$key`")
        }.toIndexedSeq
    }
  }

  class VcfHeaderUniqueContigsValidator extends UniqueValuesValidator[VcfContigHeader] {
    def toEntries(header: VcfHeader): Seq[VcfContigHeader] = header.contigs
    def toKey(entry: VcfContigHeader): String = entry.name
    def name: String = "contig"
  }

  class VcfHeaderUniqueInfosValidator extends UniqueValuesValidator[VcfInfoHeader] {
    def toEntries(header: VcfHeader): Seq[VcfInfoHeader] = header.infos
    def toKey(entry: VcfInfoHeader): String = entry.id
    def name: String = "INFO.ID"
  }

  class VcfHeaderUniqueFormatsValidator extends UniqueValuesValidator[VcfFormatHeader] {
    def toEntries(header: VcfHeader): Seq[VcfFormatHeader] = header.formats
    def toKey(entry: VcfFormatHeader): String = entry.id
    def name: String = "FORMAT.ID"
  }
}

sealed trait VcfHeaderEntryValidator extends VcfValidator {
  def validate(entry: VcfHeaderEntry): Seq[ValidationResult]
}

object VcfHeaderEntryValidator {

  val ReservedVcfHeaderEntryValidators: Seq[VcfHeaderEntryValidator] = {
    VcfHeader.ReservedVcfInfoHeaders.map(VcfInfoHeaderValidator.apply) ++ VcfHeader.ReservedVcfFormatHeaders.map(VcfFormatHeaderValidator.apply)
  }

  trait VcfInfoOrFormatHeaderValidator extends VcfHeaderEntryValidator {
    def base: VcfHeaderInfoOrFormatEntry
    protected def errorPrefix: String

    def validate(entry: VcfHeaderEntry): Seq[ValidationResult] = entry match {
      case _entry: VcfHeaderInfoOrFormatEntry if _entry.id == base.id =>
        val builder = new ListBuffer[ValidationResult]()
        // Count
        if (base.count != _entry.count) {
          builder.append(error(
            f"$errorPrefix: expected count `${base.count}` found `${_entry.count}`"
          ))
        }
        // Kind
        if (base.kind != _entry.kind) {
          builder.append(error(
            f"$errorPrefix: expected type `${base.kind}` found `${_entry.kind}`"
          ))
        }
        if (builder.isEmpty) IndexedSeq.empty else builder.toIndexedSeq
      case _ => Seq.empty
    }
  }

  case class VcfInfoHeaderValidator(base: VcfInfoHeader) extends VcfInfoOrFormatHeaderValidator {
    protected val errorPrefix: String = f"Header INFO.${base.id}"
  }

  case class VcfFormatHeaderValidator(base: VcfFormatHeader) extends VcfInfoOrFormatHeaderValidator {
    protected val errorPrefix: String = f"Header FORMAT.${base.id}"
  }
}

