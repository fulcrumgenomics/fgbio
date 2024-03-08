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

import com.fulcrumgenomics.vcf.api.{Genotype, Variant, VcfFormatHeader}
import com.fulcrumgenomics.vcf.api.VcfHeader.ReservedVcfFormatHeaders
import com.fulcrumgenomics.vcf.validation.ValidationResult.error


sealed trait GenotypeValidator extends VariantValidator {
  final def validate(variant: Variant): Seq[ValidationResult] = {
    variant.genotypes.values.flatMap { genotype => validate(variant=variant, genotype=genotype) }.toIndexedSeq
  }
  def validate(variant: Variant, genotype: Genotype): Seq[ValidationResult]
}

object GenotypeValidator {
  val VariantFormatValidators: Seq[VariantFormatValidator] = ReservedVcfFormatHeaders.map(VariantFormatValidator.apply)

  case class VariantFormatValidator(format: VcfFormatHeader) extends GenotypeValidator {
    def validate(variant: Variant, genotype: Genotype): Seq[ValidationResult] = {
      genotype.get[Any](format.id).toIndexedSeq.flatMap { value =>
        this.validate(kind=format.kind, count=format.count, variant=variant, genotype=Some(genotype), source=s"FORMAT.${format.id}", value=value)
      }
    }
  }

  case object PhaseSetGenotypeValidator extends GenotypeValidator {
    def validate(variant: Variant, genotype: Genotype): Seq[ValidationResult] = {
      (genotype.phased, genotype.get[Any]("PS")) match {
        case (false, None)            => Seq.empty // OK
        case (true, Some(value: Int)) => Seq.empty // OK
        case (false, Some(value))     => Seq(error(
          s"FORMAT.PS genotype was unphased but had a phase set value `$value`", variant=variant, genotype=genotype
        ))
        case (true, None)             => Seq(error(
          s"FORMAT.PS genotype was phased but had no phase set value", variant=variant, genotype=genotype
        ))
        case (true, Some(value))      => Seq(error(
          s"FORMAT.PS genotype was phased but had a non-integer phase set value `${value.getClass.getSimpleName}", variant=variant, genotype=genotype
        ))
      }
    }
  }

  case class VariantFormatExtraFieldValidator(formatKeys: Set[String]) extends VariantValidator {
    def validate(variant: Variant): Seq[ValidationResult] = {
      // check for FORMAT values not described in the header
      variant.genotypes
        .flatMap { case (_, genotype) => genotype.attrs.keySet }
        .toSet
        .diff(formatKeys)
        .map { id: String => error(s"FORMAT.$id found in record but missing in header", variant=variant) }
        .toIndexedSeq
    }
  }
}

