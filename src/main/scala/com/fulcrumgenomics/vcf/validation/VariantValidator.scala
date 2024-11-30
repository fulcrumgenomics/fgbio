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

import com.fulcrumgenomics.vcf.api.VcfHeader._
import com.fulcrumgenomics.vcf.api._
import com.fulcrumgenomics.vcf.validation.ValidationResult._

import scala.collection.immutable.ArraySeq
import scala.collection.mutable.ListBuffer

trait VariantValidator extends VcfValidator {
  import Variant.{FlagValue, Missing, MissingChar, MissingFloat, MissingInt}

  def validate(variant: Variant): Seq[ValidationResult]

  private def toKind(value: Any): Seq[VcfFieldType] = value match {
    case arr: ArrayAttr[_]               => arr.flatMap(toKind)
    case arr: ArraySeq[_]                => arr.flatMap(toKind)
    case MissingInt                      => Seq(VcfFieldType.Integer)
    case _ if value.isInstanceOf[Int]    => Seq(VcfFieldType.Integer)
    case MissingFloat                    => Seq(VcfFieldType.Float)
    case _ if value.isInstanceOf[Float]  => Seq(VcfFieldType.Float)
    case Missing                         => Seq(VcfFieldType.String)
    case _ if value.isInstanceOf[String] => Seq(VcfFieldType.String)
    case MissingChar                     => Seq(VcfFieldType.Character)
    case _ if value.isInstanceOf[Char]   => Seq(VcfFieldType.Character)
    case FlagValue                       => Seq(VcfFieldType.Flag) // Note: could also be VcfFieldType.Fixed(0)
  }

  protected def validate(kind: VcfFieldType,
                         count: VcfCount,
                         variant: Variant,
                         source: String,
                         value: Any,
                         genotype: Option[Genotype] = None): Seq[ValidationResult] = {
    val builder = new ListBuffer[ValidationResult]()
    val actualCount: Int = value match {
      case _ if value == Variant.FlagValue => 0
      case arr: ArrayAttr[_]               => arr.length
      case arr: ArraySeq[_]                => arr.length
      case seq: Seq[_]                     => seq.length // should not happen, but for good measure
      case _                               => 1
    }
    val expectedCount: Option[Int] = count match {
      case VcfCount.OnePerAltAllele => Some(variant.alleles.alts.length)
      case VcfCount.OnePerAllele    => Some(variant.alleles.size)
      case VcfCount.OnePerGenotype  => Some(1)
      case VcfCount.Unknown         => None
      case VcfCount.Fixed(n)        => Some(n)
    }
    if (expectedCount.exists(_ != actualCount)) {
      val _expectedCount = expectedCount.getOrElse(0)
      builder.append(error(
        s"$source expected `${_expectedCount}` values, found `$actualCount`", variant=Some(variant), genotype=genotype
      ))
    }

    val actualKind: Seq[VcfFieldType] = toKind(value=value)
    val kindOk: Boolean               = {
      if (actualCount == 0 && actualKind.forall(_ == VcfFieldType.Flag)) {
        // can also be VcfFieldType.Fixed(0)
        kind == VcfFieldType.Flag || count == VcfCount.Fixed(0)
      }
      else {
        actualKind.forall(_ == kind)
      }
    }
    if (!kindOk) {
      val _actualKind = actualKind match {
        case Seq(_kind) => _kind
        case kinds      => kinds.distinct.mkString(",")
      }
      builder.append(error(
        s"$source expected `$kind` kind, found `${_actualKind}`", variant=Some(variant), genotype=genotype
      ))
    }

    if (builder.isEmpty) IndexedSeq.empty[ValidationResult] else builder.toIndexedSeq
  }
}



object VariantValidator {
  val VariantInfoValidators  : Seq[VariantInfoValidator]   = ReservedVcfInfoHeaders.map(VariantInfoValidator.apply)

  case class VariantInfoValidator(info: VcfInfoHeader) extends VariantValidator {
    def validate(variant: Variant): Seq[ValidationResult] = {
      variant.get[Any](info.id).toIndexedSeq.flatMap { value =>
        this.validate(kind=info.kind, count=info.count, variant=variant, source=f"INFO.${info.id}", value=value)
      }
    }
  }

  case class VariantInfoExtraFieldValidator(infoKeys: Set[String]) extends VariantValidator {
    def validate(variant: Variant): Seq[ValidationResult] = {
      variant.attrs.keys
        .filterNot(infoKeys.contains)
        .map { id =>
          error(s"INFO.$id found in record but missing in header", variant = variant)
        }.toIndexedSeq
    }
  }
}

// TODO: validators across variants (ex. spanning alleles), or across genotypes (phase set)
sealed trait VariantsValidator {
  def validate(variant: Variant): Seq[ValidationResult]
  def finish(): Seq[ValidationResult]
}