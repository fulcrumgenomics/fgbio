/*
 * The MIT License
 *
 * Copyright (c) 2026 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.vcf.api._

/**
  * Helpers shared by fgbio's VCF tools for working with [[Variant]]s: the standard FORMAT header lines that tools add to
  * the VCFs they write, subsetting a variant to samples, and coercing INFO and FORMAT values (e.g. of fields not
  * described in the header) to numbers.  Predicates on a single variant (e.g. [[Variant.isSnp]]) are on [[Variant]].
  */
private[fulcrumgenomics] object VariantUtil {
  /** The definitions of the standard FORMAT fields (as HTSJDK defines them) that tools add to the VCFs they write. */
  val StandardFormatHeaderLines: Seq[VcfFormatHeader] = Seq(
    VcfFormatHeader("GT", VcfCount.Fixed(1),       VcfFieldType.String,  "Genotype"),
    VcfFormatHeader("GQ", VcfCount.Fixed(1),       VcfFieldType.Integer, "Genotype Quality"),
    VcfFormatHeader("DP", VcfCount.Fixed(1),       VcfFieldType.Integer, "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"),
    VcfFormatHeader("AD", VcfCount.OnePerAllele,   VcfFieldType.Integer, "Allelic depths for the ref and alt alleles in the order listed"),
    VcfFormatHeader("PL", VcfCount.OnePerGenotype, VcfFieldType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"),
    VcfFormatHeader("FT", VcfCount.Unknown,        VcfFieldType.String,  "Genotype-level filter")
  )

  /**
    * Subsets the variant to the genotypes of the given samples, and its alleles to the reference allele and the
    * alternate alleles called in those genotypes, as HTSJDK's `VariantContext.subContextFromSamples` does.  Other
    * per-allele values (e.g. AD) are not subset.  Returns the variant itself if it has no other samples and all of its
    * alternate alleles are called.
    */
  def subsetToSamples(variant: Variant, samples: Set[String]): Variant = {
    val genotypes = if (variant.genotypes.keysIterator.forall(samples.contains)) variant.genotypes else {
      variant.genotypes.filter { case (sample, _) => samples.contains(sample) }
    }
    val called = genotypes.valuesIterator.flatMap(_.calls).toSet
    if ((genotypes eq variant.genotypes) && variant.alleles.alts.forall(called.contains)) variant else {
      val alleles = AlleleSet(variant.alleles.ref, variant.alleles.alts.filter(called.contains))
      variant.copy(alleles=alleles, genotypes=genotypes.map { case (sample, gt) => sample -> gt.copy(alleles=alleles) })
    }
  }

  /**
    * Converts the value of a numeric INFO or FORMAT field to one or more Doubles.  A Float is converted through its
    * shortest decimal text (e.g. `0.05f` to `0.05`, not `0.05000000074505806`), and a String (e.g. for a field not
    * described in the header) is split on commas.
    *
    * @param value the value of the field
    * @param context a description of the field for error messages, e.g. `AF of sample s1 at chr1:100`
    * @throws IllegalArgumentException if the value is missing or not numeric
    */
  def toDoubles(value: Any, context: => String): IndexedSeq[Double] = value match {
    case x if Variant.isMissingValue(x) => throw new IllegalArgumentException(s"Missing value for $context")
    case f: Float                       => IndexedSeq(f.toString.toDouble)
    case d: Double                      => IndexedSeq(d)
    case i: Int                         => IndexedSeq(i.toDouble)
    case s: String                      => s.split(',').iterator.map(x => parse(x, context)(_.toDouble)).toIndexedSeq
    case xs: Iterable[_]                => xs.iterator.flatMap(toDoubles(_, context)).toIndexedSeq
    case _                              => throw new IllegalArgumentException(s"Not a numeric value for $context: $value")
  }

  /**
    * Converts the value of an Integer INFO or FORMAT field that holds a single value to an Int.  A String (e.g. for a
    * field not described in the header) is parsed, and a one-element list (as a field not described in the header may
    * be read) is unwrapped.
    *
    * @param value the value of the field
    * @param context a description of the field for error messages, e.g. `PS at chr1:100`
    * @throws IllegalArgumentException if the value is missing, not an integer, or holds more than one value
    */
  def toInt(value: Any, context: => String): Int = value match {
    case x if Variant.isMissingValue(x) => throw new IllegalArgumentException(s"Missing value for $context")
    case i: Int                         => i
    case s: String if !s.contains(',')  => parse(s, context)(_.toInt)
    case Seq(x)                         => toInt(x, context)
    case xs: Iterable[_] if xs.size > 1 => throw new IllegalArgumentException(s"Expected one value for $context but found: ${xs.mkString(",")}")
    case s: String                      => throw new IllegalArgumentException(s"Expected one value for $context but found: $s")
    case _                              => throw new IllegalArgumentException(s"Not an integer value for $context: $value")
  }

  /** Parses a String with `f`, reporting a value that cannot be parsed as an [[IllegalArgumentException]] naming the
    * context. */
  private def parse[A](s: String, context: => String)(f: String => A): A = {
    try f(s) catch {
      case _: NumberFormatException => throw new IllegalArgumentException(s"Not a numeric value for $context: $s")
    }
  }
}
