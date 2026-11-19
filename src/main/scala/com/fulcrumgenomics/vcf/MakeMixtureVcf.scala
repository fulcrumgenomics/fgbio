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
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.vcf.MakeMixtureVcf.Sample
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api._

import java.text.DecimalFormat
import scala.collection.mutable


object MakeMixtureVcf {
  /** The format field for the allele fraction to be written to. */
  val AlleleFractionField = "AF"
  /** The filter that is placed on records where the full genotype cannot be computed. */
  val UnknownGtFilter     = "unknown_gt"

  /** The FILTER and FORMAT header lines added to the mixture VCF. */
  val HeaderLines: Seq[VcfHeaderEntry] = Seq(
    VcfFilterHeader(UnknownGtFilter, "GT and AF not known due to one or more input samples being no-called."),
    VcfFormatHeader(AlleleFractionField, VcfCount.OnePerAltAllele, VcfFieldType.Float, "Alt allele fraction.")
  )

  /** Case class to store a sample and a proportion in the mixture. */
  private[vcf] case class Sample(name: String, proportion: Double = 0) {
    require(proportion >= 0 && proportion <= 1, "Sample proportions must be between 0 and 1 inclusive.")
    override def toString: String = f"${name}@${proportion}%.4f"
  }

  /**
    * Gets the genotype of the sample, substituting in a HomRef genotype for no-calls if required.
    */
  def genotypeOf(sample: String, variant: Variant, noCallIsHomRef: Boolean): Genotype = {
    val actual = variant.gt(sample)
    if (actual.isNoCall && noCallIsHomRef) Genotype(variant.alleles, sample, IndexedSeq(variant.alleles.ref, variant.alleles.ref))
    else actual
  }

  /** Makes a diploid no-call genotype with the provided sample name. */
  def noCall(sampleName: String, alleles: AlleleSet): Genotype = {
    Genotype(alleles, sampleName, IndexedSeq(NoCallAllele, NoCallAllele))
  }

  /**
    * Subsets the variant to the genotypes of the given samples, and its alleles to the reference allele and the
    * alternate alleles called in those genotypes.
    */
  private[vcf] def subsetToSamples(variant: Variant, samples: Set[String]): Variant = {
    val genotypes = variant.genotypes.filter { case (sample, _) => samples.contains(sample) }
    val called    = genotypes.valuesIterator.flatMap(_.calls).toSet
    val alleles   = AlleleSet(variant.alleles.ref, variant.alleles.alts.filter(called.contains))
    variant.copy(alleles=alleles, genotypes=genotypes.map { case (sample, gt) => sample -> gt.copy(alleles=alleles) })
  }

  /**
    * Builds a writer for a mixture VCF.  The header has the version, contigs, FILTER, FORMAT and INFO lines of the input
    * header along with the given FILTER and FORMAT lines (where the input does not already define them), and the given
    * samples.
    */
  def makeWriter(output: PathToVcf, in: VcfHeader, lines: Seq[VcfHeaderEntry], samples: String*): VcfWriter = {
    val filters = lines.collect { case f: VcfFilterHeader => f }.filterNot(f => in.filter.contains(f.id))
    val formats = lines.collect { case f: VcfFormatHeader => f }.filterNot(f => in.format.contains(f.id))
    val header  = in.copy(
      filters = in.filters ++ filters,
      formats = in.formats ++ formats,
      others  = Seq.empty,
      samples = samples.toIndexedSeq
    )
    VcfWriter(output, header)
  }

  /** Converts the value of an allele fraction field to one or more fractions. */
  private[vcf] def toFractions(value: Any): Seq[Double] = value match {
    case x if Variant.isMissingValue(x) => throw new IllegalArgumentException(s"Missing value in allele fraction: $value")
    case f: Float      => Seq(f.toString.toDouble) // via the decimal text, as the VCF has it
    case d: Double     => Seq(d)
    case i: Int        => Seq(i.toDouble)
    case s: String     => s.split(',').toSeq.map(_.toDouble)
    case xs: Seq[_]    => xs.flatMap(toFractions)
    case _             => throw new IllegalArgumentException(s"Could not parse allele fraction: $value")
  }
}


@clp(group=ClpGroups.VcfOrBcf, description=
  """Creates a VCF with one sample whose genotypes are a mixture of other samples'.
    |
    |The input VCF must contain all samples to be mixed, and may optionally contain other samples.
    |Sample mixtures can be specified in one of several ways:
    |
    |  1. `--samples s1 s2 s3`: specifies that the samples with names `s1`, `s`2 and `s3` should be mixed equally.
    |  2. ` --samples s1@0.1 s2@0.1 s3@0.8`: specifies that the three samples should be mixed at `0.1`, `0.1`, and `0.8`
    |                                        respectively
    |  3. `--samples s1@0.1 s2 s3`: specifies that `s1` should form `0.1` of the mixture and that the remaining (`0.9`)
    |                               should be equally split amongst `s2` and s`3`
    |  4. If no sample names are given, all samples in the input VCF will be used, mixing equally
    |
    |The input samples are assumed to be diploid with the allele fraction defined by the genotype (`0/0`, `0/1` or `1/1`).
    |The `--allele-fraction-field` option may be specified, in which case the allele fraction of the input genotypes
    |will be retrieved from the specified `FORMAT` field in the VCF genotypes.  All genotypes except hom-ref and
    |no-call genotypes must have this field present if supplied.
    |
    |If the --`no-call-is-hom-ref` flag is `true` (the default) then no-call genotypes in the input VCF are interpreted
    |as hom-ref genotypes.  If it is `false`, any location with a no-call genotype will be emitted as a no-call in
    |the mixture, and will be filtered.
  """)
class MakeMixtureVcf
( @arg(flag='i', doc="Input VCF containing genotypes for the samples to be mixed.") val input: PathToVcf,
  @arg(flag='o', doc="Output VCF of mixture sample.") val output: PathToVcf,
  @arg(flag='s', doc="Samples to mix. See general usage for format and examples.", minElements=0) val samples: Seq[String] = Seq.empty,
  @arg(flag='S', doc="Output sample name.") val outputSampleName: String = "mixture",
  @arg(flag='N', doc="Treat no-calls for samples as hom-ref genotypes.") val noCallIsHomRef: Boolean = true,
  @arg(flag='a', doc="Format field containing allele fraction.") val alleleFractionField: Option[String] = None,
  @arg(flag='p', doc="Digits of precision in generated allele fractions.") val precision: Int = 5
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  /**
    * If no sample names are provided creates Sample objects for all the samples in the VCF
    * in equal proportion.  If sample names are provided, creates Sample objects for each,
    * auto-sets the proportion for any that do not have it specified, then validates that
    * all specified samples exist in the VCF.
    *
    * @param header the input VCF header
    * @return the set of samples with proportions
    */
  private[vcf] def determineSamples(header: VcfHeader): Seq[Sample] = this.samples match {
    case Seq() =>
      val fraction = 1.0 / header.samples.size.toDouble
      header.samples.map(s => Sample(s, fraction))
    case sampleStrings =>
      // Extract the samples and fractions from their strings
      var samples = sampleStrings.map { s =>
        s.lastIndexOf('@') match {
          case -1 => Sample(s)
          case  i => Sample(s.substring(0, i), s.substring(i+1).toDouble)
        }
      }

      // Validate that each sample was specified only once
      require(samples.size == samples.map(_.name).distinct.size, "Each sample name can be specified at most once.")

      // Auto-set fraction for any samples that didn't have it
      if (samples.exists(_.proportion == 0)) {
        val total = samples.map(_.proportion).sum
        require(total <= 1, "Sample fractions added to more than 1.0")
        val remainder = 1 - total
        val fraction  = remainder / samples.count(_.proportion == 0)
        samples = samples.map(s => if (s.proportion == 0) s.copy(proportion = fraction) else s)
      }

      // Validate that the proportions add up to ~1
      val total = samples.map(_.proportion).sum
      require(total >= 0.999 && total <= 1.001, s"Sample fractions must add to 1.0. Total: $total")

      // Then lastly validate that all given samples exist in the VCF
      val missing = samples.filterNot(s => header.samples.contains(s.name))
      require(missing.isEmpty, s"Samples are not present in the VCF: ${missing.mkString(", ")}")

      samples
  }

  /** The main execute method that reads the inputs, mixes and writes the outputs. */
  override def execute(): Unit = {
    val in            = VcfSource(input, allowMissingFieldsInHeader=true)
    val samples       = determineSamples(in.header)
    val sampleNameSet = samples.map(_.name).toSet
    val out           = MakeMixtureVcf.makeWriter(output, in.header, MakeMixtureVcf.HeaderLines, this.outputSampleName)
    val fmt           = new DecimalFormat("0." + ("#" * precision))

    logger.info("Making mixture of: ", samples.mkString(", "))
    val progress = new ProgressLogger(logger, noun="variants", unit=50000)

    in.iterator.map(MakeMixtureVcf.subsetToSamples(_, sampleNameSet)).foreach { variant =>
      val outputGt = {
        if (!this.noCallIsHomRef && variant.gts.exists(_.isNoCall)) {
          MakeMixtureVcf.noCall(this.outputSampleName, variant.alleles)
        }
        else {
          val alleleFractions = mutable.Map[Allele, Double](variant.alleles.toSeq.map(a => a -> 0d):_*)
          samples.foreach(updateAlleleFractionsForSample(variant, _, alleleFractions))

          // The allele fractions are written as a string instead of numbers so that they are written with the
          // requested precision, rather than the three decimal places that htsjdk's VCF code would round them to.
          val afString = variant.alleles.alts.map(a => fmt.format(alleleFractions(a))).mkString(",")
          val calls    = variant.alleles.iterator.filter(a => alleleFractions(a) > 0).toIndexedSeq
          Genotype(variant.alleles, this.outputSampleName, calls, attrs=Map(MakeMixtureVcf.AlleleFractionField -> afString))
        }
      }

      // Finally make the variant, keeping the input filters and adding one if the genotype is unknown
      val filters = variant.filters.filterNot(Variant.PassingFilters.contains) ++
        (if (outputGt.isNoCall) Some(MakeMixtureVcf.UnknownGtFilter) else None)
      out.write(Variant(
        chrom     = variant.chrom,
        pos       = variant.pos,
        id        = variant.id,
        alleles   = variant.alleles,
        filters   = if (filters.isEmpty) Variant.PassingFilters else filters,
        genotypes = Map(this.outputSampleName -> outputGt)
      ))
      progress.record(variant.chrom, variant.pos)
    }

    in.safelyClose()
    out.close()
  }

  /**
    * Updates a map of allele fractions based on the contribution from the provided sample.
    * This method should not be called when one or more samples is no-called and noCallIsHomRef is false.
    */
  private[vcf] def updateAlleleFractionsForSample(variant: Variant, sample: Sample, alleleFractions: mutable.Map[Allele,Double]): Unit = {
    val gt  = MakeMixtureVcf.genotypeOf(sample.name, variant, this.noCallIsHomRef)
    val ref = variant.alleles.ref
    require(!gt.isNoCall, "Bug: method should not be called with no-call GTs when noCallIsHomRef == false.")

    this.alleleFractionField match {
      case None =>
        val count = gt.calls.length.toDouble
        gt.calls.foreach(a => alleleFractions(a) += sample.proportion / count)

      case Some(afField) =>
        getAfs(gt, afField) match {
          case Seq()   =>
            // Missing AF field is only allowed if the genotype is hom-ref
            if (gt.isHomRef) alleleFractions(ref) += sample.proportion
            else fail(s"Sample ${sample.name} missing format field ${afField} at ${variant.chrom}:${variant.pos}")
          case Seq(d)  =>
            val pos = s"${variant.chrom}:${variant.pos}"
            if (gt.isHetNonRef) fail(s"Only one value for AF for sample ${sample.name} at $pos with het-non-ref genotype.")
            else if (gt.isHomRef) alleleFractions(ref) += sample.proportion
            else if (gt.isHet)    gt.calls.foreach(a => alleleFractions(a) += sample.proportion * (if (a == ref) 1-d else d))
            else if (gt.isHomVar) {
              alleleFractions(ref)          += sample.proportion * (1-d)
              alleleFractions(gt.calls(0))  += sample.proportion * d
            }
            else fail(s"Got a single value for AF for sample ${sample.name} at $pos with genotype ${gt.gtWithBases}")
          case ds =>
            // If we got back an array of length > 1 then it _must_ have an entry for every non-ref allele
            // on the variant, not just those in the genotype
            require(ds.length == variant.alleles.alts.size,
              s"Variant at ${variant.chrom}:${variant.pos} had ${variant.alleles.alts.size} alt alleles but sample " +
                s"${sample.name} had an ${afField} with ${ds.length} entries.")

            ds.zip(variant.alleles.alts).foreach { case (af, allele) =>
              alleleFractions(allele) += sample.proportion * af
            }

            alleleFractions(ref) += sample.proportion * (1 - ds.sum)
        }
    }
  }

  /**
    * Retrieves the requested field as a sequence of fractions, empty if the field is not present.
    */
  private[vcf] def getAfs(gt: Genotype, afField: String): Seq[Double] = gt.get[Any](afField) match {
    case None        => Seq.empty
    case Some(value) => MakeMixtureVcf.toFractions(value)
  }
}
