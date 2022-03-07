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

import java.text.DecimalFormat
import java.util
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.vcf.MakeMixtureVcf.Sample
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.util.CollectionUtil
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.mutable
import scala.jdk.CollectionConverters._


object MakeMixtureVcf {
  /** The format field for the allele fraction to be written to. */
  val AlleleFractionField = "AF"
  /** The filter that is placed on records where the full genotype cannot be computed. */
  val UnknownGtFilter     = "unknown_gt"

  val HeaderLines = Seq[VCFHeaderLine](
    new VCFFilterHeaderLine(UnknownGtFilter, "GT and AF not known due to one or more input samples being no-called."),
    new VCFFormatHeaderLine(AlleleFractionField, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Alt allele fraction.")
  )

  /** Case class to store a sample and a proportion in the mixture. */
  private[vcf] case class Sample(name: String, proportion: Double = 0) {
    require(proportion >= 0 && proportion <= 1, "Sample proportions must be between 0 and 1 inclusive.")
    override def toString: String = f"${name}@${proportion}%.4f"
  }

  /**
    * Gets the genotype of the sample, substituting in a HomRef genotype for no-calls if required.
    */
  def genotypeOf(sample: String, ctx: VariantContext, noCallIsHomRef: Boolean): Genotype = {
    val actual = ctx.getGenotype(sample)
    if (actual.isNoCall && noCallIsHomRef) new GenotypeBuilder(sample, util.Arrays.asList(ctx.getReference, ctx.getReference)).make()
    else actual
  }

  /** Makes a no-call genotype object with the provided sample name. */
  def noCall(sampleName: String): Genotype = {
    new GenotypeBuilder(sampleName).alleles(util.Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).make()
  }

  /** Builds a header that can be used to write a mixture VCF. */
  def makeWriter(output: PathToVcf, in: VCFHeader, lines: Seq[VCFHeaderLine], samples: String*): VariantContextWriter = {
    val header = new VCFHeader(Collections.emptySet[VCFHeaderLine](), util.Arrays.asList(samples:_*))

    in.getFilterLines.foreach(header.addMetaDataLine)
    in.getFormatHeaderLines.foreach(header.addMetaDataLine)
    in.getInfoHeaderLines.foreach(header.addMetaDataLine)
    lines.foreach(header.addMetaDataLine)

    header.setSequenceDictionary(in.getSequenceDictionary)

    val writer = new VariantContextWriterBuilder()
      .setOutputFile(output.toFile)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setReferenceDictionary(in.getSequenceDictionary)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .build()

    writer.writeHeader(header)
    writer
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
  private[vcf] def determineSamples(header: VCFHeader): Seq[Sample] = this.samples match {
    case Seq() =>
      val ss       = header.getSampleNamesInOrder
      val fraction = 1.0 / ss.size().toDouble
      header.getSampleNamesInOrder.map(s => Sample(s, fraction)).toSeq
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
      val missing = samples.filterNot(s => header.getSampleNamesInOrder.contains(s.name))
      require(missing.isEmpty, s"Samples are not present in the VCF: ${missing.mkString(", ")}")

      samples
  }

  /** The main execute method that reads the inputs, mixes and writes the outputs. */
  override def execute(): Unit = {
    val in            = new VCFFileReader(input.toFile, false)
    val samples       = determineSamples(in.getFileHeader)
    val sampleNameSet = CollectionUtil.makeSet(samples.map(_.name):_*)
    val out           = MakeMixtureVcf.makeWriter(output, in.getFileHeader, MakeMixtureVcf.HeaderLines, this.outputSampleName)
    val fmt           = new DecimalFormat("0." + ("#" * precision))

    logger.info("Making mixture of: ", samples.mkString(", "))
    val progress = new ProgressLogger(logger, noun="variants", unit=50000)

    in.map(_.subContextFromSamples(sampleNameSet)).foreach { ctx =>
      val outputGt = {
        if (!this.noCallIsHomRef && ctx.getNoCallCount > 0) {
          MakeMixtureVcf.noCall(this.outputSampleName)
        }
        else {
          val alleleFractions = mutable.Map[Allele, Double](ctx.getAlleles.toSeq.map(a => a -> 0d):_*)
          samples.foreach(updateAlleleFractionsForSample(ctx, _, alleleFractions))

          // This horror (turning the AF into a string instead of a double[]) brought to you by the
          // fact that htsjdk's VCF code will otherwise format the doubles by rounding to 3 decimal places!
          val afString = ctx.getAlternateAlleles.map(a => fmt.format(alleleFractions(a))).mkString(",")
          val alleles = alleleFractions.filter(_._2 > 0).keys.toSeq.sortBy(ctx.getAlleleIndex)
          val builder = new GenotypeBuilder(this.outputSampleName)
          builder.alleles(CollectionUtil.makeList(alleles:_*))
          builder.attribute(MakeMixtureVcf.AlleleFractionField, afString)
          builder.make()
        }
      }

      // Finally make the variant context and genotype objects
      val builder = new VariantContextBuilder(ctx.getSource, ctx.getContig, ctx.getStart(), ctx.getEnd(), ctx.getAlleles)
      builder.id(ctx.getID)
      builder.genotypes(outputGt)
      builder.filters(new util.HashSet[String](ctx.getFilters))
      if (outputGt.isNoCall) builder.filter(MakeMixtureVcf.UnknownGtFilter)
      out.add(builder.make())
      progress.record(ctx.getContig, ctx.getStart)
    }

    out.close()
  }

  /**
    * Updates a map of allele fractions based on the contribution from the provided sample.
    * This method should not be called when one or more samples is no-called and noCallIsHomRef is false.
    */
  private[vcf] def updateAlleleFractionsForSample(ctx: VariantContext, sample: Sample, alleleFractions: mutable.Map[Allele,Double]): Unit = {
    val gt = MakeMixtureVcf.genotypeOf(sample.name, ctx, this.noCallIsHomRef)
    require(!gt.isNoCall, "Bug: method should not be called with no-call GTs when noCallIsHomRef == false.")

    this.alleleFractionField match {
      case None =>
        val alleles = gt.getAlleles
        val count   = alleles.size().toDouble
        alleles.foreach(a => alleleFractions(a) += sample.proportion / count)

      case Some(afField) =>
        getAfs(gt, afField) match {
          case Array()   =>
            // Missing AF field is only allowed if the genotype is hom-ref
            if (gt.isHomRef) alleleFractions(ctx.getReference) += sample.proportion
            else fail(s"Sample ${sample.name} missing format field ${afField} at ${ctx.getContig}:${ctx.getStart}")
          case Array(d)  =>
            val pos = s"${ctx.getContig}:${ctx.getStart}"
            if (gt.isHetNonRef) fail(s"Only one value for AF for sample ${sample.name} at $pos with het-non-ref genotype.")
            else if (gt.isHomRef) alleleFractions(ctx.getReference) += sample.proportion
            else if (gt.isHet)    gt.getAlleles.foreach(a => alleleFractions(a) += sample.proportion * (if (a.isReference) 1-d else d))
            else if (gt.isHomVar) {
              alleleFractions(ctx.getReference) += sample.proportion * (1-d)
              alleleFractions(gt.getAllele(0))  += sample.proportion * d
            }
            else fail(s"Got a single value for AF for sample ${sample.name} at $pos with genotype ${gt.getGenotypeString}")
          case ds =>
            // If we got back an array of length > 1 then it _must_ have an entry for every non-ref allele
            // on the variant context, not just those in the genotype
            require(ds.length == ctx.getAlternateAlleles.size,
              s"Variant at ${ctx.getContig}:${ctx.getStart} had ${ctx.getNAlleles-1} alt alleles but sample " +
                s"${sample.name} had an ${afField} with ${ds.length} entries.")

            ds.zip(ctx.getAlternateAlleles.asScala).foreach { case (af, allele) =>
              alleleFractions(allele) += sample.proportion * af
            }

            alleleFractions(ctx.getReference) += sample.proportion * (1 - ds.sum)
        }
    }
  }

  /**
    * Retrieves the requested field as an array of doubles.
    */
  private[vcf] def getAfs(gt: Genotype, afField: String): Array[Double] = gt.getExtendedAttribute(afField) match {
    case null      => Array[Double]()
    case s: String => s.split(',').map(_.toDouble)
  }
}
