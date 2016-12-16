/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

import java.util
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.MakeTwoSampleMixtureVcf._
import dagr.sopt
import dagr.sopt.{arg, clp}
import htsjdk.samtools.util.IntervalList
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.JavaConversions.asJavaCollection
import scala.collection.JavaConversions.iterableAsScalaIterable
import scala.collection.mutable.ListBuffer

object MakeTwoSampleMixtureVcf {
  val AlleleFractionField = "AF"
  val GermlineFilter      = "alt_allele_in_normal"
  val UnknownGtFilter     = "unknown_gt"
  val MultiAllelicFilter  = "multi_allelic"

  val TumorSampleName  = "tumor"
  val NormalSampleName = "normal"
}

/**
  * Creates a VCF by mixing two germline samples at a given proportion.
  *
  * @author Tim Fennell
  */
@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Creates a simulated tumor or tumor/normal VCF by in-silico mixing genotypes from two samples.
    |The tumor genotypes are created by mixing the incoming genotypes for the for the
    |'tumor' sample and the incoming genotypes for the 'normal' samples with the 'tumor'
    |alleles accounting for 'tumorFraction' of the resulting mixture and the 'normal' alleles
    |accounting for '1 - tumorFraction' of the resulting mixture.  E.g. if the 'tumor' genotype
    |is A/C and the 'normal' genotype is C/C and 'tumorFraction' is set at 0.5 the resulting
    |tumor genotype will be A/C with an allele fraction of 0.75.  The resulting allele fraction
    |is written to the AF info field in the VCF.
    |
    |In tumor-only mode only tumor genotypes are output.  In tumor/normal mode genotypes for
    |the 'normal' samples are also emitted, and match the genotypes from the input sample.
    |
    |All loci (potentially restricted by intervals) that are variant in one or both samples are written
    |to the output VCF, though in several cases the variants will be filtered:
    |  - If either of the tumor or normal sample is no-called the resulting locus will have the
    |    'unknown_gt' filter applied
    |  - If the tumor and the normal have more than one alternative allele between them the
    |    'multi_allelic' filter will be applied
    |  - In tumor/normal mode (as opposed to tumor-only) loci that have an alt allele in the normal
    |    sample will have the 'alt_allele_in_normal' filter applied
  """
)
class MakeTwoSampleMixtureVcf
( @arg(flag="i", doc="Input VCF file.")  val input: PathToVcf,
  @arg(flag="o", doc="Output VCF file.") val output: PathToVcf,
  @arg(flag="t", doc="Name of the 'tumor' sample in the input VCF.")  val tumor: String,
  @arg(flag="n", doc="Name of the 'normal' sample in the input VCF.") val normal: String,
  @arg(flag="f", doc="What fraction of the mixture comes from the 'tumor' sample.") val tumorFraction: Double = 0.5,
  @arg(flag="T", doc="Tumor only mode - only output tumor genotypes and don't filter sites.") val tumorOnly: Boolean = false,
  @arg(flag="N", doc="Treat no-calls for either sample as hom-ref genotypes.") val noCallIsHomRef: Boolean = true,
  @arg(flag="l", doc="Optional set of intervals to restrict to.") val intervals: Option[PathToIntervals] = None
) extends FgBioTool {
  if (tumorFraction < 0 || tumorFraction > 1) fail("Tumor fraction must be between 0 and 1.")
  if (tumor == normal) fail("Tumor and Normal samples must be different samples.")

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertCanWriteFile(_))

  override def execute(): Unit = {
    val in = new VCFFileReader(input.toFile, true)
    if (!in.getFileHeader.getSampleNamesInOrder.contains(tumor))  fail(s"VCF does not contain tumor sample ${tumor}")
    if (!in.getFileHeader.getSampleNamesInOrder.contains(normal)) fail(s"VCF does not contain normal sample ${normal}")

    val out = if (tumorOnly) makeWriter(in.getFileHeader, "tumor") else makeWriter(in.getFileHeader, "tumor", "normal")
    iterator(in, this.intervals).foreach { ctx =>
      val oldTumorGt = genotypeOf(tumor, ctx)
      val oldNormalGt = genotypeOf(normal, ctx)
      val oldGts = Seq(oldTumorGt, oldNormalGt)

      val alleles = (oldTumorGt.getAlleles ++ oldNormalGt.getAlleles ++ Some(ctx.getReference)).filterNot(_ == Allele.NO_CALL).toSet
      val isMultiAllelic = alleles.size > 2
      val isMonomorphic = alleles.size < 2

      if (!isMonomorphic) {
        val (newTumorGt, newNormalGt) = if (isMultiAllelic) {
          (noCall(TumorSampleName), if (tumorOnly) None else Some(noCall(NormalSampleName)))
        }
        else {
          val refAllele = ctx.getReference
          val altAllele = ctx.getAlternateAllele(0)

          // Build the Tumor genotype
          val builder = new GenotypeBuilder(TumorSampleName)
          val alleles = new util.ArrayList[Allele]()
          if (oldGts.forall(_.isCalled)) {
            if (oldGts.map(_.countAllele(refAllele)).sum > 0) alleles.add(refAllele)
            if (oldGts.map(_.countAllele(altAllele)).sum > 0) alleles.add(altAllele)
            if (alleles.size() == 1) alleles.add(alleles.get(0))
            builder.alleles(alleles)

            // We count alts only, because ref and no-call can be treated as ref here
            val af = (this.tumorFraction * oldTumorGt.countAllele(altAllele) / 2) + ((1 - this.tumorFraction) * oldNormalGt.countAllele(altAllele) / 2)
            builder.attribute(MakeTwoSampleMixtureVcf.AlleleFractionField, af)
          }
          else {
            builder.alleles(util.Arrays.asList(Allele.NO_CALL, Allele.NO_CALL))
          }

          (builder.make(), if (tumorOnly) None else Some(new GenotypeBuilder("normal", oldNormalGt.getAlleles).make()))
        }

        val gts = Seq(Some(newTumorGt), newNormalGt).flatten
        val builder = new VariantContextBuilder(ctx.getSource, ctx.getContig, ctx.getStart(), ctx.getEnd(), ctx.getAlleles)
        builder.genotypes(gts)

        val filters = ListBuffer[String]()
        if (isMultiAllelic) filters += MultiAllelicFilter
        if (!tumorOnly && !oldNormalGt.isHomRef) filters += GermlineFilter
        if (newTumorGt.isNoCall) filters += UnknownGtFilter
        if (filters.nonEmpty) builder.filters(filters: _*) else builder.passFilters()
        out.add(builder.make())
      }
    }

    out.close()
    in.safelyClose()
  }

  /** Gets the genotype of the sample, substituting in a HomRef genotype for no-calls if required. */
  private def genotypeOf(sample: String, ctx: VariantContext): Genotype = {
    val actual = ctx.getGenotype(sample)
    if (actual.isNoCall && this.noCallIsHomRef) new GenotypeBuilder(sample, util.Arrays.asList(ctx.getReference)).make()
    else actual
  }

  /** Makes a no-call genotype object with the provided sample name. */
  private def noCall(sampleName: String): Genotype = {
    new GenotypeBuilder(sampleName).alleles(util.Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).make()
  }

  /** Builds a header that can be used to write the mixture VCF. */
  def makeWriter(in: VCFHeader, samples: String*): VariantContextWriter = {
    val header = new VCFHeader(Collections.emptySet[VCFHeaderLine](), util.Arrays.asList(samples:_*))

    in.getFilterLines.foreach(header.addMetaDataLine)
    header.addMetaDataLine(new VCFFilterHeaderLine(GermlineFilter,  "Evidence seen in the normal sample"))
    header.addMetaDataLine(new VCFFilterHeaderLine(UnknownGtFilter, "GT and AF not known due to one or more input samples being no-called."))
    header.addMetaDataLine(new VCFFilterHeaderLine(MultiAllelicFilter, "Locus is multi-allelic in chosen samples."))

    in.getFormatHeaderLines.foreach(header.addMetaDataLine)
    header.addMetaDataLine(new VCFFormatHeaderLine("AF", 1, VCFHeaderLineType.Float, "Alt allele fraction in the tumor"))

    in.getInfoHeaderLines.foreach(header.addMetaDataLine)
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

  /** Generates an iterator over the whole file, or over the intervals if provided. */
  def iterator(in: VCFFileReader, intervals: Option[PathToIntervals]): Iterator[VariantContext] = {
    import scala.collection.JavaConversions.asScalaIterator
    intervals match {
      case None       => in.iterator()
      case Some(path) => ByIntervalListVariantContextIterator(in, IntervalList.fromFile(path.toFile))
    }
  }
}
