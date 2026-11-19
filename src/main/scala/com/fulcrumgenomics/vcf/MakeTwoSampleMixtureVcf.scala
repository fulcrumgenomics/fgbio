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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.MakeTwoSampleMixtureVcf._
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api._
import htsjdk.samtools.util.IntervalList

import scala.collection.mutable.ListBuffer

object MakeTwoSampleMixtureVcf {
  val AlleleFractionField = "AF"
  val GermlineFilter      = "alt_allele_in_normal"
  val UnknownGtFilter     = "unknown_gt"
  val MultiAllelicFilter  = "multi_allelic"

  val TumorSampleName  = "tumor"
  val NormalSampleName = "normal"

  /** The FILTER and FORMAT header lines added to the mixture VCF. */
  val HeaderLines: Seq[VcfHeaderEntry] = MakeMixtureVcf.HeaderLines ++ Seq(
    VcfFilterHeader(GermlineFilter,  "Evidence seen in the normal sample"),
    VcfFilterHeader(MultiAllelicFilter, "Locus is multi-allelic in chosen samples.")
  )
}

/**
  * Creates a VCF by mixing two germline samples at a given proportion.
  *
  * @author Tim Fennell
  */
@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Creates a simulated tumor or tumor/normal VCF by in-silico mixing genotypes from two samples.
    |
    |The tumor genotypes are created by mixing the incoming genotypes for the for the
    |'tumor' sample and the incoming genotypes for the 'normal' samples with the 'tumor'
    |alleles accounting for `tumorFraction` of the resulting mixture and the 'normal' alleles
    |accounting for `1 - tumorFraction` of the resulting mixture (see `--tumor-fraction`).  E.g. if the 'tumor' genotype
    |is `A/C`, the 'normal' genotype is `C/C`, and 'tumorFraction' is set at 0.5, then the resulting
    |tumor genotype will be `A/C` with an allele fraction of `0.75`.  The resulting allele fraction
    |is written to the `AF` info field in the VCF.
    |
    |In tumor-only mode (`--tumor-only`) only tumor genotypes are output.  In tumor/normal mode genotypes for
    |the 'normal' samples are also emitted, and match the genotypes from the input sample.
    |
    |All loci (potentially restricted by intervals) that are variant in one or both samples are written
    |to the output VCF, though in several cases the variants will be filtered:
    |
    |  * If either of the tumor or normal sample is no-called the resulting locus will have the
    |    `unknown_gt` filter applied
    |  * If the tumor and the normal have more than one alternative allele between them the
    |    `multi_allelic` filter will be applied
    |  * In tumor/normal mode (as opposed to tumor-only) loci that have an alt allele in the normal
    |    sample will have the `alt_allele_in_normal` filter applied
  """
)
class MakeTwoSampleMixtureVcf
( @arg(flag='i', doc="Input VCF file.")  val input: PathToVcf,
  @arg(flag='o', doc="Output VCF file.") val output: PathToVcf,
  @arg(flag='t', doc="Name of the 'tumor' sample in the input VCF.")  val tumor: String,
  @arg(flag='n', doc="Name of the 'normal' sample in the input VCF.") val normal: String,
  @arg(flag='f', doc="What fraction of the mixture comes from the 'tumor' sample.") val tumorFraction: Double = 0.5,
  @arg(flag='T', doc="Tumor only mode - only output tumor genotypes and don't filter sites.") val tumorOnly: Boolean = false,
  @arg(flag='N', doc="Treat no-calls for either sample as hom-ref genotypes.") val noCallIsHomRef: Boolean = true,
  @arg(flag='l', doc="Optional set of intervals to restrict to.") val intervals: Option[PathToIntervals] = None
) extends FgBioTool {
  if (tumorFraction < 0 || tumorFraction > 1) fail("Tumor fraction must be between 0 and 1.")
  if (tumor == normal) fail("Tumor and Normal samples must be different samples.")

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertCanWriteFile(_))

  override def execute(): Unit = {
    val in = VcfSource(input, allowMissingFieldsInHeader=true)
    if (!in.header.samples.contains(tumor))  fail(s"VCF does not contain tumor sample ${tumor}")
    if (!in.header.samples.contains(normal)) fail(s"VCF does not contain normal sample ${normal}")

    val outputSamples = if (tumorOnly) Seq(TumorSampleName) else Seq(TumorSampleName, NormalSampleName)
    val out = MakeMixtureVcf.makeWriter(output, in.header, MakeTwoSampleMixtureVcf.HeaderLines, outputSamples:_*)
    iterator(in, this.intervals).foreach { variant =>
      val oldTumorGt  = MakeMixtureVcf.genotypeOf(tumor, variant, this.noCallIsHomRef)
      val oldNormalGt = MakeMixtureVcf.genotypeOf(normal, variant, this.noCallIsHomRef)
      val oldGts      = Seq(oldTumorGt, oldNormalGt)

      val alleles = (oldTumorGt.calls ++ oldNormalGt.calls :+ variant.alleles.ref).filterNot(_ == NoCallAllele).toSet
      val isMultiAllelic = alleles.size > 2
      val isMonomorphic = alleles.size < 2

      if (!isMonomorphic) {
        val (newTumorGt, newNormalGt) = if (isMultiAllelic) {
          (MakeMixtureVcf.noCall(TumorSampleName, variant.alleles), if (tumorOnly) None else Some(MakeMixtureVcf.noCall(NormalSampleName, variant.alleles)))
        }
        else {
          // The locus is bi-allelic in the two samples, so this is the one alt allele that they carry (which need not be
          // the variant's first alt allele)
          val refAllele = variant.alleles.ref
          val altAllele = alleles.find(_ != refAllele).getOrElse(unreachable("A bi-allelic locus must have an alt allele."))

          // Build the Tumor genotype
          val tumorGt = if (oldGts.forall(gt => !gt.isNoCall)) {
            def count(allele: Allele): Int = oldGts.map(_.calls.count(_ == allele)).sum
            val calls = Seq(refAllele, altAllele).filter(a => count(a) > 0).toIndexedSeq match {
              case IndexedSeq(allele) => IndexedSeq(allele, allele)
              case present            => present
            }

            // We count alts only, because ref and no-call can be treated as ref here
            val af = (this.tumorFraction * oldTumorGt.calls.count(_ == altAllele) / 2) + ((1 - this.tumorFraction) * oldNormalGt.calls.count(_ == altAllele) / 2)
            val afs = variant.alleles.alts.map(a => if (a == altAllele) af else 0d)
            Genotype(variant.alleles, TumorSampleName, calls, attrs=Map(MakeTwoSampleMixtureVcf.AlleleFractionField -> afs))
          }
          else {
            MakeMixtureVcf.noCall(TumorSampleName, variant.alleles)
          }

          (tumorGt, if (tumorOnly) None else Some(Genotype(variant.alleles, NormalSampleName, oldNormalGt.calls)))
        }

        val gts = Seq(Some(newTumorGt), newNormalGt).flatten

        val filters = ListBuffer[String]()
        if (isMultiAllelic) filters += MultiAllelicFilter
        if (!tumorOnly && !oldNormalGt.isHomRef) filters += GermlineFilter
        if (newTumorGt.isNoCall) filters += UnknownGtFilter

        out.write(Variant(
          chrom     = variant.chrom,
          pos       = variant.pos,
          id        = variant.id,
          alleles   = variant.alleles,
          filters   = if (filters.nonEmpty) filters.toSet else Variant.PassingFilters,
          attrs     = variant.attrs.filter { case (key, _) => key == "END" },
          genotypes = gts.map(gt => gt.sample -> gt).toMap
        ))
      }
    }

    out.close()
    in.safelyClose()
  }

  /** Generates an iterator over the whole file, or over the intervals if provided. */
  def iterator(in: VcfSource, intervals: Option[PathToIntervals]): Iterator[Variant] = {
    intervals match {
      case None       => in.iterator
      case Some(path) => ByIntervalListVariantIterator(in, IntervalList.fromPath(path).uniqued(false))
    }
  }
}
