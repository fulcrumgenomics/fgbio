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

package com.fulcrumgenomics.vcf.filtration

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam._
import com.fulcrumgenomics.bam.api.{QueryType, SamSource}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.api.{VcfHeader, VcfSource, VcfWriter}

import scala.collection.immutable.ListMap

@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Applies one or more filters to a VCF of somatic variants. The VCF must contain genotype information
    |for the tumor sample.  If the VCF also contains genotypes for one or more other samples, the
    |`--sample` option must be provided to specify the sample whose genotypes to examine and whose reads
    |are present in the BAM file.
    |
    |Various options are available for filtering the reads coming from the BAM file, including
    |`--min-mapping-quality`, `--min-base-quality` and `--paired-reads-only`.  The latter filters to only
    |paired end reads where both reads are mapped.  Reads marked as duplicates, secondary alignments
    |and supplemental alignments are all filtered out.
    |
    |Each available filter may generate annotations in the `INFO` field of the output VCF and optionally,
    |if a threshold is specified, may apply one or more `FILTER`s to applicable variants.
    |
    |In previous versions of this tool, the only available filter was specific to A-base addition artifacts and was
    |referred to as the 'End Repair Artifact Filter.' This filter has been renamed to 'A-tailing Artifact Filter', but
    |its functionality is unchanged. The filter's associated command-line parameters, `INFO` field key, and `FILTER` tag
    |have also been renamed accordingly, as described below.
    |
    |## Available Filters
    |
    |### A-tailing Artifact Filter (previously 'End Repair Artifact Filter')
    |
    |The A-tailing artifact filter attempts to measure the probability that a single-nucleotide mismatch is the product of
    |errors in the template generated during the A-base addition steps that are common to many Illumina library
    |preparation protocols.  The artifacts occur if/when a recessed 3' end is incorrectly filled in with one or more
    |adenines during A-base addition. Incorrect adenine incorporation presents specifically as errors to T at the
    |beginning of reads (and in very short templates, as matching errors to A at the ends of reads).
    |
    |The filter adds the `INFO` field `ATAP` (previously `ERAP`) to SNVs with an A or T alternate allele. This field
    |records the p-value representing the probability of the null hypothesis that the variant is a true mutation, so
    |lower p-values indicate that the variant is more likely an A-tailing artifact. If a threshold p-value is specified,
    |the `FILTER` tag `ATailingArtifact` (previously `EndRepairArtifact`) will be applied to variants with p-values less
    |than or equal to the threshold.
    |
    |Two options are available:
    |
    |* `--a-tailing-distance`  (previously `--end-repair-distance`) allows control over how close to the ends of
    |                          reads/templates errors can be considered to be candidates for the A-tailing artifact.
    |                          Higher values decrease the power of the test, so this should be set as low as possible
    |                          given observed errors.
    |* `--a-tailing-p-value`   (previously `--end-repair-p-value`) the p-value at or below which a filter should be applied.
    |                          If no value is supplied only the `INFO` annotation is produced and no `FILTER` is applied.
    |
    |### End Repair Fill-in Artifact Filter
    |
    |The end repair fill-in artifact filter attempts to measure the probability that a single-nucleotide mismatch is the
    |product of an error in the template generated during the end repair fill-in step that is common to many Illumina
    |library preparation protocols, in which single-stranded 3' overhangs are filled in to create a blunt end. These
    |artifacts originate from single-stranded templates containing damaged bases, often as a consequence of oxidative
    |damage. These DNA lesions, for example 8-oxoguanine, undergo mismatched pairing, which after PCR appear as mutations
    |at the ends of reads.
    |
    |The filter adds the `INFO` field `ERFAP` to records SNVs. This field records the p-value representing the probability of the
    |null hypothesis (e.g. that the variant is a true mutation), so lower p-values indicate that the variant is more likely
    |an end repair fill-in artifact. If a threshold p-value is specified, then the `FILTER` tag `EndRepairFillInArtifact`
    |will be applied to variants with p-values less than or equal to the threshold.
    |
    |Two options are available:
    |
    |* `--end-repair-fill-in-distance`  allows control over how close to the ends of reads/templates errors can be
    |                                   considered to be candidates for the artifact. Higher values decrease the
    |                                   power of the test, so this should be set as low as possible given observed
    |                                   errors.
    |* `--end-repair-fill-in-p-value`   the p-value below which a filter should be applied. If no value is supplied
    |                                   only the annotation is produced and no filtering is performed.
  """)
class FilterSomaticVcf
( @arg(flag='i', doc="Input VCF of somatic variant calls.")         val input: PathToVcf,
  @arg(flag='o', doc="Output VCF of filtered somatic variants.")    val output: PathToVcf,
  @arg(flag='b', doc="BAM file for the tumor sample.")              val bam: PathToBam,
  @arg(flag='s', doc="Sample name in VCF if `> 1` sample present.") val sample: Option[String] = None,
  @arg(flag='m', doc="Minimum mapping quality for reads.")          val minMappingQuality: Int = 30,
  @arg(flag='q', doc="Minimum base quality.")                       val minBaseQuality: Int    = 20,
  @arg(flag='p', doc="Use only paired reads mapped in pairs.")      val pairedReadsOnly: Boolean = false,
  // Developer Note: filter-specific attributes should NOT be given flags as we will likely run
  //                 out of flags that way and end up with a mix of +flag and -flag options.
  @arg(doc="Distance from end of read to implicate A-base addition artifacts. Set to :none: to deactivate the filter.")    val aTailingDistance: Option[Int] = Some(2),
  @arg(doc="Minimum acceptable p-value for the A-base addition artifact test.")                                            val aTailingPValue: Option[Double] = None,
  @arg(doc="Distance from end of read to implicate end repair fill-in artifacts. Set to :none: to deactivate the filter.") val endRepairFillInDistance: Option[Int] = Some(15),
  @arg(doc="Minimum acceptable p-value for the end repair fill-in artifact test.")                                         val endRepairFillInPValue: Option[Double] = None
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(bam)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val vcfIn = VcfSource(input)
    val bamIn = SamSource(bam)
    val vcfSample = sample match {
      case Some(s) =>
        validate(vcfIn.header.samples.contains(s), s"Sample $s not present in VCF.")
        s
      case None =>
        if (vcfIn.header.samples.size == 1) vcfIn.header.samples.head
        else invalid("Must supply --sample when VCF contains more than one sample.")
    }

    val filters = Seq.empty[ReadEndSomaticVariantFilter] ++
      aTailingDistance.map { dist => new ATailingArtifactLikelihoodFilter(dist, aTailingPValue) } ++
      endRepairFillInDistance.map { dist => new EndRepairFillInArtifactLikelihoodFilter(dist, endRepairFillInPValue) }

    if (filters.isEmpty) logger.warning("All filters are disabled and this tool will no-op!")

    val builder = new PileupBuilder(bamIn.dict, mappedPairsOnly=pairedReadsOnly, minBaseQ=minBaseQuality, minMapQ=minMappingQuality)
    val out     = makeWriter(output, vcfIn.header, filters)

    vcfIn.foreach { variant =>
      val gt                = variant.genotypes(vcfSample)
      val applicableFilters = filters.filter(_.appliesTo(gt))

      val variantOut = if (applicableFilters.isEmpty) variant else {
        val iterator = bamIn.query(variant.chrom, variant.pos, variant.pos, QueryType.Overlapping)
        val pileup   = filterPileup(builder.build(iterator, variant.chrom, variant.pos).withoutOverlaps)
        iterator.close()

        val allAnnotations = variant.attrs.toBuffer
        val allFilters     = variant.filters.toBuffer

        applicableFilters.foreach { f =>
          // Generate the set of annotations, push them into the ctxBuilder, then also add any filters to the ctxBuilder
          val annotations = f.annotations(pileup, gt)
          allAnnotations ++= annotations
          allFilters     ++= f.filters(annotations)
        }

        variant.copy(attrs = ListMap(allAnnotations.toSeq:_*), filters = allFilters.toSet)
      }

      out += variantOut
    }

    out.close()
    vcfIn.safelyClose()
    bamIn.safelyClose()
  }

  /**
    * Filters the pileup to remove entries that are from PE reads who's insert size indicates that the
    * aligned base we're looking at is _outside_ the insert!
    */
  private def filterPileup[A <: PileupEntry](pileup: Pileup[A]): Pileup[A] = {
    val filtered = pileup.pile.filter { e =>
      if (e.rec.isFrPair) {
        val (start, end) = Bams.insertCoordinates(e.rec)
        if (pileup.pos < start || pileup.pos > end) {
          logger.warning(
            s"""
              |Ignoring read ${e.rec.name} mapped over variant site ${pileup.refName}:${pileup.pos} that has
              |incompatible insert size yielding insert coordinates of ${pileup.refName}:$start-$end which do not overlap the variant.
            """.stripMargin.trim.linesIterator.mkString(" "))
          false
        }
        else {
          true
        }
      }
      else {
        true
      }
    }

    pileup.copy(pile=filtered)
  }

  /** Builds a VCF writer that had all the necessary headers. */
  private def makeWriter(output: PathToVcf, in: VcfHeader, filters: Seq[SomaticVariantFilter]): VcfWriter = {
    val header = in.copy(
      filters = (in.filters ++ filters.flatMap(_.VcfFilterLines)).distinct,
      infos   = (in.infos   ++ filters.flatMap(_.VcfInfoLines)).distinct
    )

    VcfWriter(output, header)
  }
}
