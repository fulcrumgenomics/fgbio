/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf

import java.io.{Closeable, InputStream}
import java.util
import java.util.NoSuchElementException

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Io
import dagr.commons.CommonsDef._
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriter, VariantContextWriterBuilder}
import htsjdk.variant.vcf._

import scala.collection.JavaConversions._
import scala.collection.JavaConverters._
import scala.io.Source

@clp(
  description =
    """
      |Converts the output of HapCut1/HapCut2 to a VCF.
      |
      |The output of HAPCUT does not include all input variants, but simply those variants that are in phased blocks.
      |This tool takes the original VCF and the output of HapCut, and produces a VCF containing both the variants that
      |were phased and the variants that were not phased, such that all variants in the original file are in the output.
      |
      |The original VCF provided to HAPCUT is assumed to contain a single sample, as HAPCUT only supports a single
      |sample.
      |
      |By default, all phased genotypes are annotated with the "PS" (phase set) FORMAT tag, which by convention is the
      |position of the first variant in the phase set (see the VCF specification).  Furthermore, this tool formats the
      |alleles of a phased genotype using the '|' separator instead of the '/' separator, where the latter indicates the
      |genotype is unphased.  If the option to output phased variants in GATK's ReadBackedPhasing format is used, then
      |the first variant in a phase set will have '/' instead of '|' separating its alleles.  Also, all phased variants
      |will have "PASS" set in its FILTER column, while unphased variants will have "NotPhased" set in their FILTER
      |column.  Unlike GATK's ReadBackedPhasing, homozygous variants will always be unphased.
      |
      |More information about the purpose and operation of GATK's Read-backed phasing, including its output format, can
      |be found here:
      |  http://gatkforums.broadinstitute.org/gatk/discussion/45/purpose-and-operation-of-read-backed-phasing
      |
      |Additional FORMAT fields for phased variants are provided corresponding to per-genotype information produced by
      |HAPCUT1:
      |  1. The "RC" tag gives the counts of calls supporting allele0 and allele1 respectively.
      |  2. The "LC" tag gives the change in likelihood if this SNP is made homozygous or removed.
      |  3. The "MCL" tag gives the maximum change in likelihood if this SNP is made homozygous or removed.
      |  4. The "RMEC" tag gives the reduction in MEC score if we remove this variant altogether.
      |
      |Additional FORMAT fields for phased variants are provided corresponding to per-genotype information produced by
      |HAPCUT2:
      |  1. The "PR" tag is 1 if HapCut2 pruned this variant, 0 otherwise.
      |  2. The "SE" tag gives the confidence (log10) that there is not a switch error occurring immediately before the SNV
      |  3. The "NE" tag gives the confidence (log10) that the SNV is not a mismatch (single SNV) error.
      |HapCut2 should not be run with `--call_homozygous` as the genotypes may be different than the input and is not
      |currently supported.
      |
      |For more information about HAPCUT1, see the source code or paper below.
      |  source code: https://github.com/vibansal/hapcut
      |  HapCut1 paper: An efficient and accurate algorithm for the haplotype assembly problem Bioinformatics. 2008 Aug
      |    15;24(16):i153-9.
      |
      |For more information about HAPCUT2, see the source code below.
      |   source code: https://github.com/pjedge/hapcut2
    """,
  group=ClpGroups.VcfOrBcf
)
class HapCutToVcf
( @arg(flag="v", doc="The original VCF provided to HAPCUT1/HAPCUT2.") val vcf: PathToVcf,
  @arg(flag="i", doc="The output file from HAPCUT1/HAPCUT2.") val input: FilePath,
  @arg(flag="o", doc="The output VCF with both phased and unphased variants.") val output: PathToVcf,
  @arg(flag="r", doc="Output phased variants in GATK's ReadBackedPhasing format.") val gatkPhasingFormat: Boolean = false
) extends FgBioTool with LazyLogging {
  import HapCutType._

  Io.assertReadable(vcf)
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val vcfReader = new VCFFileReader(vcf.toFile, false)
    val iterator  = new HapCutAndVcfMergingIterator(input, vcfReader, gatkPhasingFormat)
    val vcfWriter = makeWriter(output, vcfReader, iterator.hapCutType)

    iterator.foreach(vcfWriter.add)

    vcfReader.safelyClose()
    iterator.safelyClose()
    vcfWriter.close()
  }

  /** Creates a VCF writer, adding extra header lines if the output is the phased VCF. */
  private def makeWriter(path: PathToVcf, vcfReader: VCFFileReader, hapCutType: HapCutType): VariantContextWriter = {
    val inputHeader = vcfReader.getFileHeader
    val builder = new VariantContextWriterBuilder()
      .setOutputFile(path.toFile)
      .setReferenceDictionary(inputHeader.getSequenceDictionary)
      .setOption(Options.INDEX_ON_THE_FLY)
      .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
    val writer: VariantContextWriter = builder.build

    // get the header lines in the input header that we wish to skip/replace with our own definitions
    val headerLinesToSkip = HeaderLines.formatHeaderKeys(hapCutType).flatMap(key => Option(inputHeader.getFormatHeaderLine(key)))
    val headerLines: util.Set[VCFHeaderLine] = new util.HashSet[VCFHeaderLine](
      inputHeader.getMetaDataInSortedOrder.filterNot(headerLinesToSkip.contains).toSet
    )

    // add standard header lines
    VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS)

    // add the new format header lines
    headerLines.addAll(HeaderLines.formatHeaderLines(hapCutType))

    // add the new filter header line if we are to set a filter (not phased) on the unphased variants
    if (gatkPhasingFormat) headerLines.add(HapCut1VcfHeaderLines.NotPhasedFilterHeaderLine)

    // create it
    val outHeader = new VCFHeader(headerLines, inputHeader.getSampleNamesInOrder) // create the header

    // write it and return
    writer.writeHeader(outHeader)
    writer
  }
}

object HeaderLines {
  import HapCutType._
  def formatHeaderKeys(hapCutType: HapCutType): Seq[String] = {
    hapCutType match {
      case HapCut1 => HapCut1VcfHeaderLines.formatHeaderKeys
      case HapCut2 => HapCut2VcfHeaderLines.formatHeaderKeys
      case _       => Seq.empty // empty file
    }
  }

  def formatHeaderLines(hapCutType: HapCutType): Seq[VCFHeaderLine] = {
    hapCutType match {
      case HapCut1 => HapCut1VcfHeaderLines.formatHeaderLines
      case HapCut2 => HapCut2VcfHeaderLines.formatHeaderLines
      case _       => Seq.empty // empty file
    }
  }
}

trait HeaderLines {
  val PhaseSetFormatTag = "PS"
  val PhaseSetFormatDescription = "Phase set for the genotype (position of the first variant)"
  val PhaseSetFormatHeaderLine = new VCFFormatHeaderLine(PhaseSetFormatTag, 1, VCFHeaderLineType.Integer, PhaseSetFormatDescription)

  val NotPhasedFilterName        = "NotPhased"
  val NotPhasedFilterDescription = "The variant was not phased by HapCut."
  val NotPhasedFilterHeaderLine  = new VCFFilterHeaderLine(NotPhasedFilterName, NotPhasedFilterDescription)
}

object HapCut1VcfHeaderLines extends HeaderLines {
  val ReadCountFormatTag         = "RC"
  val ReadCountFormatDescription = "Counts of calls supporting allele0 and allele1 respectively"
  val ReadCountFormatHeaderLine  = new VCFFormatHeaderLine(ReadCountFormatTag, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, ReadCountFormatDescription)

  val LikelihoodChangeFormatTag         = "LC"
  val LikelihoodChangeFormatDescription = "Change in likelihood if this SNP is made homozygous or removed"
  val LikelihoodChangeFormatHeaderLine  = new VCFFormatHeaderLine(LikelihoodChangeFormatTag, 3, VCFHeaderLineType.Float, LikelihoodChangeFormatDescription)

  val MaxLikelihoodChangeFormatTag         = "MLC"
  val MaxLikelihoodChangeFormatDescription = "Maximum change in likelihood if this SNP is made homozygous or removed"
  val MaxLikelihoodChangeFormatHeaderLine  = new VCFFormatHeaderLine(MaxLikelihoodChangeFormatTag, 1, VCFHeaderLineType.Float, MaxLikelihoodChangeFormatDescription)

  val MecReductionFormatTag = "RMEC"
  val MecReductionFormatDescription = "Reduction in MEC score if we remove this variant altogether"
  val MecReductionFormatHeaderLine  = new VCFFormatHeaderLine(MecReductionFormatTag, 1, VCFHeaderLineType.Float, MecReductionFormatDescription)

  /** The VCF header FORMAT keys that will be added for HapCut1 specific-genotype information. */
  val formatHeaderKeys: Seq[String] = {
    Seq(ReadCountFormatTag, LikelihoodChangeFormatTag, MaxLikelihoodChangeFormatTag, MecReductionFormatTag, PhaseSetFormatTag)
  }

  /** The VCF header FORMAT lines that will be added for HapCut1 specific-genotype information. */
  val formatHeaderLines: Seq[VCFHeaderLine] = {
    Seq(ReadCountFormatHeaderLine, LikelihoodChangeFormatHeaderLine, MaxLikelihoodChangeFormatHeaderLine, MecReductionFormatHeaderLine, PhaseSetFormatHeaderLine)
  }
}

object HapCut2VcfHeaderLines extends HeaderLines {
  val PrunedFormatTag         = "PR"
  val PrunedFormatDescription = "1 if HapCut2 pruned this variant, 0 otherwise"
  val PrunedFormatHeaderLine  = new VCFFormatHeaderLine(PrunedFormatTag, 1, VCFHeaderLineType.Integer, PrunedFormatDescription)

  val Log10SwitchErrorFormatTag         = "SE"
  val Log10SwitchErrorFormatDescription = "The confidence (log10) that there is not a switch error occurring immediately before the SNV"
  val Log10SwitchErrorFormatHeaderLine  = new VCFFormatHeaderLine(Log10SwitchErrorFormatTag, 1, VCFHeaderLineType.Float, Log10SwitchErrorFormatDescription)

  val Log10NoErrorFormatTag         = "NE"
  val Log10NoErrorFormatDescription = "The confidence (log10) that the SNV is not a mismatch (single SNV) error."
  val Log10NoErrorFormatHeaderLine  = new VCFFormatHeaderLine(Log10NoErrorFormatTag, 1, VCFHeaderLineType.Float, Log10NoErrorFormatDescription)

  /** The VCF header FORMAT keys that will be added for HapCut2 specific-genotype information. */
  val formatHeaderKeys: Seq[String] = {
    Seq(Log10SwitchErrorFormatTag, Log10NoErrorFormatTag, PhaseSetFormatTag)
  }

  /** The VCF header FORMAT lines that will be added for HapCut2 specific-genotype information. */
  val formatHeaderLines: Seq[VCFHeaderLine] = {
    Seq(Log10SwitchErrorFormatHeaderLine, Log10NoErrorFormatHeaderLine, PhaseSetFormatHeaderLine)
  }
}


/** An iterator over all variants in the provided VCFFileReader.
  *
  * @param hapCutPath the output file with phased variants produced by HapCut.
  * @param vcfReader the VCF file with the original variants provided to HapCut.
  * @param gatkPhasingFormat true to output in GATK's ReadBackedPhasing format, false if to use the recommendations in the VCF spec.
  */
private class HapCutAndVcfMergingIterator(hapCutPath: FilePath,
                                          vcfReader: VCFFileReader,
                                          gatkPhasingFormat: Boolean)
  extends Iterator[VariantContext] with Closeable {
  import HapCutType.HapCutType

  private val sourceIterator = vcfReader.toStream.zipWithIndex.iterator.buffered
  private val hapCutReader   = HapCutReader(path=hapCutPath)
  private val sampleName     = vcfReader.getFileHeader.getSampleNamesInOrder.head

  def hasNext(): Boolean = {
    if (sourceIterator.isEmpty && hapCutReader.hasNext()) throw new IllegalStateException("HapCut has more phased variants but no more variants in the input")
    sourceIterator.hasNext
  }

  /** Returns either the phased variant with added HapCut-specific genotype information (Left), or the original variant (Right). */
  def next(): VariantContext = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false")
    val (sourceContext, sourceOffset) = sourceIterator.next()
    if (!hapCutReader.hasNext()) formatSourceContext(sourceContext)
    else {
      val HapCutOffsetAndCall(offset, callOption) = hapCutReader.next()
      if (offset != sourceOffset+1) throw new IllegalStateException("BUG: calls are out of order")
      callOption match {
        case None =>
          formatSourceContext(sourceContext)
        case Some(hapCutCall) =>
          require(hapCutCall.pos == sourceContext.getStart)
          require(offset == hapCutCall.offset)
          val hapCutContext = hapCutCall.toVariantContext(sampleName)
          replaceGenotypes(source=sourceContext, genotype=hapCutContext.getGenotype(0), isPhased = !gatkPhasingFormat || !hapCutCall.firstInBlock)
      }
    }
  }

  /** Returns the the HapCut type (HapCut1 or HapCut2), or Unknown if the file was empty. */
  def hapCutType: HapCutType = hapCutReader.hapCutType

  def close(): Unit = { this.hapCutReader.close }

  /** Returns a new variant context with the phase unset (if it was set) and appropriately formatted (filtered as
    * not phased if `gatkPhasingFormat` is true).
    */
  private def formatSourceContext(sourceContext: VariantContext): VariantContext = {
    val hasPhasedGenotype = sourceContext.getGenotypes.exists(_.isPhased)
    val hasPhasingSetId   = sourceContext.getGenotypes.exists(_.hasExtendedAttribute(HapCut1VcfHeaderLines.PhaseSetFormatTag))

    if (!hasPhasedGenotype && !hasPhasingSetId && !gatkPhasingFormat) sourceContext
    else {
      val builder = new VariantContextBuilder(sourceContext)
      if (hasPhasedGenotype || hasPhasingSetId) {
        // unset the phase and remove the phasing set ID if the input has phase set
        builder.genotypes(sourceContext.getGenotypes().map { g =>
          val builder = new GenotypeBuilder(g).phased(false)
          val attrs = g.getExtendedAttributes.filterNot { case (tag, value) =>  tag == HapCut1VcfHeaderLines.PhaseSetFormatTag }
          builder.noAttributes()
          builder.attributes(attrs)
          builder.make()
        })
      }
      if (gatkPhasingFormat) {
        // set the variant as filtered due to not being phased
        builder.filter(HapCut1VcfHeaderLines.NotPhasedFilterName)
      }
      builder.make()
    }
  }

  /** Replaces the original genotype with HapCut's genotype and adds any HapCut-specific genotype information. */
  private def replaceGenotypes(source: VariantContext, genotype: Genotype, isPhased: Boolean): VariantContext = {
    val builder = new VariantContextBuilder(source)
    val sourceAlleles = source.getAlleles.toSeq
    val genotypeAlleles = genotype.getAlleles.map {
      allele => sourceAlleles.find(a => a.toString == allele.toString) match {
        case None =>
          throw new IllegalStateException(s"Could not find allele '$allele' in source alleles: " + sourceAlleles.map{_.toString}.mkString(", "))
        case Some(a) => a
      }
    }
    val sourceGenotype = source.getGenotype(0)
    val genotypeBuilder = new GenotypeBuilder(sourceGenotype).alleles(genotypeAlleles).phased(isPhased)
    genotypeBuilder.attributes(genotype.getExtendedAttributes)
    builder.genotypes(genotypeBuilder.make()).make()
  }
}

sealed trait BlockInfo {
  /** The 1-based variant offset from the first variant given to HapCut. */
  def offset: Int
  /** The total number of variants that this block spans, including unphased variants. */
  def len: Int
  /** The total number of variants phased. */
  def phased: Int
  /** The number of bases this block spans. */
  def span: Int
  /** The number of fragments contributing to variants in this block. */
  def fragments: Int
  /** The 1-based variant offset (from the first variant given to HapCut) of the last genotype in this block. */
  def offsetOfLastVariant: Int = this.offset + this.len - 1
  /** The phase set identifier for this block (typically the position of the first variant). Throws an exception if
    * the `phaseSetOption` was None when this was constructed. */
  def phaseSet: Int = phaseSetOption.get
  /** Updates the phase set to which this block belongs. */
  def withPhaseSet(phaseSet: Option[Int]): BlockInfo
  /** Stores the phase set, or None if it has not been set. */
  protected def phaseSetOption: Option[Int]
}

/** Haplotype block information produced by HapCut.
  *
  * @param phaseSetOption the phase set identifer, or None if not yet available.
  * @param offset the 1-based variant offset from the first variant given to HapCut.
  * @param len the total number of variants that this block spans, including unphased variants.
  * @param phased the total number of variants phased.
  * @param span the number of bases this block spans.
  * @param mecScore the minimum error correction (MEC) score.  For HAPCUT2, this will be -1.
  * @param fragments the number of fragments contributing to variants in this block.
  */
private case class HapCut1BlockInfo (protected val phaseSetOption: Option[Int] = None,
                                     offset: Int,
                                     len: Int,
                                     phased: Int,
                                     span: Int,
                                     mecScore: Float,
                                     fragments: Int) extends BlockInfo {
  def withPhaseSet(phaseSet: Option[Int]): HapCut1BlockInfo = this.copy(phaseSetOption=phaseSet)
}

/** Haplotype block information produced by HapCut2.
  *
  * @param phaseSetOption the phase set identifer, or None if not yet available.
  * @param offset the 1-based variant offset from the first variant given to HapCut.
  * @param len the total number of variants that this block spans, including unphased variants.
  * @param phased the total number of variants phased.
  * @param span the number of bases this block spans.
  * @param fragments the number of fragments contributing to variants in this block.
  */
private case class HapCut2BlockInfo (protected val phaseSetOption: Option[Int] = None,
                                     offset: Int,
                                     len: Int,
                                     phased: Int,
                                     span: Int,
                                     fragments: Int) extends BlockInfo {
  def withPhaseSet(phaseSet: Option[Int]): HapCut2BlockInfo = this.copy(phaseSetOption=phaseSet)
}

private object GenotypeInfo {
  import HapCutType._
  def apply(info: String, hapCutType: HapCutType, thresholdPruning: Boolean = false): GenotypeInfo = {
    hapCutType match {
      case HapCut1 => HapCut1GenotypeInfo(info=info)
      case HapCut2 => HapCut2GenotypeInfo(info=info, thresholdPruning=thresholdPruning)
      case _       => unreachable("Unknown HapCut type when building a GenotypeInfo.")
    }
  }
}

/** Genotype-level information specific to HapCut1 or HapCut2. */
private trait GenotypeInfo {
  def addTo(builder: GenotypeBuilder): GenotypeBuilder
}

/** Genotype-level information produced by HapCut1 */
private class HapCut1GenotypeInfo private(val readCounts: List[Int],
                                          val likelihoods: List[Float],
                                          val delta: Float,
                                          val rMEC: Float) extends GenotypeInfo {
  import HapCut1VcfHeaderLines._
  /** Adds the HapCut-specific genotype information to the given genotype builder. */
  def addTo(builder: GenotypeBuilder): GenotypeBuilder = {
    builder
      .attribute(ReadCountFormatTag, this.readCounts.asJava)
      .attribute(LikelihoodChangeFormatTag, this.likelihoods.asJava)
      .attribute(MaxLikelihoodChangeFormatTag, this.delta)
      .attribute(MecReductionFormatTag, this.rMEC)
  }
}

/** Values and methods for HapCut-specific genotype information. */
private object HapCut1GenotypeInfo {
  // A developer's notes
  //  R0,R1:L00,L01,L11:delta:rMEC
  //  "%d,%d:%0.1f,%0.1f,%0.1f:%0.1f:%0.1f"
  //  R0,R1: counts of calls supporting allele0 and allele1 respectively
  //  L00,L01,L11: change in likelihood if this SNP is made homozygous or removed
  //  delta: maximum change in likelihood if this SNP is made homozygous or removed
  //  rMEC: reduction in MEC score if we remove this variant altogether
  val HapCut1CallInfoPattern = """(\d+),(\d+):(\d+\.\d+),(\d+\.\d+),(\d+\.\d+):(\d+\.\d+):(\d+\.\d+)""".r

  /** Parses the genotype information produced by HapCut. */
  def apply(info: String): GenotypeInfo = {
    val tokens = info.split(":")
    new HapCut1GenotypeInfo(
      readCounts  = tokens(0).split(",").map(_.toInt).toList,
      likelihoods = tokens(1).split(",").map(_.toFloat).toList,
      delta       = tokens(2).toFloat,
      rMEC        = tokens(3).toFloat
    )
  }
}

/** Genotype-level information produced by HapCut2 */
private class HapCut2GenotypeInfo private(val pruned: Boolean, val log10SwitchError: Double, val log10NoError: Double) extends GenotypeInfo {
  import HapCut2VcfHeaderLines._
  def addTo(builder: GenotypeBuilder): GenotypeBuilder = {
    builder
      .attribute(PrunedFormatTag, if (this.pruned) 1 else 0)
      .attribute(Log10SwitchErrorFormatTag, this.log10SwitchError)
      .attribute(Log10NoErrorFormatTag, this.log10NoError)
  }
}

/** Values and methods for HapCut2-specific genotype information. */
private object HapCut2GenotypeInfo {
  /** Parses the genotype information produced by HapCut2. */
  def apply(info: String, thresholdPruning: Boolean): GenotypeInfo = {
    val tokens = info.split("\t")
    new HapCut2GenotypeInfo(
      pruned           = "1" == tokens(0) || thresholdPruning,
      log10SwitchError = tokens(1).toDouble,
      log10NoError     = tokens(2).toDouble
    )
  }
}

/** Information about a variant phased by HapCut.
  *
  * @param block the haplotype block to which this variant is associated.
  * @param offset the 1-based offset of the variant relative to the variants given to HapCut.
  * @param hap1Allele the first allele.
  * @param hap2Allele the second allele.
  * @param contig the contig.
  * @param pos the position.
  * @param ref the reference allele.
  * @param alts the alternate alleles.
  * @param genotype the genotype.
  * @param info the HapCut-specific genotype information.
  * @param phaseSet the phase set identifier for the phase block to which this call belongs
  */
private case class HapCutCall private(block: BlockInfo,
                                      offset: Int,
                                      hap1Allele: Int, hap2Allele: Int,
                                      contig: String, pos: Int,
                                      ref: String, alts: String,
                                      genotype: String, info: GenotypeInfo,
                                      phaseSet: Int
                                     ) {
  import HapCut1VcfHeaderLines._

  /** true if this variant is the first phased variant in a phased block, false otherwise. */
  def firstInBlock: Boolean = phaseSet == pos

  /** Converts the HapCut variant representation to a VariantContext. */
  def toVariantContext(sampleName: String): VariantContext = {
    // Parse the alleles
    val refAllele = Allele.create(this.ref, true)
    val altAlleles = this.alts.split(",").map { alt =>
      Allele.create(alt, false)
    }.toSeq
    val alleles = refAllele +: altAlleles
    val allelesCollection = alleles.asJavaCollection

    // Assumes the genotype is before the first ":" token
    val genotypeAlleles = this.genotype.split(":").head.split("[/|]").map(_.toInt).map(i => alleles(i))
    val genotype = this.info.addTo(new GenotypeBuilder(sampleName, genotypeAlleles.toList))
      .attribute(PhaseSetFormatTag, phaseSet)
      .make()

    // build the context, and make sure to recompute the end
    new VariantContextBuilder(
      s"${this.offset}",
      this.contig,
      this.pos.toLong,
      this.pos.toLong,
      allelesCollection
    )
      .computeEndFromAlleles(seqAsJavaList(alleles), this.pos)
      .genotypes(genotype).make()
  }
}

private object HapCutCall {
  import HapCutType._

  /** Parse a variant line.
    *
    * @param callLine the line to parse.
    * @param firstVariantPosition the position of the first variant in the block to which this variant belongs, or None
    *                             if this variant is the first variant in the block.
    **/
  def toHapCutCall(callLine: String, block: BlockInfo, firstVariantPosition: Option[Int], hapCutType: HapCutType): HapCutCall = {
    val tokens = callLine.split('\t')
    if (9 != tokens.length && 11 != tokens.length) {
      throw new IllegalStateException(s"Did not parse 9 or 12 fields (parsed ${tokens.length}) in the HapCut/HapCut2 line: " + callLine.trim)
    }
    val offset = tokens(0).toInt
    val pos    = tokens(4).toInt
    val hap1Allele = if ("-" == tokens(1)) -1 else tokens(1).toInt
    val hap2Allele = if ("-" == tokens(2)) -1 else tokens(2).toInt
    val thresholdPruning = hapCutType == HapCut2 && (hap1Allele < 0 || hap2Allele < 0)
    new HapCutCall(
      block        = block,
      offset       = offset,
      hap1Allele   = hap1Allele,
      hap2Allele   = hap2Allele,
      contig       = tokens(3),
      pos          = pos,
      ref          = tokens(5),
      alts         = tokens(6),
      genotype     = tokens(7),
      info         = GenotypeInfo(tokens.drop(8).mkString("\t"), hapCutType=hapCutType, thresholdPruning=thresholdPruning),
      phaseSet     = firstVariantPosition.getOrElse(pos)
    )
  }
}

object HapCutType extends Enumeration {
  type HapCutType = Value
  val HapCut1, HapCut2, Unknown = Value
}

private[vcf] case class HapCutOffsetAndCall(offset: Int, call: Option[HapCutCall] = None)

object HapCutReader {
  val HapCut1BlockLinePattern = """^BLOCK: offset: (\d+) len: (\d+) phased: (\d+) SPAN: (\d+) MECscore (\d+.\d+) fragments (\d+)$""".r
  val HapCut2BlockLinePattern = """^BLOCK: offset: (\d+) len: (\d+) phased: (\d+) SPAN: (\d+) fragments (\d+)$""".r

  def apply(stream: InputStream): HapCutReader = new HapCutReader(Source.fromInputStream(stream).getLines(), Some(stream))

  def apply(path: FilePath): HapCutReader = apply(Io.toInputStream(path))
}

/** Reads phased blocks and variants from a HapCut output file.  For each variant call given to HapCut, returns either
  * some `HapCutCall` if that variant was phased, or None if it was not phased.  HapCut does not include unphased
  * variants in its output, but does provide the 1-based offset of every phased variant relative to the variant calls
  * given to HapCut.  Therefore, we are able to determine if None should be returned (unphased, not in the HapCut file)
  * or if a `HapCutCall` should be returned (phased, in the HapCut file).  Nonetheless, variants that were not phased
  * after the last HapCut phased variant will not be returned (hasNext will be false), as HapCut doesn't provide the
  * total # of variants that were given to HapCut, so these trailing unphased variant calls must be treated specially
  * for now.
  */
private[vcf] class HapCutReader(iterator: Iterator[String],
                                private[this] val source: Option[{ def close(): Unit }] = None)
  extends Iterator[HapCutOffsetAndCall] with Closeable {
  import HapCutType._
  import scala.collection.mutable.ListBuffer

  private val lineIterator = iterator.buffered
  val hapCutType: HapCutType = {
    if (lineIterator.isEmpty) Unknown // Default to true if the file is empty
    else {
      lineIterator.head match {
        case HapCutReader.HapCut1BlockLinePattern(offset, len, phased, span, mecScore, fragments) => HapCut1
        case HapCutReader.HapCut2BlockLinePattern(offset, len, phased, span, fragments)           => HapCut2
        case blockLine => throw new IllegalStateException("Could not parse block line: " + blockLine)
      }
    }
  }
  private var nextBlockInfo: Option[BlockInfo] = this.forNextBlockInfo()
  private val calls = new CallBuffer


  /** Simple class to store calls within a range of offsets. */
  private class CallBuffer {
    private val callBuffer = new ListBuffer[HapCutOffsetAndCall]()
    private var firstOffset = 1

    /** Add a call to the buffer.  This may add empty calls with less than offset of this call relative to the call
      * with the greatest offset ever present in this buffer. */
    def add(call: HapCutCall): Unit = {
      require(firstOffset <= call.offset)
      val offsetAndCall = new HapCutOffsetAndCall(offset=call.offset, call=Some(call))
      // add empty values *prior* to this call
      while (firstOffset + this.callBuffer.length < call.offset) {
        this.callBuffer.append(new HapCutOffsetAndCall(offset=firstOffset+this.callBuffer.length, call=None))
      }
      // add the call
      if (this.callBuffer.nonEmpty && call.offset <= firstOffset + this.callBuffer.length - 1) {
        this.callBuffer(call.offset-firstOffset) = offsetAndCall
      }
      else {
        this.callBuffer.append(offsetAndCall)
      }
    }

    /** Removes the next call. */
    def removeFirst(): HapCutOffsetAndCall = {
      val call = this.callBuffer.remove(0)
      this.firstOffset = call.offset + 1
      call
    }

    def isEmpty: Boolean = this.callBuffer.isEmpty

    def nonEmpty: Boolean = this.callBuffer.nonEmpty
  }

  def hasNext(): Boolean = {
    if (calls.nonEmpty) true
    else this.advance()
  }

  def next(): HapCutOffsetAndCall = {
    if (!hasNext()) throw new NoSuchElementException("Calling next() when hasNext() is false")
    calls.removeFirst()
  }

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  def close(): Unit = this.source.foreach(_.close())

  /** Reads in the next block info. */
  private def forNextBlockInfo(): Option[BlockInfo] = {
    // find the next block
    while (this.lineIterator.hasNext && !this.lineIterator.head.startsWith("BLOCK:")) {
      this.lineIterator.next()
    }
    if (this.lineIterator.hasNext) {
      // Read the block information
      val blockLine = this.lineIterator.next().trim
      val currentBlock: BlockInfo = blockLine match {
        case HapCutReader.HapCut1BlockLinePattern(offset, len, phased, span, mecScore, fragments) =>
          if (hapCutType != HapCut1) throw new IllegalStateException(s"Found a HapCut1 block line but the file type was $hapCutType")
          new HapCut1BlockInfo(None, offset.toInt, len.toInt, phased.toInt, span.toInt, mecScore.toFloat, fragments.toInt)
        case HapCutReader.HapCut2BlockLinePattern(offset, len, phased, span, fragments) =>
          if (hapCutType != HapCut2) throw new IllegalStateException(s"Found a HapCut2 block line but the file type was $hapCutType")
          new HapCut2BlockInfo(None, offset.toInt, len.toInt, phased.toInt, span.toInt, fragments.toInt)
        case _ => throw new IllegalStateException("Could not parse block line: " + blockLine)
      }
      Some(currentBlock)
    }
    else None
  }

  /** Reads in the next set of variants for the given block. */
  private def forNextBlockCalls(block: BlockInfo): Unit = {
    // Read in the phased variants and add them to hte set of calls
    var firstVariantPosition: Option[Int] = None
    var i = 0
    while (i < block.phased) {
      if (!lineIterator.hasNext) throw new IllegalStateException("Reached the end of file when searching for a call")
      val call = HapCutCall.toHapCutCall(callLine=lineIterator.next(), block=this.nextBlockInfo.get, firstVariantPosition=firstVariantPosition, hapCutType=hapCutType)
      if (firstVariantPosition.isEmpty) firstVariantPosition = Some(call.pos)
      this.calls.add(call=call)
      i += 1
    }
  }

  /** Gets the variants for the next block, and any subsequent block that overlaps this block or others subsequently
    * read in.  Returns true if any block was read, false otherwise. */
  private def advance(): Boolean = {
    // Read in the variants for the next block
    val firstBlock = this.nextBlockInfo match {
      case None => return false
      case Some(block) =>
        forNextBlockCalls(block=block)
        block
    }
    require(this.calls.nonEmpty)

    // read in blocks until the start offset of the next block is greater than the maximum offset seen so far
    var maxOffset = firstBlock.offsetOfLastVariant
    this.nextBlockInfo = forNextBlockInfo()
    while (this.nextBlockInfo.exists(_.offset <= maxOffset)) {
      forNextBlockCalls(block=this.nextBlockInfo.get)
      maxOffset = Math.max(maxOffset, this.nextBlockInfo.get.offsetOfLastVariant)
      this.nextBlockInfo = forNextBlockInfo()
    }

    true
  }
}
