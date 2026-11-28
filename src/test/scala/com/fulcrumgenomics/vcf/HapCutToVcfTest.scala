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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.HapCutType.{HapCut1, HapCut2, HapCutType}
import com.fulcrumgenomics.vcf.api.{Variant, VcfSource}
import org.scalatest.{Assertion, ParallelTestExecution}

/**
  * Tests for HapCutToVcf.
  */
class HapCutToVcfTest extends UnitSpec with ParallelTestExecution {

  private val dir            = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf/testdata")
  private val originalVcf    = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.vcf")
  private val hapCut1Out     = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut")
  private val hapCut1Vcf     = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut.vcf")
  private val hapCut1GatkVcf = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut.gatk.vcf")
  private val hapCut2Out     = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut2")
  private val hapCut2Vcf     = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut2.vcf")
  private val hapCut2GatkVcf = dir.resolve("NA12878.GIABPedigreev0.2.17.41100000.41300000.hapcut2.gatk.vcf")

  // For testing HapCut2 producing phased blocks overlapping other phased blocks.
  private val outOfOrderIn     = dir.resolve("blocks_out_of_order.vcf")
  private val outOfOrderOut    = dir.resolve("blocks_out_of_order.hapcut2")
  private val outOfOrderOutVcf = dir.resolve("blocks_out_of_order.hapcut2.vcf")

  // For testing HapCut2 with missing variants in the input VCF
  private val missingVariantsIn  = dir.resolve("missing_leading_variants.vcf")
  private val missingVariantsOut = dir.resolve("missing_leading_variants.hapcut2")

  // For testing HapCut2 with missing genotype info
  private val missingGenotyeInfoIn = dir.resolve("hapcut2_for_missing_genotype_info.vcf")
  private val noSwitchErrorsIn     = dir.resolve("no_switch_errors.hapcut2")
  private val skipPruneIn          = dir.resolve("skip_prune.hapcut2")
  private val withSwitchErrors     = dir.resolve("with_switch_errors.hapcut2")

  // For testing HapCutToVcf with IUPAC codes
  private val withIupacIn     = dir.resolve("with_iupac.vcf")
  private val withIupacOut    = dir.resolve("with_iupac.hapcut")

  private def countVcfRecords(vcf: PathToVcf): Int = {
    val vcfReader = VcfSource(vcf, allowMissingFieldsInHeader=true)
    yieldAndThen(vcfReader.iterator.length)(vcfReader.close())
  }

  private def compareVcfs(newVcf: PathToVcf, originalVcf: PathToVcf): Unit = {
    val newVcfReader      = VcfSource(newVcf, allowMissingFieldsInHeader=true)
    val originalVcfReader = VcfSource(originalVcf, allowMissingFieldsInHeader=true)

    for (newVariant <- newVcfReader) {
      originalVcfReader.exists { originalVariant =>
        originalVariant.chrom == newVariant.chrom &&
        originalVariant.pos   == newVariant.pos &&
        originalVariant.end   == newVariant.end
      } shouldBe true
    }
  }

  /** Checks that the records (not the header) of the VCF are the same text as those of the expected VCF. */
  private def compareRecords(vcf: PathToVcf, expected: PathToVcf): Assertion = {
    def records(path: PathToVcf): Seq[String] = Io.readLines(path).filterNot(_.startsWith("#")).toSeq
    records(vcf) should contain theSameElementsInOrderAs records(expected)
  }

  private def isPhased(variant: Variant, gatkPhasingFormat: Boolean): Boolean = {
    if (gatkPhasingFormat) !variant.isFiltered // are marked as passed filter
    else variant.gts.exists(_.phased) // are marked as phased
  }

  private def getNumPhasedFromVcf(path: PathToVcf, gatkPhasingFormat: Boolean): Int = {
    val vcfReader = VcfSource(path, allowMissingFieldsInHeader=true)
    val numPhased = vcfReader.iterator.count { variant => isPhased(variant, gatkPhasingFormat) }
    vcfReader.close()
    numPhased
  }

  private def hasPhasingSetFormatTagButUnphased(path: PathToVcf, gatkPhasingFormat: Boolean): Boolean = {
    val vcfReader = VcfSource(path, allowMissingFieldsInHeader=true)
    val hasPhasingSetTag = vcfReader
      .iterator
      .filterNot { variant => isPhased(variant, gatkPhasingFormat) }
      .exists { variant => variant.gts.exists(_.attrs.contains(HapCut1VcfHeaderLines.PhaseSetFormatTag)) }
    vcfReader.close()
    hasPhasingSetTag
  }

  private def checkHapCutReader(path: FilePath, hapCutType: HapCutType) = {
    val reader = HapCutReader(path)
    reader.hapCutType shouldBe hapCutType
    val allCalls = reader.toSeq
    val calls = allCalls.flatMap(_.call)
    allCalls.length shouldBe 342 // 342 total variants
    calls.length shouldBe 237 // 237 phased variants
    calls.map(_.phaseSet).distinct.length shouldBe 8 // eight phased blocks

    // Check the second block (two variants).  The first variants is 1/0 while the second is 0/1.
    {
      val call = calls(3)
      call.phaseSet shouldBe 41106449
      call.hap1Allele shouldBe 1
      call.hap2Allele shouldBe 0
      val gt = call.toGenotype("Sample")
      gt.phased shouldBe true
      gt.calls.map(_.value) should contain theSameElementsInOrderAs Seq("CT", "C")
    }
    {
      val call = calls(4)
      call.phaseSet shouldBe 41106449
      call.hap1Allele shouldBe 0
      call.hap2Allele shouldBe 1
      val gt = call.toGenotype("Sample")
      gt.phased shouldBe true
      gt.calls.map(_.value) should contain theSameElementsInOrderAs Seq("T", "G")
    }
  }

  "HapCutReader" should "read in a HAPCUT1 file" in {
    checkHapCutReader(hapCut1Out, HapCut1)
  }

  it should "read in a HAPCUT2 file" in {
    checkHapCutReader(hapCut2Out, HapCut2)
  }

  it should "read in a HAPCUT1 file that has phased genotypes" in {
    val input = dir.resolve("block_has_phased_genotypes.hapcut")
    val reader = HapCutReader(input)
    val allCalls = reader.toSeq
    val calls = allCalls.flatMap(_.call)
    allCalls.length shouldBe 8 // 8 total variants
    calls.length shouldBe 3 // 3 phased variants
    calls.map(_.phaseSet).distinct.length shouldBe 1 // a single phased block
    reader.close()  }

  "HapCutToVcf" should "convert a HAPCUT1 file to VCF in both GATK and VCF-spec phasing format" in {
    Iterator(true, false).foreach { gatkPhasingFormat =>
      val expectedOutput = if (gatkPhasingFormat) hapCut1GatkVcf else hapCut1Vcf
      val out = makeTempFile("hap_cut_to_vcf.hapcut", ".vcf")

      new HapCutToVcf(
        vcf               = originalVcf,
        input             = hapCut1Out,
        output            = out,
        gatkPhasingFormat = gatkPhasingFormat
      ).execute()

      // check that we have the same # of records in the output as the input
      countVcfRecords(out) shouldBe countVcfRecords(originalVcf)

      // check that all records in the output are found in the input
      compareVcfs(out, originalVcf)

      // get the # of phased variants from the output
      val numPhasedFromOut = getNumPhasedFromVcf(out, gatkPhasingFormat)

      // check that the # of variants phased in the output agrees with the # of phased calls produced by HapCut
      val hapCutReader = HapCutReader(hapCut1Out)
      val numPhasedFromHapCut = hapCutReader.flatMap(_.call).length
      numPhasedFromOut shouldBe numPhasedFromHapCut
      hapCutReader.close()

      // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
      numPhasedFromOut shouldBe getNumPhasedFromVcf(expectedOutput, gatkPhasingFormat)

      // check that if a variant is not phased it does not have a PS tag
      hasPhasingSetFormatTagButUnphased(out, gatkPhasingFormat) shouldBe false

      // check that every record is as expected
      compareRecords(out, expectedOutput)
    }
  }

  it should "convert a HAPCUT2 file to VCF in both GATK and VCF-spec phasing format" in {
    Iterator(true, false).foreach { gatkPhasingFormat =>
      val expectedOutput = if (gatkPhasingFormat) hapCut2GatkVcf else hapCut2Vcf
      val out = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

      new HapCutToVcf(
        vcf               = originalVcf,
        input             = hapCut2Out,
        output            = out,
        gatkPhasingFormat = gatkPhasingFormat
      ).execute()

      val source = VcfSource(out)
      try {
        source.header.format("NE").description shouldBe "The confidence (log10) that the SNV is not a mismatch (single SNV) error."
        source.header.format("SE").description shouldBe "The confidence (log10) that there is not a switch error occurring immediately before the SNV"
      }
      finally source.close()

      // check that we have the same # of records in the output as the input
      countVcfRecords(out) shouldBe countVcfRecords(originalVcf)

      // check that all records in the output are found in the input
      compareVcfs(out, originalVcf)

      // get the # of phased variants from the output
      val numPhasedFromOut = getNumPhasedFromVcf(out, gatkPhasingFormat)

      // check that the # of variants phased in the output agrees with the # of phased calls produced by HapCut
      val hapCutReader = HapCutReader(hapCut2Out)
      val numPhasedFromHapCut = hapCutReader.flatMap(_.call).length
      numPhasedFromOut shouldBe numPhasedFromHapCut
      hapCutReader.close()

      // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
      numPhasedFromOut shouldBe getNumPhasedFromVcf(expectedOutput, gatkPhasingFormat)

      // check that if a variant is not phased it does not have a PS tag
      hasPhasingSetFormatTagButUnphased(out, gatkPhasingFormat) shouldBe false

      // check that every record is as expected
      compareRecords(out, expectedOutput)
    }
  }

  it should "convert an HAPCUT2 file to VCF when there are overlapping phase blocks" in {
    val out       = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

    new HapCutToVcf(
      vcf               = outOfOrderIn,
      input             = outOfOrderOut,
      output            = out,
      gatkPhasingFormat = false
    ).execute()

    // check that we have the same # of records in the output as the input
    countVcfRecords(out) shouldBe countVcfRecords(outOfOrderIn)

    // check that all records in the output are found in the input
    compareVcfs(out, outOfOrderIn)

    // get the # of phased variants from the output
    val numPhasedFromOut = getNumPhasedFromVcf(out, false)

    // check that the # of variants phased in the output agrees with the # of phased calls produced by HapCut
    val hapCutReader = HapCutReader(outOfOrderOut)
    val numPhasedFromHapCut = hapCutReader.flatMap(_.call).length
    numPhasedFromOut shouldBe numPhasedFromHapCut
    hapCutReader.close()

    // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
    numPhasedFromOut shouldBe getNumPhasedFromVcf(outOfOrderOutVcf, false)

    // check that if a variant is not phased it does not have a PS tag
    hasPhasingSetFormatTagButUnphased(out, false) shouldBe false
  }

  it should "convert an empty HAPCUT1/HAPCUT2 file to VCF in both GATK and VCF-spec phasing format" in {
    Iterator(true, false).foreach { gatkPhasingFormat =>
      val out       = makeTempFile("hap_cut_to_vcf.hapcut", ".vcf")
      val hapCutOut = makeTempFile("hap_cut_to_vcf.hapcut", ".hapcut")

      Io.writeLines(hapCutOut, Seq.empty)

      new HapCutToVcf(
        vcf               = originalVcf,
        input             = hapCutOut,
        output            = out,
        gatkPhasingFormat = gatkPhasingFormat
      ).execute()

      // check that we have the same # of records in the output as the input
      countVcfRecords(out) shouldBe countVcfRecords(originalVcf)

      // check that all records in the output are found in the input
      compareVcfs(out, originalVcf)

      // get the # of phased variants from the output
      getNumPhasedFromVcf(out, gatkPhasingFormat) shouldBe 0

      // check that if a variant is not phased it does not have a PS tag
      hasPhasingSetFormatTagButUnphased(out, gatkPhasingFormat) shouldBe false
    }
  }

  Seq(".vcf.gz", ".bcf").foreach { extension =>
    it should s"preserve an SVLEN-only deletion's indexed span with empty HapCut input in $extension" in {
      val input = makeTempFile("hapcut-svlen", ".vcf")
      Io.writeLines(input, Seq(
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=100000>",
        "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Structural variant length\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1").mkString("\t"),
        Seq("chr1", "16380", ".", "A", "<DEL>", ".", "PASS", "SVLEN=-50", "GT", "0/1").mkString("\t")
      ))
      val hapcut = makeTempFile("hapcut-svlen", ".hapcut")
      Io.writeLines(hapcut, Seq.empty)
      val output = makeTempFile("hapcut-svlen", extension)
      new HapCutToVcf(vcf=input, input=hapcut, output=output).execute()
      val source = VcfSource(output)
      try {
        val IndexedSeq(variant) = source.iterator.toIndexedSeq
        variant.end shouldBe 16430
        source.query("chr1", 16400, 16400).map(_.pos).toIndexedSeq shouldBe IndexedSeq(16380)
      }
      finally source.close()
    }
  }

  it should "fail when there are missing variants in the input VCF" in {
    val out       = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

    an[Exception] should be thrownBy new HapCutToVcf(
      vcf               = missingVariantsIn,
      input             = missingVariantsOut,
      output            = out,
      gatkPhasingFormat = false
    ).execute()
  }

  it should "fail when the input VCF has no samples" in {
    val vcf = makeTempFile("hap_cut_to_vcf.sites_only", ".vcf")
    Io.writeLines(vcf, Seq(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=1000>",
      Seq("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO").mkString("\t"),
      Seq("chr1", "10", ".", "A", "C", ".", ".", ".").mkString("\t")
    ))
    val hapCutOut = makeTempFile("hap_cut_to_vcf.sites_only", ".hapcut")
    Io.writeLines(hapCutOut, Seq.empty)
    val out = makeTempFile("hap_cut_to_vcf.sites_only", ".vcf")

    val ex = intercept[FailureException] {
      new HapCutToVcf(vcf=vcf, input=hapCutOut, output=out, gatkPhasingFormat=false).execute()
    }
    ex.message.value should include ("The input VCF has no samples")
  }

  it should "support missing genotype info in HapCut2" in {
    Seq(noSwitchErrorsIn, skipPruneIn, withSwitchErrors).foreach { hapCut2In =>
      val expectedOutput = PathUtil.replaceExtension(hapCut2In, ".vcf")
      val out = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

      new HapCutToVcf(
        vcf               = missingGenotyeInfoIn,
        input             = hapCut2In,
        output            = out,
        gatkPhasingFormat = false
      ).execute()

      // check that we have the same # of records in the output as the input
      countVcfRecords(out) shouldBe countVcfRecords(missingGenotyeInfoIn)

      // check that all records in the output are found in the input
      compareVcfs(out, missingGenotyeInfoIn)

      // get the # of phased variants from the output
      val numPhasedFromOut = getNumPhasedFromVcf(out, false)

      // check that the # of variants phased in the output agrees with the # of phased calls produced by HapCut
      val hapCutReader = HapCutReader(hapCut2In)

      val numPhasedFromHapCut = hapCutReader.flatMap(_.call).length
      numPhasedFromOut shouldBe numPhasedFromHapCut
      hapCutReader.close()

      // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
      //Files.copy(out, expectedOutput, StandardCopyOption.REPLACE_EXISTING)
      numPhasedFromOut shouldBe getNumPhasedFromVcf(expectedOutput, false)

      // check that if a variant is not phased it does not have a PS tag
      hasPhasingSetFormatTagButUnphased(out, false) shouldBe false
    }
  }

  it should "fail when IUPAC codes are found in the VCF" in {
    val out = makeTempFile("hap_cut_to_vcf.hapcut", ".vcf")
    an[Exception] should be thrownBy new HapCutToVcf(
      vcf                          = withIupacIn,
      input                        = withIupacOut,
      output                       = out,
      gatkPhasingFormat            = false,
      fixAmbiguousReferenceAlleles = false
    ).execute()
  }

  it should "succeed when IUPAC codes are found in the VCF and --fix-ambiguous-reference-alleles is specified" in {
    val out = makeTempFile("hap_cut_to_vcf.hapcut", ".vcf")

    new HapCutToVcf(
      vcf                          = withIupacIn,
      input                        = withIupacOut,
      output                       = out,
      gatkPhasingFormat            = false,
      fixAmbiguousReferenceAlleles = true
    ).execute()

    // check the reference alleles are not IUPAC
    val referenceBases = Io.readLines(out).filterNot(_.startsWith("#")).flatMap(_.split('\t')(3)).mkString("")
    referenceBases shouldBe "CCTGCCTG" // last G changed from K
  }

  // Every compressed VCF name HTSJDK accepts, including those Io does not decompress by name
  Seq(".vcf.gz", ".vcf.bgz", ".vcf.bgzf", ".vcf.gzip", ".VCF.GZ").foreach { suffix =>
    it should s"fix ambiguous reference alleles in a compressed input VCF named $suffix" in {
      val input = makeTempFile("hapcut-compressed", suffix)
      val gzip = new java.util.zip.GZIPOutputStream(java.nio.file.Files.newOutputStream(input))
      try gzip.write(java.nio.file.Files.readAllBytes(withIupacIn)) finally gzip.close()
      val output = makeTempFile("hapcut-compressed", ".vcf")
      new HapCutToVcf(vcf=input, input=withIupacOut, output=output, fixAmbiguousReferenceAlleles=true).execute()
      val bases = Io.readLines(output).filterNot(_.startsWith("#")).flatMap(_.split('\t')(3)).mkString
      bases shouldBe "CCTGCCTG"
    }
  }

  it should "write HapCut's phased genotype and genotype information onto the phased variants" in {
    Seq((hapCut1Out, HapCut1), (hapCut2Out, HapCut2)).foreach { case (hapCutOut, hapCutType) =>
      val out = makeTempFile("hap_cut_to_vcf.", ".vcf")
      new HapCutToVcf(vcf=originalVcf, input=hapCutOut, output=out).execute()

      val source  = VcfSource(out, allowMissingFieldsInHeader=true)
      val variant = source.find(_.pos == 41106449).value
      source.safelyClose()

      val gt = variant.gt("NA12878")
      gt.calls.map(_.value) should contain theSameElementsInOrderAs Seq("CT", "C") // haplotype 1 then haplotype 2
      gt.phased shouldBe true
      gt[Int]("PS") shouldBe 41106449 // replaces the input's phase set
      gt[Int]("DP") shouldBe 69 // other values on the input genotype are kept
      hapCutType match {
        case HapCut1 =>
          gt[Seq[Int]]("RC") should contain theSameElementsInOrderAs Seq(1, 0)
          gt[Seq[Float]]("LC") should contain theSameElementsInOrderAs Seq(-0.0f, -0.0f, -2.0f)
          gt[Float]("MLC") shouldBe 0f
          gt[Float]("RMEC") shouldBe 0f
        case _ =>
          gt[Int]("PR") shouldBe 0
          gt[Float]("SE") shouldBe 0f
          gt[Float]("NE") shouldBe -0.004408f
      }
    }
  }

  it should "write the trailing missing FORMAT values of each genotype" in {
    val out = makeTempFile("hap_cut_to_vcf.", ".vcf")
    new HapCutToVcf(vcf=dir.resolve("skip_prune.vcf"), input=skipPruneIn, output=out).execute()
    val records = Io.readLines(out).filterNot(_.startsWith("#")).map(_.split('\t')).toSeq
    records.foreach { fields => fields(9).split(':').length shouldBe fields(8).split(':').length }
  }

  it should "unphase the genotypes of samples other than the one HapCut phased" in {
    // add a second sample, phased with a phase set, to every record of a single-sample input
    val in = makeTempFile("hap_cut_to_vcf.", ".vcf")
    Io.writeLines(in, Io.readLines(dir.resolve("skip_prune.vcf")).map { line =>
      if (line.startsWith("##")) line
      else if (line.startsWith("#")) line + "\tOTHER"
      else {
        val fields = line.split('\t')
        val keys   = fields(8).split(':').toIndexedSeq
        val values = fields(9).split(':').toIndexedSeq.padTo(keys.length, ".")
        val (format, first) = if (keys.contains("PS")) (keys, values) else (keys :+ "PS", values :+ ".")
        val other  = format.map { case "GT" => "0|1"; case "PS" => "100"; case _ => "." }
        (fields.take(8) ++ Seq(format.mkString(":"), first.mkString(":"), other.mkString(":"))).mkString("\t")
      }
    }.toSeq)

    Iterator(true, false).foreach { gatkPhasingFormat =>
      val out = makeTempFile("hap_cut_to_vcf.", ".vcf")
      new HapCutToVcf(vcf=in, input=skipPruneIn, output=out, gatkPhasingFormat=gatkPhasingFormat).execute()

      val source   = VcfSource(out, allowMissingFieldsInHeader=true)
      val variants = source.toIndexedSeq
      source.safelyClose()

      // HapCut phased the first sample at 12902 (not the first variant in its block) but not at 178
      val phased   = variants.find(_.pos == 12902).value
      val unphased = variants.find(_.pos == 178).value
      phased.gt("SIM_INDIVIDUAL").phased shouldBe true
      phased.gt("SIM_INDIVIDUAL").get[Int]("PS").isDefined shouldBe true
      unphased.gt("SIM_INDIVIDUAL").phased shouldBe false
      Seq(phased, unphased).foreach { v =>
        v.gt("OTHER").phased shouldBe false
        v.gt("OTHER").get[Int]("PS") shouldBe None
        v.gt("OTHER").gtWithBases shouldBe s"${v.alleles.ref}/${v.alleles.alts.head}"
      }
    }
  }

  "HapCutToVcf.HapCut2GenotypeInfo" should "ignore the values of pruned, SE, and NE when they are '.'" in {
    // HapCut2 run with default options: switch error scores are blank
    HapCut2GenotypeInfo(info="0\t.\t16.00", missingAlleles=false, thresholdPruning=false).asInstanceOf[HapCut2GenotypeInfo] shouldBe HapCut2GenotypeInfo(pruned=Some(false), None, Some(16.00))

    // HapCut2 run with error analysis: switch error scores are present
    HapCut2GenotypeInfo(info="0\t1.00\t16.00", missingAlleles=false, thresholdPruning=false).asInstanceOf[HapCut2GenotypeInfo] shouldBe HapCut2GenotypeInfo(pruned=Some(false), Some(1.00), Some(16.00))

    // HapCut2 with skip prune: all info are blank
    HapCut2GenotypeInfo(info=".\t.\t.", missingAlleles=false, thresholdPruning=false).asInstanceOf[HapCut2GenotypeInfo] shouldBe HapCut2GenotypeInfo(pruned=None, None, None)
  }
}
