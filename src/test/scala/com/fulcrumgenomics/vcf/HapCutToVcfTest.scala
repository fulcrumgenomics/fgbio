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

import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Io
import dagr.commons.io.PathUtil
import dagr.commons.CommonsDef._
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.JavaConversions._

/**
  * Tests for HapCutToVcf.
  */
class HapCutToVcfTest extends UnitSpec {

  private val dir            = PathUtil.pathTo("src/test/resources/com/fulcrumgenomics/vcf")
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

  private def countVcfRecords(vcf: PathToVcf): Int = {
    val vcfReader = new VCFFileReader(vcf.toFile, false)
    yieldAndThen(vcfReader.iterator().length)(vcfReader.close())
  }

  private def compareVcfs(newVcf: PathToVcf, originalVcf: PathToVcf): Unit = {
    val newVcfReader      = new VCFFileReader(newVcf.toFile, false)
    val originalVcfReader = new VCFFileReader(originalVcf.toFile, false)

    for (newVariantCtx <- newVcfReader) {
      originalVcfReader.exists { originalVariantCtx =>
        originalVariantCtx.getContig == newVariantCtx.getContig &&
        originalVariantCtx.getStart  == newVariantCtx.getStart &&
        originalVariantCtx.getEnd    == newVariantCtx.getEnd
      } shouldBe true
    }
  }

  private def isPhased(ctx: VariantContext, gatkPhasingFormat: Boolean): Boolean = {
    if (gatkPhasingFormat) ctx.isNotFiltered // are marked as passed filter
    else ctx.getGenotypes.exists(_.isPhased) // are marked as phased
  }

  private def getNumPhasedFromVcf(path: PathToVcf, gatkPhasingFormat: Boolean): Int = {
    val vcfReader = new VCFFileReader(path.toFile, false)
    val numPhased = vcfReader.iterator().count { ctx => isPhased(ctx, gatkPhasingFormat) }
    vcfReader.close()
    numPhased
  }

  private def hasPhasingSetFormatTagButUnphased(path: PathToVcf, gatkPhasingFormat: Boolean): Boolean = {
    val vcfReader = new VCFFileReader(path.toFile, false)
    val hasPhasingSetTag = vcfReader
      .iterator()
      .filterNot { ctx => isPhased(ctx, gatkPhasingFormat) }
      .exists { ctx => ctx.getGenotypes.exists(_.hasExtendedAttribute(HapCut1VcfHeaderLines.PhaseSetFormatTag)) }
    vcfReader.close()
    hasPhasingSetTag
  }

  "HapCutReader" should "read in a HAPCUT1 file" in {
    val reader = HapCutReader(hapCut1Out)
    val allCalls = reader.toSeq
    val calls = allCalls.flatMap(_.call)
    allCalls.length shouldBe 342 // 342 total variants
    calls.length shouldBe 237 // 237 phased variants
    calls.map(_.phaseSet).distinct.length shouldBe 8 // eight phased blocks
    reader.close()
  }

  it should "read in a HAPCUT2 file" in {
    val reader = HapCutReader(hapCut2Out)
    val allCalls = reader.toSeq
    val calls = allCalls.flatMap(_.call)
    allCalls.length shouldBe 342 // 342 total variants
    calls.length shouldBe 237 // 237 phased variants
    calls.map(_.phaseSet).distinct.length shouldBe 8 // eight phased blocks
    reader.close()
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
    Stream(true, false).foreach { gatkPhasingFormat =>
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
    }
  }

  it should "convert a HAPCUT2 file to VCF in both GATK and VCF-spec phasing format" in {
    Stream(true, false).foreach { gatkPhasingFormat =>
      val expectedOutput = if (gatkPhasingFormat) hapCut2GatkVcf else hapCut2Vcf
      val out = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

      new HapCutToVcf(
        vcf               = originalVcf,
        input             = hapCut2Out,
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
      val hapCutReader = HapCutReader(hapCut2Out)
      val numPhasedFromHapCut = hapCutReader.flatMap(_.call).length
      numPhasedFromOut shouldBe numPhasedFromHapCut
      hapCutReader.close()

      // check that the # of variants phased in the output agrees with the # of phased calls in the expected output
      numPhasedFromOut shouldBe getNumPhasedFromVcf(expectedOutput, gatkPhasingFormat)

      // check that if a variant is not phased it does not have a PS tag
      hasPhasingSetFormatTagButUnphased(out, gatkPhasingFormat) shouldBe false
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
    Stream(true, false).foreach { gatkPhasingFormat =>
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

  it should "fail when there are missing variants in the input VCF" in {
    val out       = makeTempFile("hap_cut_to_vcf.hapcut2", ".vcf")

    an[Exception] should be thrownBy new HapCutToVcf(
      vcf               = missingVariantsIn,
      input             = missingVariantsOut,
      output            = out,
      gatkPhasingFormat = false
    ).execute()
  }
}
