package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef.SafelyClosable
import com.fulcrumgenomics.testing.{UnitSpec, VcfBuilder}

/** Unit tests for [[VcfSource]]. */
class VcfSourceTest extends UnitSpec {

  /** The name for test sample number 1. */
  private val sample1 = "sample1"

  /** The name for test sample number 2. */
  private val sample2 = "sample2"

  "VcfUtil.onlySample" should "return the only sample in a VCF source" in {
    val singleSample = VcfBuilder(samples = Seq(sample1)).toSource
    val doubleSample = VcfBuilder(samples = Seq(sample1, sample2)).toSource
    VcfSource.onlySample(singleSample) shouldBe sample1
    singleSample.safelyClose()
    an[IllegalArgumentException] shouldBe thrownBy { VcfSource.onlySample(doubleSample) }
  }
}
