package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.vcf.api.VcfSource

/** Namespace for common VCF and variant call utilities. */
object VcfUtil {

  /** Return the only sample in this VCF source otherwise raise an exception. */
  def onlySample(source: VcfSource): String = {
    validateHasSingleSample(source)
    source.header.samples.head
  }

  /** Validate that a VCF source has a single sample's genotypes.  */
  def validateHasSingleSample(source: VcfSource): Unit = {
    val names = source.header.samples
    val n     = names.length
    require(n == 1, s"VCF is not single-sample and has $n samples: ${names.mkString(", ")}")
  }
}
