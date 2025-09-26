---
title: MakeTwoSampleMixtureVcf
---

# MakeTwoSampleMixtureVcf

## Overview
**Group:** VCF/BCF

Creates a simulated tumor or tumor/normal VCF by in-silico mixing genotypes from two samples.

The tumor genotypes are created by mixing the incoming genotypes for the for the
'tumor' sample and the incoming genotypes for the 'normal' samples with the 'tumor'
alleles accounting for `tumorFraction` of the resulting mixture and the 'normal' alleles
accounting for `1 - tumorFraction` of the resulting mixture (see `--tumor-fraction`).  E.g. if the 'tumor' genotype
is `A/C`, the 'normal' genotype is `C/C`, and 'tumorFraction' is set at 0.5, then the resulting
tumor genotype will be `A/C` with an allele fraction of `0.75`.  The resulting allele fraction
is written to the `AF` info field in the VCF.

In tumor-only mode (`--tumor-only`) only tumor genotypes are output.  In tumor/normal mode genotypes for
the 'normal' samples are also emitted, and match the genotypes from the input sample.

All loci (potentially restricted by intervals) that are variant in one or both samples are written
to the output VCF, though in several cases the variants will be filtered:

  * If either of the tumor or normal sample is no-called the resulting locus will have the
    `unknown_gt` filter applied
  * If the tumor and the normal have more than one alternative allele between them the
    `multi_allelic` filter will be applied
  * In tumor/normal mode (as opposed to tumor-only) loci that have an alt allele in the normal
    sample will have the `alt_allele_in_normal` filter applied

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToVcf|Input VCF file.|Required|1||
|output|o|PathToVcf|Output VCF file.|Required|1||
|tumor|t|String|Name of the 'tumor' sample in the input VCF.|Required|1||
|normal|n|String|Name of the 'normal' sample in the input VCF.|Required|1||
|tumor-fraction|f|Double|What fraction of the mixture comes from the 'tumor' sample.|Optional|1|0.5|
|tumor-only|T|Boolean|Tumor only mode - only output tumor genotypes and don't filter sites.|Optional|1|false|
|no-call-is-hom-ref|N|Boolean|Treat no-calls for either sample as hom-ref genotypes.|Optional|1|true|
|intervals|l|PathToIntervals|Optional set of intervals to restrict to.|Optional|1||

