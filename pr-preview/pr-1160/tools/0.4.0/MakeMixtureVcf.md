---
title: MakeMixtureVcf
---

# MakeMixtureVcf

## Overview
**Group:** VCF/BCF

Creates a VCF with one sample whose genotypes are a mixture of other samples'.

The input VCF must contain all samples to be mixed, and may optionally contain other samples.
Sample mixtures can be specified in one of several ways:

  1. `--samples s1 s2 s3`: specifies that the samples with names `s1`, `s`2 and `s3` should be mixed equally.
  2. ` --samples s1@0.1 s2@0.1 s3@0.8`: specifies that the three samples should be mixed at `0.1`, `0.1`, and `0.8`
                                        respectively
  3. `--samples s1@0.1 s2 s3`: specifies that `s1` should form `0.1` of the mixture and that the remaining (`0.9`)
                               should be equally split amongst `s2` and s`3`
  4. If no sample names are given, all samples in the input VCF will be used, mixing equally

The input samples are assumed to be diploid with the allele fraction defined by the genotype (`0/0`, `0/1` or `1/1`).
The `--allele-fraction-field` option may be specified, in which case the allele fraction of the input genotypes
will be retrieved from the specified `FORMAT` field in the VCF genotypes.  All genotypes except hom-ref and
no-call genotypes must have this field present if supplied.

If the --`no-call-is-hom-ref` flag is `true` (the default) then no-call genotypes in the input VCF are interpreted
as hom-ref genotypes.  If it is `false`, any location with a no-call genotype will be emitted as a no-call in
the mixture, and will be filtered.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToVcf|Input VCF containing genotypes for the samples to be mixed.|Required|1||
|output|o|PathToVcf|Output VCF of mixture sample.|Required|1||
|samples|s|String|Samples to mix. See general usage for format and examples.|Optional|Unlimited||
|output-sample-name|S|String|Output sample name.|Optional|1|mixture|
|no-call-is-hom-ref|N|Boolean|Treat no-calls for samples as hom-ref genotypes.|Optional|1|true|
|allele-fraction-field|a|String|Format field containing allele fraction.|Optional|1||
|precision|p|Int|Digits of precision in generated allele fractions.|Optional|1|5|

