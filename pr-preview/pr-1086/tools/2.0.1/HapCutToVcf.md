---
title: HapCutToVcf
---

# HapCutToVcf

## Overview
**Group:** VCF/BCF

Converts the output of `HAPCUT` (`HapCut1`/`HapCut2`) to a VCF.

The output of `HAPCUT` does not include all input variants, but simply those variants that are in phased blocks.
This tool takes the original VCF and the output of `HAPCUT`, and produces a VCF containing both the variants that
were phased and the variants that were not phased, such that all variants in the original file are in the output.

The original VCF provided to `HAPCUT` is assumed to contain a single sample, as `HAPCUT` only supports a single
sample.

By default, all phased genotypes are annotated with the `PS` (phase set) `FORMAT` tag, which by convention is the
position of the first variant in the phase set (see the VCF specification).  Furthermore, this tool formats the
alleles of a phased genotype using the `|` separator instead of the `/` separator, where the latter indicates the
genotype is unphased.  If the option to output phased variants in GATK's `ReadBackedPhasing` format is used, then
the first variant in a phase set will have `/` instead of `|` separating its alleles.  Also, all phased variants
will have `PASS` set in its `FILTER` column, while unphased variants will have `NotPhased` set in their `FILTER`
column.  Unlike GATK's `ReadBackedPhasing`, homozygous variants will always be unphased.

More information about the purpose and operation of GATK's Read-backed phasing, including its output format, can
be found in the [GATK Forum](http://gatkforums.broadinstitute.org/gatk/discussion/45/purpose-and-operation-of-read-backed-phasing)

Additional `FORMAT` fields for phased variants are provided corresponding to per-genotype information produced by
`HapCut1`:

  1. The `RC` tag gives the counts of calls supporting `allele0` and `allele1` respectively.
  2. The `LC` tag gives the change in likelihood if this SNP is made homozygous or removed.
  3. The `MCL` tag gives the maximum change in likelihood if this SNP is made homozygous or removed.
  4. The `RMEC` tag gives the reduction in `MEC` score if we remove this variant altogether.

Additional `FORMAT` fields for phased variants are provided corresponding to per-genotype information produced by
`HapCut2`:

  1. The `PR` tag is `1` if `HapCut2` pruned this variant, `0` otherwise.
  2. The `SE` tag gives the confidence (`log10`) that there is not a switch error occurring immediately before
     the SNV (closer to zero means that a switch error is more likely).
  3. The `NE` tag gives the confidence (`log10`) that the SNV is not a mismatch (single SNV) error (closer to
     zero means a mismatch is more likely).

HapCut2 should not be run with `--call-homozygous` as the genotypes may be different than the input and so is not
currently supported.

For more information about `HapCut1`, see the source code or paper below.
  source code: https://github.com/vibansal/hapcut
  HapCut1 paper: An efficient and accurate algorithm for the haplotype assembly problem Bioinformatics. 2008 Aug
    15;24(16):i153-9.

For more information about `HapCut2`, see the [source code](https://github.com/vibansal/HapCUT2) below.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|vcf|v|PathToVcf|The original VCF provided to `HapCut1`/`HapCut2`.|Required|1||
|input|i|FilePath|The file produced by `HapCut1`/`HapCut2`.|Required|1||
|output|o|PathToVcf|The output VCF with both phased and unphased variants.|Required|1||
|gatk-phasing-format|r|Boolean|Output phased variants in GATK's `ReadBackedPhasing` format.|Optional|1|false|
|fix-ambiguous-reference-alleles||Boolean|Fix IUPAC codes in the original VCF to be VCF 4.3 spec-compliant (ex 'R' -> 'A').  Does not support BCF inputs.|Optional|1|false|

