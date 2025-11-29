---
title: AssessPhasing
---

# AssessPhasing

## Overview
**Group:** VCF/BCF

Assess the accuracy of phasing for a set of variants.

All phased genotypes should be annotated with the `PS` (phase set) `FORMAT` tag, which by convention is the
position of the first variant in the phase set (see the VCF specification).  Furthermore, the alleles of a phased
genotype should use the `|` separator instead of the `/` separator, where the latter indicates the genotype is
unphased.

The input VCFs are assumed to be single sample: the genotype from the first sample is used.

Only bi-allelic heterozygous SNPs are considered.

The input known phased variants can be subsetted using the known interval list, for example to keep only variants
from high-confidence regions.

If the intervals argument is supplied, only the set of chromosomes specified will be analyzed.  Note that the full
chromosome will be analyzed and start/stop positions will be ignored.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|called-vcf|c|PathToVcf|The VCF with called phased variants.|Required|1||
|truth-vcf|t|PathToVcf|The VCF with known phased variants.|Required|1||
|output|o|PathPrefix|The output prefix for all output files.|Required|1||
|known-intervals|k|PathToIntervals|The interval list over which known phased variants should be kept.|Optional|1||
|allow-missing-fields-in-vcf-header|m|Boolean|Allow missing fields in the VCF header.|Optional|1|true|
|skip-mismatching-alleles|s|Boolean|Skip sites where the truth and call are both called but do not share the same alleles.|Optional|1|true|
|intervals|l|PathToIntervals|Analyze only the given chromosomes in the interval list.  The entire chromosome will be analyzed (start and end ignored).|Optional|1||
|modify-blocks|b|Boolean|Remove enclosed phased blocks and truncate overlapping blocks.|Optional|1|true|
|debug-vcf|d|Boolean|Output a VCF with the called variants annotated by if their phase matches the truth|Optional|1|false|

