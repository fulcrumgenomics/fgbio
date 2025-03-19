---
title: EstimatePoolingFractions
---

# EstimatePoolingFractions

## Overview
**Group:** SAM/BAM

Examines sequence data generated from a pooled sample and estimates the fraction of sequence data
coming from each constituent sample. Uses a VCF of known genotypes for the samples within the
mixture along with a BAM of sequencing data derived from the pool.  Performs a multiple regression
for the alternative allele fractions at each SNP locus, using as inputs the individual sample's genotypes.
Only SNPs that are bi-allelic within the pooled samples are used.

Each sample's contribution of REF vs. ALT alleles at each site is derived in one of two ways: (1) if
the sample's genotype in the VCF has the `AF` attribute then the value from that field will be used, (2) if the
genotype has no `AF` attribute then the contribution is estimated based on the genotype (e.g. 0/0 will be 100%
ref, 0/1 will be 50% ref and 50% alt, etc.).

Various filtering parameters can be used to control which loci are used:

- _--intervals_ will restrict analysis to variants within the described intervals
- _--min-genotype-quality_ will filter out any site with any genotype with GQ < n
- _--min-mean-sample-coverage_ requires that the coverage of a site in the BAM be >= `min-mean-sample-coverage * n_samples`
- _--min-mapping-quality_ filters out reads in the BAM with MQ < n
- _--min-base-quality_ filters out bases in the BAM with Q < n

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|vcf|v|PathToVcf|VCF of individual sample genotypes.|Required|1||
|bam|b|PathToBam|Path to BAM file of sequencing data.|Required|1||
|output|o|FilePath|Output file to write with pooling fractions.|Required|1||
|intervals|l|PathToIntervals|Zero or more set of regions to restrict analysis to.|Optional|Unlimited||
|samples|s|String|Optional subset of samples from VCF to use.|Optional|Unlimited||
|non-autosomes|n|String|Non-autosomal chromosomes to avoid.|Required|Unlimited|M, chrM, MT, X, chrX, Y, chrY|
|min-genotype-quality|g|Int|Minimum genotype quality. Use -1 to disable.|Optional|1|30|
|min-mean-sample-coverage|c|Int|Minimum (sequencing coverage @ SNP site / n_samples).|Optional|1|6|
|min-mapping-quality|m|Int|Minimum mapping quality.|Optional|1|20|
|min-base-quality|q|Int|Minimum base quality.|Optional|1|5|

