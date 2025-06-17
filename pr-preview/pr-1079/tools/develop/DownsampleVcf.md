---
title: DownsampleVcf
---

# DownsampleVcf

## Overview
**Group:** VCF/BCF

Re-genotypes a VCF after downsampling the allele counts.

The input VCF must have at least one sample.

If the input VCF contains a single sample, the downsampling target may be specified as a
proportion of the original read depth using `--proportion=(0..1)`, or as the combination of
the original and target _number of sequenced bases_ (`--originalBases` and
`--downsampleToBases`). For multi-sample VCFs, the downsampling target must be specified using
`--downsampleToBases`, and a metadata file with the total number of sequenced bases per sample
is required as well. The metadata file must follow the
[[https://www.internationalgenome.org/category/meta-data/] 1000 Genomes index format], but the
only required columns are `SAMPLE_NAME` and `BASE_COUNT`. A propportion for each sample is
calculated by dividing the _target number of sequenced bases_ by the _original number of
sequenced bases_.

The tool first (optionally) winnows the VCF file to remove variants within a distance to each
other specified by `--window-size` (the default value of `0` disables winnowing). Next, each
sample at each variant is examined independently. The allele depths per-genotype are randoml
downsampled given the proportion. The downsampled allele depths are then used to re-compute
allele likelhoods and produce a new genotype.

The tool outputs a downsampled VCF file with the winnowed variants removed, and with the
genotype calls and `DP`, `AD`, and `PL` tags updated for each sample at each retained variant.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToVcf|The vcf to downsample.|Required|1||
|proportion|p|Double|Proportion of bases to retain (for single-sample VCF).|Optional|1||
|original-bases|b|Double|Original number of bases (for single-sample VCF).|Optional|1||
|metadata|m|FilePath|Index file with bases per sample.|Optional|1||
|downsample-to-bases|n|Double|Target number of bases to downsample to.|Optional|1||
|output|o|PathToVcf|Output file name.|Required|1||
|window-size|w|Int|Winnowing window size.|Optional|1|0|
|epsilon|e|Double|Sequencing Error rate for genotyping.|Optional|1|0.01|
|min-ad-homvar|v|Int|Minimum allele depth to call HOMVAR (otherwise NO-CALL)|Optional|1|3|
|min-ad-homref|r|Int|Minimum allele depth to call HOMREF (otherwise NO-CALL)|Optional|1|1|
|write-no-call|c|Boolean|True to write out no-calls.|Optional|1|false|
|seed|s|Int|Random seed value.|Optional|1|42|

