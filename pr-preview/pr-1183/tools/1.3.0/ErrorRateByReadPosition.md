---
title: ErrorRateByReadPosition
---

# ErrorRateByReadPosition

## Overview
**Group:** SAM/BAM

Calculates the error rate by read position on coordinate sorted mapped BAMs. The output file contains
a row per read (first of pair, second of pair and unpaired), per position in read, with the total number
of bases observed, the number of errors observed, the overall error rate, and the rate of each kind of
substitution error.

Substitution types are collapsed based on the reference or expected base, with only six substitution
types being reported: `A>C`, `A>G`, `A>T`, `C>A`, `C>G` and `C>T`.  For example, `T>G` is grouped in
with `A>C`.

Analysis can be restricted to a set of intervals via the `--intervals` option. Genomic positions can be
excluded from analysis by supplying a set of variants (either known variants in the sample or a catalog
of known variants such as dbSNP).  For data believed to have low error rates it is recommended to use
both the `--intervals` and `--variants` options to restrict analysis to only regions expected to be
homozygous reference in the data.

The following are reads / bases are excluded from the analysis:

- Unmapped reads
- Reads marked as failing vendor quality
- Reads marked as duplicates (unless `--include-duplicates` is specified)
- Secondary and supplemental records
- Soft-clipped bases in records
- Reads with MAPQ < `--min-mapping-quality` (default: 20)
- Bases with base quality < `--min-base-quality` (default: 0)
- Bases where either the read base or the reference base is non-ACGT

An output text file is generated with the extension `.error_rate_by_read_position.txt`

If R's `Rscript` utility is on the path and `ggplot2` is installed in the R distribution then a PDF
of error rate plots will also be generated with extension `.error_rate_by_read_position.pdf`.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathPrefix|Output metrics prefix. If not given, will use the input BAM basename.|Optional|1||
|ref|r|PathToFasta|Reference sequence fasta file.|Required|1||
|variants|v|PathToVcf|Optional file of variant sites to ignore.|Optional|1||
|intervals|l|PathToIntervals|Optional list of intervals to restrict analysis to.|Optional|1||
|include-duplicates|d|Boolean|Include duplicate reads, otherwise ignore.|Optional|1|false|
|min-mapping-quality|m|Int|The minimum mapping quality for a read to be included.|Optional|1|20|
|min-base-quality|q|Int|The minimum base quality for a base to be included.|Optional|1|0|
|collapse||Boolean|Collapse substitution types based on the reference or expected base, with only six substitution types being reported: `A>C`, `A>G`, `A>T`, `C>A`, `C>G` and `C>T`.For example, `T>G` is grouped in with `A>C`. Otherwise, all possible substitution types will be reported.|Optional|1|true|

