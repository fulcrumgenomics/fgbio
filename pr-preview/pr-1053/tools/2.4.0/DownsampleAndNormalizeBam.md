---
title: DownsampleAndNormalizeBam
---

# DownsampleAndNormalizeBam

## Overview
**Group:** SAM/BAM

Downsamples a BAM in a biased way to a uniform coverage across regions.

Attempts to downsample a BAM such that every base in the genome (or in the target `regions` if provided)
is covered by at least `coverage` reads.  When computing coverage:
  - Reads marked as secondary, duplicate or unmapped are not used
  - A base can receive coverage from only one read with the same queryname (i.e. mate overlaps are not counted)
  - Coverage is counted if a read _spans_ a base, even if that base is deleted in the read

Reads are first sorted into a random order (by hashing read names).  Reads are then consumed one template
at a time, and if any read adds coverage to base that is under the target coverage, _all_ reads (including
secondary, unmapped, etc.) for that template are emitted into the output.

Given the procedure used for downsampling, it is likely the output BAM will have coverage up to 2X the requested
coverage at regions in the input BAM that are i) well covered and ii) are close to regions that are poorly
covered.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input SAM or BAM file.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Required|1||
|coverage|c|Int|Desired minimum coverage.|Required|1||
|min-map-q|m|Int|Minimum mapping quality to count a read as covering.|Optional|1|0|
|seed|s|Int|Random seed to use when randomizing order of reads/templates.|Optional|1|42|
|regions|l|PathToIntervals|Optional set of regions for coverage targeting.|Optional|1||
|max-in-memory|M|Int|Maximum records to be held in memory while sorting.|Optional|1|1000000|

