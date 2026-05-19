---
title: TrimFastq
---

# TrimFastq

## Overview
**Group:** FASTQ

Trims reads in one or more line-matched fastq files to a specific read length. The
individual fastq files are expected to have the same set of reads, as would be the
case with an `r1.fastq` and `r2.fastq` file for the same sample.

Optionally supports dropping of reads across all files when one or more reads
is already shorter than the desired trim length.

Input and output fastq files may be gzipped.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToFastq|One or more input fastq files.|Required|Unlimited||
|output|o|PathToFastq|A matching number of output fastq files.|Required|Unlimited||
|length|l|Int|Length to trim reads to.|Required|1||
|exclude|x|Boolean|Exclude reads below the trim length.|Optional|1|false|

