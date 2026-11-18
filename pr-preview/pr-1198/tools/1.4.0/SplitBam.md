---
title: SplitBam
---

# SplitBam

## Overview
**Group:** SAM/BAM

Splits a BAM into multiple BAMs, one per-read group (or library).

The resulting BAMs will be named `<output-prefix>.<read-group-id>.bam`, or `<output-prefix>.<library-name>.bam`
when splitting by the library.  All reads without a read group, or without a library when splitting by library,
will be written to `<output-prefix>.unknown.bam`.  If no such reads exist, then no such file will exist.

By default, async writing of BAM files is controlled by the `--async-io` common tool option to
increase performance.  If the input BAM has significantly more read groups (or libraries) than your system has CPUs
it is recommended to disable this feature for this tool using `--no-async-writing`.  Asynchronous reading is not
affected.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM file.|Required|1||
|output|o|PathPrefix|Output prefix for all SAM or BAM files (ex. output/sample-name).|Required|1||
|split-by|s|SplitType|Split by library instead of read group|Optional|1|ReadGroup|
|unknown|u|String|The name to use for the unknown file|Optional|1|unknown|
|no-async-writing||Boolean|Do not write the records asynchronously. Use this to reduce memory usage when many read groups/libraries are present.|Optional|1|false|

