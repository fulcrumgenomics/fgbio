---
title: FilterBam
---

# FilterBam

## Overview
**Group:** SAM/BAM

Filters reads out of a BAM file. Removes reads that may not be useful in downstream processing, in order
to reduce the size of the file. By default will remove unmapped reads, reads with MAPQ=0, reads
marked as secondary alignments, reads marked as duplicates, and if a set of Intervals are provided,
reads that do not overlap any of the intervals.

NOTE: this will usually produce a BAM file in which some mate-pairs are orphaned (i.e. read 1 or
read 2 is included, but not both), but does not update any flag fields.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|intervals|l|PathToIntervals|Optionally remove reads not overlapping intervals.|Optional|1||
|remove-duplicates|D|Boolean|If true remove all reads that are marked as duplicates.|Optional|1|true|
|remove-unmapped-reads|U|Boolean|Remove all unmapped reads.|Optional|1|true|
|min-map-q|M|Int|Remove all mapped reads with MAPQ lower than this number.|Optional|1|1|
|remove-secondary-alignments|S|Boolean|Remove all reads marked as secondary alignments.|Optional|1|true|

