---
title: ClipOverlappingReads
---

# ClipOverlappingReads

## Overview
**Group:** SAM/BAM

DEPRECATED: Use `ClipBam` instead.

Clips reads from the same template to eliminate overlap between the reads. Ensures that downstream
processes, particularly variant calling, cannot double-count evidence from the same template when
both reads span a variant site in the same template.

Clipping is only performed on `FR` read pairs, and is implemented by clipping approximately half the
overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be
substituted using the `--soft-clip` parameter.

Secondary alignments and supplemental alignments are not clipped, but are passed through into the
output.

If the input BAM is neither `queryname` sorted nor `query` grouped, it will be sorted into queryname
order so that clipping can be performed on both ends of a pair simultaneously and so that mate
pair information can be reset across all reads for the template.  Post-clipping the reads are
resorted into coordinate order, any existing `NM`, `UQ` and `MD` tags are repaired, and the output is
written in coordinate order.

DEPRECATED: Use `ClipBam` instead.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM file of aligned reads in coordinate order.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Required|1||
|ref|r|PathToFasta|Reference sequence fasta file.|Required|1||
|soft-clip|s|Boolean|Soft clip reads instead of hard clipping.|Optional|1|false|
|auto-clip-attributes|a|Boolean|Automatically clip extended attributes that are the same length as bases.|Optional|1|false|

