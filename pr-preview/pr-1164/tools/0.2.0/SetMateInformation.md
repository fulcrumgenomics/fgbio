---
title: SetMateInformation
---

# SetMateInformation

## Overview
**Group:** SAM/BAM

Adds and/or fixes mate information on paired-end reads. Sets the MQ (mate mapping quality),
`MC` (mate cigar string), ensures all mate-related flag fields are set correctly, and that
the mate reference and mate start position are correct.

Supplementary records are handled correctly (updated with their mate's non-supplemental
attributes).  Secondary alignments are passed through but are not updated.

The input file must be query-name sorted or query-name grouped (i.e. all records from the same
query sequence must be adjacent in the file, though the ordering between queries is unspecified).

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM/BAM/CRAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|Output SAM/BAM/CRAM file.|Optional|1|/dev/stdout|
|ref|r|PathToFasta|Reference fasta, only needed if writing CRAM.|Optional|1||
|allow-missing-mates|x|Boolean|If specified, do not fail when reads marked as paired are missing their mate pairs.|Optional|1|false|

