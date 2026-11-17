---
title: LeftAlignIndels
---

# LeftAlignIndels

## Overview
**Group:** SAM/BAM

Left-aligns indels in a SAM/BAM file.

This tool is preferred over GATK's LeftAlignIndels, as the latter can only left-align reads with at most one indel.

Mate information will be reset.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|ref|r|PathToFasta|Reference sequence FASTA file, otherwise the MD tag is required.|Optional|1||
|allow-missing-mates|x|Boolean|If specified, do not fail when reads marked as paired are missing their mate pairs.|Optional|1|false|
|use-equals-and-x|X|Boolean|Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).|Optional|1|false|

