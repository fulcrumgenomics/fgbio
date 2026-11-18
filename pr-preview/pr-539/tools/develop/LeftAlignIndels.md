---
title: LeftAlignIndels
---

# LeftAlignIndels

## Overview
**Group:** SAM/BAM

Left-aligns indels in a SAM/BAM file.

An indel is shifted left only across repeated bases, such that every aligned read base remains paired with an
identical reference base.  Indels are not shifted into soft-clipped bases or across reference skips (`N` CIGAR
operators); each block of the alignment between reference skips is left-aligned independently.  Insertions and
deletions that end up adjacent to each other are written with the insertion first (e.g. `2I1D`).  Indels may be
shifted to the start of an alignment block, producing a leading insertion or deletion.

Unmapped reads, reads without stored bases (`SEQ` of `*`), and reads with padding (`P`) CIGAR operators are left
unchanged.

The reference bases are retrieved from the reference FASTA (`--ref`) when given.  An indexed FASTA (with a `.fai`)
may be used with any input sort order, while a FASTA without an index may only be used with coordinate sorted
input.  Otherwise, the reference bases are reconstructed from each read's MD tag, which is then required for every
mapped read that may be changed.

When a read's CIGAR changes, its MD and NM tags are re-computed if present.  The CIGAR uses match (`=`) and
mismatch (`X`) operators when `--use-equals-and-x` is given or the read's CIGAR already uses them, otherwise
match-or-mismatch (`M`) operators.

Mate information (e.g. the MC tag and insert size) is reset for all primary and supplementary reads, the MC tag of
secondary reads is updated with their mate's primary CIGAR, and the CIGARs in SA tags are updated.  To do so, the
input is re-sorted into query name order when not already query name sorted or grouped, then re-sorted back into
the input's sort order (if known) when written.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|ref|r|PathToFasta|Reference sequence FASTA file, otherwise the MD tag is required.|Optional|1||
|allow-missing-mates|x|Boolean|If specified, do not fail when reads marked as paired are missing their mate pairs.|Optional|1|false|
|use-equals-and-x|X|Boolean|Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).|Optional|1|false|

