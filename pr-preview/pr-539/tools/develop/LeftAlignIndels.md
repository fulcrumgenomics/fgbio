---
title: LeftAlignIndels
---

# LeftAlignIndels

## Overview
**Group:** SAM/BAM

Left-aligns indels in a SAM/BAM file.

An indel is shifted left only across repeated bases, such that every aligned read base remains paired with an
identical reference base.  Indels are not shifted into clipped bases or across reference skips (`N` CIGAR
operators); each block of the alignment between reference skips is left-aligned independently.  Within a run of
adjacent insertions and deletions, the insertions are written first (e.g. `2I1D`), so alignments that differ only
in the order of adjacent insertions and deletions are written identically.  Shifting may place an insertion directly
next to a deletion where the input had aligned bases between them; GATK's `GOOD_CIGAR` read filter (on by default in
HaplotypeCaller and Mutect2) removes reads with such adjacent indels.

An indel at the start of the alignment (after any clipping), whether shifted there or already there in the input,
is not supported by any aligned bases before it, and GATK's `GOOD_CIGAR` read filter removes reads that start with a
deletion.  `--leading-indel-mode` controls what is done with such an indel:

- `Trim` (the default): the indel is removed.  Inserted bases are soft-clipped, and the alignment start is moved
  past deleted bases.  The MD, NM, and UQ tags are re-computed if present.
- `Keep`: the indel is kept at the start of the alignment.
- `Filter`: all records of the read's template (the read, its mate, and their secondary and supplementary records)
  are removed from the output, so that no remaining record refers to a removed one.  A read with an unavailable
  query name (`*`) is not grouped with other reads, so only the read itself is removed.

Indels at the start of a later block of the alignment (after a reference skip) are always kept.

Unmapped reads, reads without stored bases (`SEQ` of `*`), and reads with padding (`P`) CIGAR operators are left
unchanged.  Read bases stored as `=` (equal to the reference) are compared using the reference base, and are written
as the reference base when the read's CIGAR changes.

The reference bases are retrieved from the reference FASTA (`--ref`) when given.  An indexed FASTA (with a `.fai`)
may be used with any input sort order, while a FASTA without an index may only be used with coordinate sorted
input whose sequence dictionary lists the contigs in the same order as the FASTA.  The FASTA contig lengths must
match the input's sequence dictionary.  Without `--ref`, the reference bases are reconstructed from each read's MD
tag, which is then required for every mapped read that may be changed, and such reads may not have `=` read bases
that the MD tag reports as matching the reference.

When a read's alignment changes, its MD, NM, and UQ tags are re-computed if present, counting a base as a mismatch
unless the read and reference bases are the same and not `N` (as samtools does).  UQ is removed if the read has no
base qualities.  The alignment score (AS) is not changed, even when indels are merged or trimmed.  The CIGAR uses
match (`=`) and mismatch (`X`) operators, with the same definition of a mismatch, when `--use-equals-and-x` is
given or the read's CIGAR already uses them, otherwise match-or-mismatch (`M`) operators.  A read whose CIGAR
already uses `=` and `X` operators keeps its CIGAR when no indel is moved.

Mate information (e.g. the MC tag and insert size) is reset for all primary and supplementary reads, including
those whose alignment did not change.  The mate start, insert size, and MC tag of secondary reads are updated from
their mate's primary alignment, and the positions and CIGARs in SA tags are updated.  An SA tag entry is updated only when
exactly one alignment of the template matches it, otherwise it is left unchanged and counted.  To do so, input that
is not query name sorted or grouped is sorted into query name order, then sorted back into the input's sort order
(if known) when written; input of unknown order is written in query name order.  Query name sorted or grouped input
is written in the input order, except that template-coordinate grouped input is sorted again when trimming
leading indels may have changed the template's coordinates.  Reads with an unavailable query name (`*`) are not
grouped with other reads, so their mate information and SA tags are not updated.  Unpaired reads that share a query
name with another unpaired read are written without updating their SA tags.

If the tool fails, the partially written output is removed.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|ref|r|PathToFasta|Reference sequence FASTA file, indexed unless the input is coordinate sorted.  If not given, the MD tag is used.|Optional|1||
|allow-missing-mates|x|Boolean|If specified, do not fail when reads marked as paired are missing their mate pairs.|Optional|1|false|
|use-equals-and-x||Boolean|Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).|Optional|1|false|
|leading-indel-mode|l|LeadingIndelMode|What to do with an indel at the start of a read's alignment: trim it from the alignment, keep it, or filter the read's template.|Optional|1|Trim|

