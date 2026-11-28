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
  past deleted bases.  The MD, NM, and UQ tags are re-computed if present.  The indel is kept if no aligned bases
  follow it before the first reference skip (e.g. `2I10N2M`).
- `Keep`: the indel is kept at the start of the alignment.
- `Filter`: all records of the read's template (the read, its mate, and their secondary and supplementary records)
  are removed from the output, so that no remaining record refers to a removed one.  A read with an unavailable
  query name (`*`) is not grouped with other reads, so only the read itself is removed.

Indels at the start of a later block of the alignment (after a reference skip) are always kept.

Unmapped reads, reads without stored bases (`SEQ` of `*`), and reads with padding (`P`) CIGAR operators are left
unchanged.  Read bases stored as `=` (equal to the reference) are compared using the reference base, and are written
as the reference base when the read's CIGAR changes.  If a resolved base cannot be stored unchanged in SAM/BAM
(e.g. `U`, which is normalized to `T`), the tool fails with the read's name.

When `--ref` is given, the entire reference FASTA is loaded into memory before processing reads.  Plain, gzip, and
BGZF FASTAs are supported; no index is needed, and any `.fai`, `.gzi`, or `.dict` sidecars are ignored.  Every contig
in the input's sequence dictionary must be present with the same length, even if it has no reads.  Additional FASTA
contigs (which are read but not kept in memory) and a different contig order are allowed.  A gzip or BGZF FASTA must
have a `.gz` or `.gzip` extension (e.g. not `.bgz`); a plain FASTA need not use a standard filename extension.  The
reference is loaded and checked even for an empty input or one with no alignments that need changing.

A human reference needs about 3.1 GB for its bases alone, plus memory for parsing, reads, and sorting buffers.
The launcher's default heap may be too small.  Use an explicit larger heap (e.g. `fgbio -Xmx8g ...` or
`java -Xmx8g -jar fgbio.jar ...`) as a starting point; larger references, reads, or tags may need more.

Without `--ref`, the reference bases are reconstructed from each read's MD tag, which is required for every mapped
read that may be changed.  Such reads may not have `=` read bases that the MD tag reports as matching the reference.

When a read's alignment changes, its MD, NM, and UQ tags are re-computed if present, counting a base as a mismatch
unless the read and reference bases are the same and not `N` (as samtools does).  UQ is removed if the read has no
base qualities.  The BQ, minimap2 cs, ds, and de tags are removed when the alignment edit path changes, but kept
when only M versus =/X notation changes.  The alignment score (AS) is not changed, even when indels are merged or
trimmed.  The CIGAR uses match (`=`) and mismatch (`X`) operators, with the same definition of a mismatch, when
`--use-equals-and-x` is given or the read's CIGAR already uses them, otherwise match-or-mismatch (`M`) operators.
A read whose CIGAR already uses `=` and `X` operators keeps its CIGAR when no indel is moved.

Mate information (e.g. the MC tag and insert size) is reset for all primary and supplementary reads, including
those whose alignment did not change; existing ms scores are preserved and absent ms tags are not added.  The mate
start, insert size, and MC tag of secondary reads are updated from their mate's primary alignment only when its
original reference, start, strand, and optional MC agree.  The positions, CIGARs, and edit distances in SA tags are
updated.  An SA tag entry is matched against the original start and CIGAR, accepting exact CIGARs and minimap2-style summaries
(end clips, the smaller of the aligned query and reference spans as M, then any net I or D).  Hard and soft clips,
and M versus =/X, are treated equivalently.  An entry is updated only when exactly one alignment of the template
matches either form, otherwise it is left unchanged and counted.

Records are grouped by query name before any alignments are changed.  Coordinate and template-coordinate input is
sorted into query name order, then sorted into coordinate or template-coordinate order after repairs, in every
leading indel mode.  These two sorts require additional memory and temporary disk space.  The template-coordinate
order is fgbio's, so may order ties differently than the input (e.g. if sorted by samtools).  Template-coordinate
`GO:query` is not trusted because secondary or supplementary alignments may be separated from their primary pair.
Other query name sorted input, or `SO:unsorted GO:query` input, keeps its input order and ordering metadata.  Random
input is re-sorted into random order.  Unknown or ordinary unsorted input is written in query name order, including
`SO:unknown GO:query` (or `GO:query` without `SO`): the unknown sort order takes precedence over the grouping
declaration.

With `--leading-indel-mode Trim`, trimming changes the coordinates of a read's template, so in template-coordinate
output (e.g. from `GroupReadsByUmi`) the reads of a molecule (`MI` tag) may no longer be adjacent, as required by
the consensus callers.  Left-align indels before grouping reads by UMI, or use `Keep` or `Filter`, which do not
change the coordinates of the reads that are written.

For template-coordinate output, every retained paired record with a mapped mate must have an MC tag after repairs,
including in `Keep` and `Filter` modes.  If one does not, the tool fails with a diagnostic instead of changing the
output order.  Add mate CIGARs (e.g. with `samtools fixmate -m` or fgbio `SetMateInformation` on query name grouped
input), or query name sort the input to obtain query name output.  `--allow-missing-mates` does not waive this
requirement.

Reads with an unavailable query name (`*`) are not grouped with other reads, so their mate information and SA tags
are not updated.  Unpaired reads that share a query name with another unpaired read are written without updating
their SA tags, and with `--leading-indel-mode Filter` only those reads with a leading indel are removed.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|ref|r|PathToFasta|Reference sequence FASTA file, loaded entirely into memory without requiring an index.  If not given, the MD tag is used.|Optional|1||
|allow-missing-mates|x|Boolean|If specified, do not fail when reads marked as paired are missing their mate pairs.|Optional|1|false|
|use-equals-and-x||Boolean|Use match (=) and mismatch (X) SAM CIGAR operators instead of match-or-mismatch (M).|Optional|1|false|
|leading-indel-mode|l|LeadingIndelMode|What to do with an indel at the start of a read's alignment: trim it from the alignment, keep it, or filter the read's template.|Optional|1|Trim|

