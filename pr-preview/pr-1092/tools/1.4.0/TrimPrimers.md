---
title: TrimPrimers
---

# TrimPrimers

## Overview
**Group:** SAM/BAM

Trims primers from reads post-alignment.  Takes in a BAM file of aligned reads
and a tab-delimited file with five columns (`chrom`, `left_start`, `left_end`,
`right_start`, and `right_end`) which provide the 1-based inclusive start and end
positions of the primers for each amplicon.  The primer file must include headers, e.g:

```
chrom  left_start  left_end  right_start right_end
chr1   1010873     1010894   1011118     1011137
```

Paired end reads that map to a given amplicon position are trimmed so that the
alignment no-longer includes the primer sequences. All other aligned reads have the
_maximum primer length trimmed_!

Reads that are trimmed will have the `NM`, `UQ` and `MD` tags cleared as they are no longer
guaranteed to be accurate.  If a reference is provided the reads will be re-sorted
by coordinate after trimming and the `NM`, `UQ` and `MD` tags recalculated.

If the input BAM is not `queryname` sorted it will be sorted internally so that mate
information between paired-end reads can be corrected before writing the output file.

The `--first-of-pair` option will cause only the first of pair (R1) reads to be trimmed
based solely on the primer location of R1.  This is useful when there is a target
specific primer on the 5' end of R1 but no primer sequenced on R2 (eg. single gene-specific
primer target enrichment).  In this case, the location of each target specific primer should
be specified in an amplicons left or right primer exclusively.  The coordinates of the
non-specific-target primer should be `-1` for both start and end, e.g:

```
chrom  left_start  left_end  right_start right_end
chr1   1010873     1010894   -1          -1
chr2   -1          -1        1011118     1011137
```

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|primers|p|FilePath|File with primer locations.|Required|1||
|hard-clip|H|Boolean|If true, hard clip reads, else soft clip.|Optional|1|false|
|slop|S|Int|Match to primer locations +/- this many bases.|Optional|1|5|
|sort-order|s|SamOrder|Sort order of output BAM file (defaults to input sort order).|Optional|1||
|ref|r|PathToFasta|Optional reference fasta for recalculating NM, MD and UQ tags.|Optional|1||
|auto-trim-attributes|a|Boolean|Automatically trim extended attributes that are the same length as bases.|Optional|1|false|
|first-of-pair||Boolean|Trim only first of pair reads (R1s), otherwise both ends of a pair|Optional|1|false|

