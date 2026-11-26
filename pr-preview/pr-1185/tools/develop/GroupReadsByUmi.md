---
title: GroupReadsByUmi
---

# GroupReadsByUmi

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Groups reads together that appear to have come from the same original molecule. Reads
are grouped by template, and then templates are sorted by the 5' mapping positions of
the reads from the template, used from earliest mapping position to latest. Reads that
have the same end positions are then sub-grouped by UMI sequence.

Accepts reads in any order (including `unsorted`) and outputs reads sorted by:

   1. The lower genome coordinate of the two outer ends of the templates (strand-aware)
   2. Strand orientation
   3. The cellular barcode (CB tag)
   4. The assigned UMI tag
   5. Read Name
   6. The sequencing library

It is recommended to sort the reads into template-coordinate (i.e. `SO:unsorted GO:query SS:unsorted:template-coordinate`)
prior to running this tool to avoid this tool re-sorting the input.  It is recommended to use
`samtools sort --template-coordinate --threads $(nrpoc)` for the pre-sorting.
The output will always be written in template-coordinate order.

During grouping, reads and templates are filtered out as follows:

1. Templates are filtered if all reads for the template are unmapped
2. Templates are filtered if any non-secondary, non-supplementary read has mapping quality < `min-map-q`
3. Templates are filtered if any UMI sequence contains one or more `N` bases
4. Templates are filtered if `--min-umi-length` is specified and the UMI does not meet the length requirement
5. Records are filtered out if flagged as either secondary or supplementary

Grouping of UMIs is performed by one of four strategies:

1. **identity**:  only reads with identical UMI sequences are grouped together. This strategy
                  may be useful for evaluating data, but should generally be avoided as it will
                  generate multiple UMI groups per original molecule in the presence of errors.
2. **edit**:      reads are clustered into groups such that each read within a group has at least
                  one other read in the group with <= edits differences and there are inter-group
                  pairings with <= edits differences. Effective when there are small numbers of
                  reads per UMI, but breaks down at very high coverage of UMIs.
3. **adjacency**: a version of the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755)
                  that allows for errors between UMIs but only when there is a count gradient.
4. **paired**:    similar to adjacency but for methods that produce template such that a read with A-B is related
                  to but not identical to a read with B-A. Expects the UMI sequences to be stored in a single SAM
                  tag separated by a hyphen (e.g. `ACGT-CCGG`) and allows for one of the two UMIs to be absent
                  (e.g. `ACGT-` or `-ACGT`). The molecular IDs produced have more structure than for single
                  UMI strategies and are of the form `{base}/{A|B}`. E.g. two UMI pairs would be mapped as
                  follows AAAA-GGGG -> 1/A, GGGG-AAAA -> 1/B.

Strategies `edit`, `adjacency`, and `paired` make use of the `--edits` parameter to control the matching of
non-identical UMIs.

UMIs may have different lengths, but only UMIs with the same length are compared, so reads whose UMIs have
different lengths are never assigned to the same molecule.  For UMIs with multiple segments separated by `-`, each
segment must have the same length (e.g. `ACGT-CCGGA` is never grouped with `ACGTC-CGGA` or `ACGT-CCGG`).  For the
`paired` strategy, the two orientations of a UMI pair (e.g. `ACGT-CCGGA` and `CCGGA-ACGT`) are compared with
each other as usual.

If `--min-umi-length=len` is specified then reads that have fewer than `len` bases in their UMI are discarded.  The
number of bases is the number of A, C, G and T characters in the UMI: dashes and other characters are not counted,
and reads with an N in their UMI are discarded earlier.  The UMIs of the remaining reads are not changed, so UMIs
of different lengths are still never grouped together. This filter supports every strategy, including `paired`,
where bases are counted across both UMI segments.

To group UMIs of different lengths together, specify `--truncate` together with `--min-umi-length=len`. Each UMI
is truncated to the shortest prefix containing `len` bases, retaining dashes before the last kept base. Only UMIs
whose truncated segments have the same lengths are compared. For example, at six bases `ACGTAC` and `ACGTACGT`
both become `ACGTAC`, but `AC-GTAC` and `ACG-TAC` stay separate. A dash immediately after the last kept base is
dropped: at eight bases `ACGTACGT-` and `ACGTACGT-AAAA` both become `ACGTACGT`. Truncation can also remove later
segments entirely: at four bases `ACGT-ACGTAC` and `ACGTAC-GTAC` both become `ACGT`.

Truncation requires `--raw-tag RX` so that the working UMIs are available to consensus callers. RX is truncated on
both mates, and QX, when present, must match the original RX length and is truncated to the same character prefix
(including spaces corresponding to RX dashes). When a record's RX changes and OX is absent, the original RX is
saved in OX and the original QX, if present, in BZ (replacing any existing BZ). If OX already exists, both OX and BZ
are left untouched. Missing qualities are never synthesized. Only `--truncate` is rejected with the `paired`
strategy; `--min-umi-length` alone is supported. Without `--truncate`, a single warning is logged when
`--min-umi-length` retains UMIs with different segment lengths at the same position: `--min-umi-length` is a filter only.

If the `--mark-duplicates` option is given, reads will also have their duplicate flag set in the BAM file.
Each tag-family is treated separately, and a single template within the tag family is chosen to be the "unique"
template and marked as non-duplicate, while all other templates in the tag family are then marked as duplicate.
There are a few limitations of duplicate-marking mode (vs. e.g. Picard MarkDuplicates):

1. read pairs with one unmapped read are duplicate-marked independently from read pairs with both reads mapped
2. secondary and supplementary records are discarded

Note: the `--min-map-q` parameter defaults to 0 in duplicate marking mode and 1 otherwise, and is directly settable
on the command line.

Multi-threaded operation is supported via the `--threads/-@` option. This only applies to the Adjacency and Paired
strategies. Additionally the only operation that is multi-threaded is the comparisons of UMIs at the same genomic
position.  Running with e.g. `--threads 8` can provide a _substantial_ reduction in runtime when there are many
UMIs observed at the same genomic location, such as can occur in amplicon sequencing or ultra-deep coverage data.

## Single-Cell / Cell Barcode Support

When processing data with cell barcodes (e.g. single-cell or single-nuclei sequencing), the `--cell-tag` option
controls how cell identity is incorporated into grouping. Reads at the same genomic coordinates are partitioned
by cell barcode _before_ molecular index assignment, so reads from different cells will never be grouped together
even if they share a UMI and mapping position. No correction is performed on the cell barcode itself; it is
treated as a known/fixed value. The default tag is `CB`.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input BAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|The output BAM file.|Optional|1|/dev/stdout|
|family-size-histogram|f|FilePath|Optional output of tag family size counts.|Optional|1||
|grouping-metrics|g|FilePath|Optional output of UMI grouping metrics.|Optional|1||
|raw-tag|t|String|The tag containing the raw UMI.|Optional|1|RX|
|assign-tag|T|String|The output tag for UMI grouping.|Optional|1|MI|
|cell-tag|c|String|The tag containing the cell barcode.|Optional|1|CB|
|mark-duplicates|d|Boolean|Turn on duplicate marking mode.|Optional|1|false|
|min-map-q|m|Int|Minimum mapping quality for mapped reads.|Optional|1||
|include-non-pf-reads|n|Boolean|Include non-PF reads.|Optional|1|false|
|strategy|s|Strategy|The UMI assignment strategy.|Required|1||
|edits|e|Int|The allowable number of edits between UMIs.|Optional|1|1|
|min-umi-length|l|Int|The minimum UMI length, in bases.  Discard reads whose UMI has fewer bases (A, C, G, T) than this; dashes and other characters are not counted.  UMIs are not truncated unless `--truncate` is also given, so UMIs of different lengths are still never grouped together. Supported with every strategy; paired UMIs are filtered by the total base count across both segments.|Optional|1||
|truncate||Boolean|Truncate RX to the shortest prefix containing `--min-umi-length` bases, keeping dashes before the last kept base. Truncate matching QX qualities too, preserving originals in OX/BZ when OX is absent. Requires `--min-umi-length` and `--raw-tag RX`; not supported with the paired strategy.|Optional|1|false|
|threads|@|Int|Number of threads to use when comparing UMIs. Only recommended for amplicon or similar data.|Optional|1|1|

