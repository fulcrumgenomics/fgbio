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

   1. The lower genome coordinate of the two outer ends of the templates
   2. The sequencing library
   3. The assigned UMI tag
   4. Read Name

If the input is not template-coordinate sorted (i.e. `SO:unsorted GO:query SS:unsorted:template-coordinate`), then
this tool will re-sort the input. The output will be written in template-coordinate order.

During grouping, reads and templates are filtered out as follows:

1. Templates are filtered if all reads for the template are unmapped
2. Templates are filtered if any non-secondary, non-supplementary read has mapping quality < `min-map-q`
3. Templates are filtered if R1 and R2 are mapped to different chromosomes and `--allow-inter-contig` is false
4. Templates are filtered if any UMI sequence contains one or more `N` bases
5. Templates are filtered if `--min-umi-length` is specified and the UMI does not meet the length requirement
6. Reads are filtered out if flagged as secondary and `--include-secondary` is false
7. Reads are filtered out if flagged as supplementary and `--include-supplementary` is false

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

By default, all UMIs must be the same length. If `--min-umi-length=len` is specified then reads that have a UMI
shorter than `len` will be discarded, and when comparing UMIs of different lengths, the first len bases will be
compared, where `len` is the length of the shortest UMI. The UMI length is the number of [ACGT] bases in the UMI
(i.e. does not count dashes and other non-ACGT characters). This option is not implemented for reads with UMI pairs
(i.e. using the paired assigner).

If the `--mark-duplicates` option is given, reads will also have their duplicate flag set in the BAM file.
Each tag-family is treated separately, and a single template within the tag family is chosen to be the "unique"
template and marked as non-duplicate, while all other templates in the tag family are then marked as duplicate.
One limitation of duplicate-marking mode, vs. e.g. Picard MarkDuplicates, is that read pairs with one unmapped read
are duplicate-marked independently from read pairs with both reads mapped.

Several parameters have different defaults depending on whether duplicates are being marked or not (all are
directly settable on the command line):

  1. `--min-map-q` defaults to 0 in duplicate marking mode and 1 otherwise
  2. `--include-secondary` defaults to true in duplicate marking mode and false otherwise
  3. `--include-supplementary` defaults to true in duplicate marking mode and false otherwise

Multi-threaded operation is supported via the `--threads/-@` option. This only applies to the Adjacency and Paired
strategies. Additionally the only operation that is multi-threaded is the comparisons of UMIs at the same genomic
position.  Running with e.g. `--threads 8` can provide a _substantial_ reduction in runtime when there are many
UMIs observed at the same genomic location, such as can occur in amplicon sequencing or ultra-deep coverage data.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input BAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|The output BAM file.|Optional|1|/dev/stdout|
|family-size-histogram|f|FilePath|Optional output of tag family size counts.|Optional|1||
|raw-tag|t|String|The tag containing the raw UMI.|Optional|1|RX|
|assign-tag|T|String|The output tag for UMI grouping.|Optional|1|MI|
|mark-duplicates|d|Boolean|Turn on duplicate marking mode.|Optional|1|false|
|include-secondary|S|Boolean|Include secondary reads.|Optional|1||
|include-supplementary|U|Boolean|Include supplementary reads.|Optional|1||
|min-map-q|m|Int|Minimum mapping quality for mapped reads.|Optional|1||
|include-non-pf-reads|n|Boolean|Include non-PF reads.|Optional|1|false|
|strategy|s|Strategy|The UMI assignment strategy.|Required|1||
|edits|e|Int|The allowable number of edits between UMIs.|Optional|1|1|
|min-umi-length|l|Int|The minimum UMI length. If not specified then all UMIs must have the same length, otherwise discard reads with UMIs shorter than this length and allow for differing UMI lengths.|Optional|1||
|allow-inter-contig|x|Boolean|DEPRECATED: this option will be removed in future versions and inter-contig reads will be automatically processed.|Optional|1|true|
|threads|@|Int|Number of threads to use when comparing UMIs. Only recommended for amplicon or similar data.|Optional|1|1|

