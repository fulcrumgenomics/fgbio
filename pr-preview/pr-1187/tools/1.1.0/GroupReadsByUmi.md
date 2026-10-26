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

Reads are aggressively filtered out so that only high quality reads/mappings are taken forward. Single-end
reads must have mapping quality >= `min-map-q`.  Paired-end reads must have both reads mapped to the same
chromosome with both reads having mapping quality >= `min-mapq`.  (Note: the `MQ` tag is required on reads
with mapped mates).

This is done with the expectation that the next step is building consensus reads, where
it is undesirable to either:

   1. Assign reads together that are really from different source molecules
   2. Build two groups from reads that are really from the same molecule

Errors in mapping reads could lead to both and therefore are minimized.

Grouping of UMIs is performed by one of three strategies:

1. **identity**:  only reads with identical UMI sequences are grouped together. This strategy
                  may be useful for evaluating data, but should generally be avoided as it will
                  generate multiple UMI groups per original molecule in the presence of errors.
2. **edit**:      reads are clustered into groups such that each read within a group has at least
                  one other read in the group with <= edits differences and there are inter-group
                  pairings with <= edits differences. Effective when there are small numbers of
                  reads per UMI, but breaks down at very high coverage of UMIs.
3. **adjacency**: a version of the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755)
                  that allows for errors between UMIs but only when there is a count gradient.
4. **paired**:    similar to adjacency but for methods that produce template with a pair of UMIs
                  such that a read with A-B is related to but not identical to a read with B-A.
                  Expects the pair of UMIs to be stored in a single tag, separated by a hyphen
                  (e.g. `ACGT-CCGG`).  The molecular IDs produced have more structure than for single
                  UMI strategies, and are of the form `{base}/{AB|BA}`. E.g. two UMI pairs would be
                  mapped as follows AAAA-GGGG -> 1/AB, GGGG-AAAA -> 1/BA.

`edit`, `adjacency` and `paired` make use of the `--edits` parameter to control the matching of
non-identical UMIs.

By default, all UMIs must be the same length. If `--min-umi-length=len` is specified then reads that have a UMI
shorter than `len` will be discarded, and when comparing UMIs of different lengths, the first len bases will be
compared, where `len` is the length of the shortest UMI. The UMI length is the number of [ACGT] bases in the UMI
(i.e. does not count dashes and other non-ACGT characters). This option is not implemented for reads with UMI pairs
(i.e. using the paired assigner).

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|The input BAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|The output BAM file.|Optional|1|/dev/stdout|
|family-size-histogram|f|FilePath|Optional output of tag family size counts.|Optional|1||
|raw-tag|t|String|The tag containing the raw UMI.|Optional|1|RX|
|assign-tag|T|String|The output tag for UMI grouping.|Optional|1|MI|
|min-map-q|m|Int|Minimum mapping quality.|Optional|1|30|
|include-non-pf-reads|n|Boolean|Include non-PF reads.|Optional|1|false|
|strategy|s|Strategy|The UMI assignment strategy.|Required|1||
|edits|e|Int|The allowable number of edits between UMIs.|Optional|1|1|
|min-umi-length|l|Int|The minimum UMI length. If not specified then all UMIs must have the same length, otherwise discard reads with UMIs shorter than this length and allow for differing UMI lengths.|Optional|1||

