---
title: AssignPrimers
---

# AssignPrimers

## Overview
**Group:** SAM/BAM

Assigns reads to primers post-alignment. Takes in a BAM file of aligned reads and a tab-delimited file with five columns
(`chrom`, `left_start`, `left_end`, `right_start`, and `right_end`) which provide the 1-based inclusive start and
end positions of the primers for each amplicon.  The primer file must include headers, e.g:

```
chrom  left_start  left_end  right_start right_end
chr1   1010873     1010894   1011118     1011137
```

Optionally, a sixth column column `id` may be given with a unique name for the amplicon.  If not given, the
coordinates of the amplicon's primers will be used:
  `<chrom>:<left_start>-<left_end>,<chrom>:<right_start>:<right_end>`

Each read is assigned independently of its mate (for paired end reads). The primer for a read is assumed to be
located at the start of the read in 5' sequencing order.  Therefore, a positive strand
read will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
strand read will use its aligned end position to match against the amplicon's right-most coordinate.

For paired end reads, the assignment for mate will also be stored in the current read, using the same procedure as
above but using the mate's coordinates.  This requires the input BAM have the mate-cigar ("MC") SAM tag.  Read
pairs must have both ends mapped in forward/reverse configuration to have an assignment.  Furthermore, the amplicon
assignment may be different for a read and its mate.  This may occur, for example, if tiling nearby amplicons and
a large deletion occurs over a given primer and therefore "skipping" an amplicon.  This may also occur if there are
translocations across amplicons.

The output will have the following tags added:
- ap: the assigned primer coordinates (ex. `chr1:1010873-1010894`)
- am: the mate's assigned primer coordinates (ex. `chr1:1011118-1011137`)
- ip: the assigned amplicon id
- im: the mate's assigned amplicon id (or `=` if the same as the assigned amplicon)

The read sequence of the primer is not checked against the expected reference sequence at the primer's genomic
coordinates.

In some cases, large deletions within one end of a read pair may cause a primary and supplementary alignments to be
produced by the aligner, with the supplementary alignment containing the primer end of the read (5' sequencing order).
In this case, the primer may not be assigned for this end of the read pair.  Therefore, it is recommended to prefer
or choose the primary alignment that has the closest aligned read base to the 5' end of the read in sequencing order.
For example, from `bwa` version `0.7.16` onwards, the `-5` option may be used.  Consider also using the `-q` option 
for `bwa` `0.7.16` as well, which is standard in `0.7.17` onwards when the `-5` option is used.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|metrics|m|FilePath|Output metrics file.|Required|1||
|primers|p|FilePath|File with primer locations.|Required|1||
|slop|S|Int|Match to primer locations +/- this many bases.|Optional|1|5|
|unclipped-coordinates|U|Boolean|True to based on the unclipped coordinates (adjust based on hard/soft clipping), otherwise the aligned bases|Optional|1|true|
|primer-coordinates-tag||String|The SAM tag for the assigned primer coordinate.|Optional|1|rp|
|mate-primer-coordinates-tag||String|The SAM tag for the mate's assigned primer coordinate.|Optional|1|mp|
|amplicon-identifier-tag||String|The SAM tag for the assigned amplicon identifier.|Optional|1|ra|
|mate-amplicon-identifier-tag||String|The SAM tag for the mate's assigned amplicon identifier.|Optional|1|ma|

