---
title: ClipBam
---

# ClipBam

## Overview
**Group:** SAM/BAM

Clips reads from the same template. Ensures that at least N bases are clipped from any end of the read (i.e.
R1 5' end, R1 3' end, R2 5' end, and R2 3' end).  Optionally clips reads from the same template to eliminate overlap
between the reads.  This ensures that downstream processes, particularly variant calling, cannot double-count
evidence from the same template when both reads span a variant site in the same template.

Clipping overlapping reads is only performed on `FR` read pairs, and is implemented by clipping approximately half
the overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be substituted
using the `--soft-clip` parameter.

Secondary alignments and supplemental alignments are not clipped, but are passed through into the
output.

In order to correctly clip reads by template and update mate information, the input BAM must be either
`queryname` sorted or `query` grouped.  If your input BAM is not in an appropriate order the sort can be
done in streaming fashion with, for example:

```
samtools sort -n -u in.bam | fgbio ClipBam -i /dev/stdin ...
```

The output sort order may be specified with `--sort-order`.  If not given, then the output will be in the same
order as input.

Any existing `NM`, `UQ` and `MD` tags are repaired, and mate-pair information updated.

Three clipping modes are supported:
1. `Soft` - soft-clip the bases and qualities.
2. `SoftWithMask` - soft-clip and mask the bases and qualities (make bases Ns and qualities the minimum).
3. `Hard` - hard-clip the bases and qualities.

The `--upgrade-clipping` parameter will convert all existing clipping in the input to the given more stringent mode:
from `Soft` to either `SoftWithMask` or `Hard`, and `SoftWithMask` to `Hard`. In all other cases, clipping remains
the same prior to applying any other clipping criteria.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input SAM or BAM file of aligned reads in coordinate order.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Required|1||
|metrics|m|FilePath|Optional output of clipping metrics.|Optional|1||
|ref|r|PathToFasta|Reference sequence fasta file.|Required|1||
|clipping-mode|c|ClippingMode|The type of clipping to perform.|Optional|1|Hard|
|auto-clip-attributes|a|Boolean|Automatically clip extended attributes that are the same length as bases.|Optional|1|false|
|upgrade-clipping|H|Boolean|Upgrade all existing clipping in the input to the given clipping mode prior to applying any other clipping criteria.|Optional|1|false|
|read-one-five-prime||Int|Require at least this number of bases to be clipped on the 5' end of R1|Optional|1|0|
|read-one-three-prime||Int|Require at least this number of bases to be clipped on the 3' end of R1|Optional|1|0|
|read-two-five-prime||Int|Require at least this number of bases to be clipped on the 5' end of R2|Optional|1|0|
|read-two-three-prime||Int|Require at least this number of bases to be clipped on the 3' end of R2|Optional|1|0|
|clip-overlapping-reads||Boolean|Clip overlapping reads.|Optional|1|false|
|clip-bases-past-mate||Boolean|Clip reads in FR pairs that sequence past the far end of their mate.|Optional|1|false|
|sort-order|S|SamOrder|The sort order of the output. If not given, output will be in the same order as input if the input.|Optional|1||
