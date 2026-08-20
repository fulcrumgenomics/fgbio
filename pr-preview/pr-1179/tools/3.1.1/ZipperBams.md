---
title: ZipperBams
---

# ZipperBams

## Overview
**Group:** SAM/BAM

Zips together an unmapped and mapped BAM to transfer metadata into the output BAM.

Both the unmapped and mapped BAMs _must_ be a) queryname sorted or grouped (i.e. all records with the same
name are grouped together in the file), and b) have the same ordering of querynames.  If either of these are
violated the output is undefined!

All tags present on the unmapped reads are transferred to the mapped reads.  The options `--tags-to-reverse`
and `--tags-to-revcomp` will cause tags on the unmapped reads to be reversed or reverse complemented before
being copied to reads mapped to the negative strand.  These options can take a mixture of two-letter tag names
and the names of tag sets, which will be expanded into sets of tag names.  Currently the only named tag set
is "Consensus" which contains all the per-base consensus tags produced by fgbio consensus callers.

By default the mapped BAM is read from standard input (stdin) and the output BAM is written to standard
output (stdout). This can be changed using the `--input/-i` and `--output/-o` options.

By default the output BAM file is emitted in the same order as the input BAMs.  This can be overridden
using the `--sort` option, though in practice it may be faster to do the following:

```
fgbio --compression 0 ZipperBams -i mapped.bam -u unmapped.bam -r ref.fa | samtools sort -@ $(nproc)
```

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Mapped SAM or BAM.|Optional|1|/dev/stdin|
|unmapped|u|PathToBam|Unmapped SAM or BAM.|Required|1||
|ref|r|PathToFasta|Path to the reference used in alignment. Must have accompanying .dict file.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Optional|1|/dev/stdout|
|tags-to-remove||String|Tags to remove from the mapped BAM records.|Optional|Unlimited||
|tags-to-reverse||String|Set of optional tags to reverse on reads mapped to the negative strand.|Optional|Unlimited||
|tags-to-revcomp||String|Set of optional tags to reverse complement on reads mapped to the negative strand.|Optional|Unlimited||
|sort|s|SamOrder|Sort the output BAM into the given order.|Optional|1||
|buffer|b|Int|Buffer this many read-pairs while reading the input BAMs.|Optional|1|5000|

