---
title: FindSwitchbackReads
---

# FindSwitchbackReads

## Overview
**Group:** SAM/BAM

Finds reads where a template switch occurred during library construction.  Some library construction methods,
notably ultra-low-input shotgun methods, are prone to template switching events that create molecules
(templates, inserts) that instead of being a linear copy of a segment of genomic DNA, instead are chimeras formed
by starting on one strand of DNA and then, switching to the opposite strand.  Frequently when the switch occurs
there may be a small offset between where the first strand was departed and the opposite strand joined.

## Algorithm

Templates that contain strand switch events (switch-backs) are found by this tool in two different ways:

1. By looking at reads that contain soft-clipped bases at their 5' end that, when reverse complemented, matches the
   genome proximal to the 5'-most mapped base of the read.  We call these matches
   "read based switchbacks".  Finding read based switchbacks is based on several parameters:

   1. `max-offset` controls how far away to search for the reverse-complemented sequence.  The default value of
      `35` allows matches to be found when the soft-clipped sequence matches the genome _starting_ at most 35bp
      from the 5' mapped position of the read, and reading in the opposite direction.
   2. `min-length` controls the minimum number of soft-clipped bases that must exist to trigger the search.
      Given that the search looks at `2 * max-offset` locations (default=70) it is important that `min-length`
      is set such that `4^min-length >> 2 * `max-offset` in order to avoid false positives.
   3. `max-error-rate` allows for some mismatches to exist between the soft-clipped sequence and the genome when matching.

2. By identifying templates with `FF` or `RR` (aka tandem) orientation where it is surmised that the template
   switch occurred in the un-sequenced region of the template between R1 and R2.  We call these `tandem based
   switchbacks`.  This is controlled by a single parameter, `max-gap`, which causes the tool to only identify a
   tandem read pair as a switch-back _if_ the gap between the end of the first read and the start of the second
   read is `+/- max-gap`.

By default, when a switch-back template is identified, the primary reads are made unmapped (and the original
alignment stored in the OA tag) and all secondary and supplementary alignments are discarded.  This can be
disabled with the `--dont-unmap` or `-d` option.

All reads from a switch-back template are also tagged with an `sb` tag that describes the nature of the
switchback.  If the template was identified base on soft-clipped sequence within a read the format is:

```sb:Z:r,[read|mate],{offset},{length}```

If the template is identified due to it's tandem pair orientation then the format is:

```sb:Z:t,{gap}```

## Inputs and Outputs

The tool takes as input a SAM or BAM file, and by default consumes from `stdin`.  The primary output is also a
SAM or BAM file, and defaults to compressed BAM on `stdout`.  This allows the tool to be run immediately after
an aligner in a pipe, e.g. `bwa mem ref.fa r1.fq r2.fq | fgbio -Xmx8g --ref=ref.fa | ...`.

If the input is neither `queryname` sorted nor `queryname` grouped (i.e. all reads with the same name grouped
together) it will be sorted into `queryname` order by the tool before processing.

By default the output BAM is produced in the order the reads were processed (i.e. the input ordering _or_
queryname sorted if sorting was required).  This can be overridden with the `--sort-order` option.

A number of text files are also produced if the `--metrics` option is specified.  E.g. when specifying
`--metrics=s1.switchback` the following files are produced:

1. `s1.switchback.summary.txt`: A table of summary metrics describing the number of reads, switchbacks, etc.
2. `s1.switchback.lengths.txt`: A table of the distribution of observed switchback lengths in read-based switchbacks.
3. `s1.switchback.offsets.txt`: A table of the distribution of observed offsets in read-based switchbacks.
4. `s1.switchback.gaps.txt`: A table of the distribution of gap lengths in tampl
5. `s1.switchback.plots.pdf`: A PDF containing plots of the distributions from 2-4.

Note: because this tool accesses the reference genome in a random manner it pre-loads the entire reference fasta
into memory.  As a result the tool is best run with `-Xmx8g` to give it sufficient memory.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|Output BAM file.|Optional|1|/dev/stdout|
|metrics|m|PathPrefix|Metrics output file.|Optional|1||
|sort-order|s|SamOrder|Output sort order.|Optional|1||
|ref|r|PathToFasta|Reference genome fasta file.|Required|1||
|max-offset|O|Int|Maximum offset between end the two segments of the read on the reference. Set to 0 to disable read-based checks.|Optional|1|35|
|max-gap|g|Int|Maximum gap between R1 and R2 of tandem reads to call a template a switchback. Set to 0 to disable tandem-based checks.|Optional|1|500|
|min-length|l|Int|Minimum match length of the switched back segment.|Optional|1|6|
|max-error-rate|e|Double|Maximum mismatch error rate of switchback match to genome.|Optional|1|0.1|
|dont-unmap|d|Boolean|IF true, do NOT unmap reads from switchback templates.|Optional|1|false|

