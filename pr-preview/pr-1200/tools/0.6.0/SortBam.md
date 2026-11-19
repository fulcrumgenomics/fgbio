---
title: SortBam
---

# SortBam

## Overview
**Group:** SAM/BAM

Sorts a SAM or BAM file.  Several sort orders are available:

1. **Coordinate**: sorts reads by their reference sequence and left-most aligned coordinate
2. **Queryname**: sort the reads by their query (i.e. read) name
3. **Random**: sorts the reads into a random order. The output is deterministic for any given input.
and several
4. **RandomQuery**: sorts the reads into a random order but keeps reads with the same
   queryname together. The ordering is deterministic for any given input.

Uses a temporary directory to buffer sets of sorted reads to disk. The number of reads kept in memory
affects memory use and can be changed with the `--max-records-in-ram` option.  The temporary directory
to use can be set with the fgbio global option `--tmp-dir`.

An example invocation might look like:

```bash
java -Xmx4g -jar fgbio.jar --tmp-dir=/my/big/scratch/volume \
  SortBam --input=queryname.bam --sort-order=Coordinate --output coordinate.bam
```

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM.|Required|1||
|output|o|PathToBam|Output SAM or BAM.|Required|1||
|sort-order|s|SamOrder|Order into which to sort the records.|Optional|1|Coordinate|
|max-records-in-ram|m|Int|Max records in RAM.|Optional|1|1000000|

