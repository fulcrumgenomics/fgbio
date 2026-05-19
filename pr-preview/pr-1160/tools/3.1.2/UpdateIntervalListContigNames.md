---
title: UpdateIntervalListContigNames
---

# UpdateIntervalListContigNames

## Overview
**Group:** FASTA

Updates the sequence names in an Interval List file.

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.

Use `--skip-missing` to ignore intervals where a contig name could not be updated (i.e. missing from the sequence dictionary).

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToIntervals|Input interval list.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|output|o|PathToIntervals|Output interval list.|Required|1||
|skip-missing||Boolean|Skip contigs in the interval list that are not found in the sequence dictionary.|Optional|1|false|

