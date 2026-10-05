---
title: UpdateGffContigNames
---

# UpdateGffContigNames

## Overview
**Group:** Utilities

Updates the contig names in a GFF.

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.

Please note: the output GFF will be in the same order as the input GFF.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|FilePath|Input GFF.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|output|o|FilePath|Output GFF.|Required|1||
|skip-missing||Boolean|Skip contigs in the GFF that are not found in the sequence dictionary.|Optional|1|false|

