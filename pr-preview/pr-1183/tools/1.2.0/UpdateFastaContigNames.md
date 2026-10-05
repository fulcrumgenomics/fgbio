---
title: UpdateFastaContigNames
---

# UpdateFastaContigNames

## Overview
**Group:** FASTA

Updates the sequence names in a FASTA.

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToFasta|Input FASTA.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|output|o|PathToFasta|Output FASTA.|Required|1||
|line-length|l|Int|Line length or sequence lines.|Optional|1|100|
|skip-missing||Boolean|Skip missing source contigs.|Optional|1|false|

