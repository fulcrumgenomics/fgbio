---
title: HardMaskFasta
---

# HardMaskFasta

## Overview
**Group:** FASTA

Converts soft-masked sequence to hard-masked in a FASTA file. All lower case bases are
converted to Ns, all other bases are left unchanged.  Line lengths are also standardized
to allow easy indexing with `samtools faidx`"

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToFasta|Input FASTA file.|Required|1||
|output|o|PathToFasta|Output FASTA file.|Required|1||
|line-length|l|Int|Line length or sequence lines.|Optional|1|100|

