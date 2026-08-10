---
title: UpdateVcfContigNames
---

# UpdateVcfContigNames

## Overview
**Group:** VCF/BCF

Updates then contig names in a VCF.

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToVcf|Input VCF.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|output|o|PathToVcf|Output VCF.|Required|1||
|skip-missing||Boolean|Skip missing contigs.|Optional|1|false|

