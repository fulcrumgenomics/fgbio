---
title: GenerateRegionsFromFasta
---

# GenerateRegionsFromFasta

## Overview
**Group:** Personal

Generates a list of `freebayes`/`bamtools` region specifiers.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input||PathToFasta|The input FASTA.|Required|1||
|output||Path|The output.|Optional|1|/dev/stdout|
|region-size||Int|The size of the regions to output.|Optional|1|100000|

