---
title: RandomizeBam
---

# RandomizeBam

## Overview
**Group:** SAM/BAM

Randomizes the order of reads in a SAM or BAM file. Randomization is done by sorting
on a hash of the `queryname` (and bases and quals if not query-grouping). By default
reads with the same query name are grouped together in the output file; this can be
turned off by specifying --query-group=false.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input SAM or BAM file.|Optional|1|/dev/stdin|
|output|o|PathToBam|The output SAM or BAM file.|Optional|1|/dev/stdout|
|seed|s|Int|Random seed.|Optional|1|42|
|query-group|q|Boolean|Group together reads by queryname.|Optional|1|true|

