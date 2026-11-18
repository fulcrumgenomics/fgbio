---
title: SortFastq
---

# SortFastq

## Overview
**Group:** FASTQ

Sorts a FASTQ file.  Sorts the records in a FASTQ file based on the lexicographic ordering
of their read names.  Input and output files can be either uncompressed or gzip-compressed.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToFastq|Input fastq file.|Required|1||
|output|o|PathToFastq|Output fastq file.|Required|1||
|max-records-in-ram|m|Int|Maximum records to keep in RAM at one time.|Optional|1|500000|

