---
title: StripFastqReadNumbers
---

# StripFastqReadNumbers

## Overview
**Group:** Personal

Removes trailing /# from read names in fastq.  Read names that end in a slash and
a single digit will have the slash and digit removed.  Read names that do not end
with a slash and a number will be untouched.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToFastq|Input fastq file.|Required|1||
|output|o|PathToFastq|Output fastq file.|Required|1||

