---
title: UpdateDelimitedFileContigNames
---

# UpdateDelimitedFileContigNames

## Overview
**Group:** Utilities

Updates the contig names in columns of a delimited data file (e.g. CSV, TSV).

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.  Use `--skip-missing` to ignore lines
where a contig name could not be updated (i.e. missing from the sequence dictionary).

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|FilePath|Input delimited data file.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|columns|c|Int|The column indices for the contig names (0-based).|Required|Unlimited||
|delimiter|T|Char|The delimiter|Optional|1|	|
|comment|H|String|Treat lines with this starting string as comments (always printed)|Optional|1|#|
|output|o|FilePath|Output delimited data file.|Required|1||
|output-first-num-lines|n|Int|Output the first `N` lines as-is (always printed).|Optional|1|0|
|skip-missing||Boolean|Skip lines where a contig name could not be updated (i.e. missing from the sequence dictionary).|Optional|1|false|

