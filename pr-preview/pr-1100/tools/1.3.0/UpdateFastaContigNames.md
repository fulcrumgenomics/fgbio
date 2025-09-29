---
title: UpdateFastaContigNames
---

# UpdateFastaContigNames

## Overview
**Group:** FASTA

Updates the sequence names in a FASTA.

The name of each sequence must match one of the names (including aliases) in the given sequence dictionary.  The
new name will be the primary (non-alias) name in the sequence dictionary.

By default, the sort order of the contigs will be the same as the input FASTA.  Use the `--sort-by-dict` option to
sort by the input sequence dictionary.  Furthermore, the sequence dictionary may contain **more** contigs than the
input FASTA, and they wont be used.

Use the `--skip-missing` option to skip contigs in the input FASTA that cannot be renamed (i.e. who are not present
in the input sequence dictionary); missing contigs will not be written to the output FASTA.  Finally, use the
`--default-contigs` option to specify an additional FASTA which will be queried to locate contigs not present in
the input FASTA but present in the sequence dictionary.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToFasta|Input FASTA.|Required|1||
|dict|d|PathToSequenceDictionary|The path to the sequence dictionary with contig aliases.|Required|1||
|output|o|PathToFasta|Output FASTA.|Required|1||
|line-length|l|Int|Line length or sequence lines.|Optional|1|100|
|skip-missing||Boolean|Skip missing source contigs (will not be outputted).|Optional|1|false|
|sort-by-dict||Boolean|Sort the contigs based on the input sequence dictionary.|Optional|1|false|
|default-contigs||PathToFasta|Add sequences from this FASTA when contigs in the sequence dictionary are missing from the input FASTA.|Optional|1||

